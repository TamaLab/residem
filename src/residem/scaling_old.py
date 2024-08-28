import mmtbx.f_model
import math
from cctbx import miller, maptbx
import numpy as np
from residem.master import DataInput
from cctbx.array_family import flex
from mmtbx import map_tools
import os
from collections import defaultdict
from residem.multiprocess_func import make_parallel
import tqdm
from iotbx.map_model_manager import map_model_manager
import pandas as pd
import matplotlib.pyplot as plt
from cctbx.xray import structure
import iotbx.map_tools
from scipy.signal import find_peaks
from scipy.optimize import minimize, least_squares
from iotbx.data_manager import DataManager
import pathlib
from functools import reduce
import operator
from scipy.optimize import curve_fit, minimize_scalar
from sklearn.metrics import r2_score
import matplotlib.backends.backend_pdf
import warnings
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

matplotlib.use('Agg')
font = {'weight': 'bold', 'size': 8}

matplotlib.rc('font', **font)


def deg_convert(theta):
    """ convert negative thetha to positive radians and converts as theta to radians
    """
    theta = math.fmod(math.fmod(theta, 360) + 360, 360) if theta < 0 else math.radians(theta)
    return theta


class ScaleModel(DataInput):
    def __init__(self, reference_pdb, reference_mtz, Triggered_mtz, reference_label=None, triggered_label=None,
                 high_resolution=None,
                 low_resolution=None, sigma=3.0, grid_step=0.5, multiscaling=False, occupancy=False,
                 weight="default", write=False, path=os.getcwd(), inverse_scaling=False, scale="aniso"):
        super().__init__(reference_pdb, reference_mtz, Triggered_mtz, reference_label, triggered_label, high_resolution,
                         low_resolution)
        self.scale = scale
        self.sigma = sigma
        self.write = write
        self.inverse_scaling = inverse_scaling
        ### input
        self.path = path
        self.model = mmtbx.model.manager(model_input=self.pdb_in)
        self.model.setup_scattering_dictionaries(
            scattering_table="n_gaussian",
            d_min=self.high_resolution)
        self.grid_step = grid_step
        self.resolution_factor = self.grid_step / self.high_resolution

        self.F_dark_miller = self.dark_dict[self.map_obs_label()[0]]
        self.F_light_miller = self.light_dict[self.map_obs_label()[1]]
        self.F_cal = self.dark_dict[self.map_obs_label()[2]]
        ### check crystal symmetry
        self.assertion = self.check_sym()

        info0 = self.F_dark_miller.info()
        info1 = self.F_light_miller.info()

        self.F_dark_miller = self.F_dark_miller.resolution_filter(d_min=self.high_resolution,
                                                                  d_max=self.low_resolution).set_info(info0)
        self.F_light_miller = self.F_light_miller.resolution_filter(d_min=self.high_resolution,
                                                                    d_max=self.low_resolution).set_info(info1)
        
        self.F_light_miller = self.F_light_miller.customized_copy(crystal_symmetry=self.F_dark_miller.crystal_symmetry())

        sel_dark = self.F_dark_miller.data() > self.F_dark_miller.sigmas() * sigma
        self.F_dark_cut = self.F_dark_miller.select(sel_dark).set_info(info0)
        sel_light = self.F_light_miller.data() > self.F_light_miller.sigmas() * sigma
        self.F_light_cut = self.F_light_miller.select(sel_light).set_info(info0)

        self.xrs = self.pdb_in.xray_structure_simple()
        self.light_xrs = structure(scatterers=self.xrs.scatterers(),
                                   crystal_symmetry=self.F_light_miller.crystal_symmetry())

        self.F_dark_common = self.F_dark_cut.common_set(self.F_light_cut,
                                                        assert_is_similar_symmetry=self.assertion).average_bijvoet_mates()
        self.F_light_common = self.F_light_cut.common_set(self.F_dark_cut,
                                                          assert_is_similar_symmetry=self.assertion).average_bijvoet_mates()

        ## f_model
        self.f_model_dark = mmtbx.f_model.manager(f_obs=self.F_dark_common, xray_structure=self.xrs)
        self.f_model_dark.update_all_scales()
        self.f_model_light = mmtbx.f_model.manager(f_obs=self.F_light_common, xray_structure=self.light_xrs)
        self.f_model_light.update_all_scales()
        self.f_model_dark_s = mmtbx.f_model.manager(
            f_obs=self.f_model_dark.f_obs().common_set(self.f_model_light.f_obs(),
                                                       assert_is_similar_symmetry=self.assertion),
            xray_structure=self.xrs)
        self.f_model_light_s = mmtbx.f_model.manager(
            f_obs=self.f_model_light.f_obs().common_set(self.f_model_dark.f_obs(),
                                                        assert_is_similar_symmetry=self.assertion),
            xray_structure=self.light_xrs)
        self.F_cal_common = self.F_cal.common_set(self.f_model_dark_s.f_obs(),
                                                  assert_is_similar_symmetry=self.assertion)

        ## cell symmetry
        self.indices = self.f_model_dark_s.f_obs().indices()
        self.cs = self.F_dark_common.crystal_symmetry()
        self.uc = self.cs.unit_cell()
        self.miller_set = miller.set(self.cs, self.indices, anomalous_flag=False)

        ### model map scaling
        self.f_model_dark_scaled = self.map_model_scaling(self.f_model_dark_s)
        self.f_model_light_scaled = self.map_model_scaling(self.f_model_light_s)
        self.F_dark = self.f_model_dark_scaled.data().as_numpy_array()
        self.SIGF_dark = self.f_model_dark_scaled.sigmas().as_numpy_array()
        self.F_light = self.f_model_light_scaled.data().as_numpy_array()
        self.SIGF_light = self.f_model_light_scaled.sigmas().as_numpy_array()
        self.F_calc = self.F_cal_common.data().as_numpy_array()
        ## weight and occupancy
        self.weight = None if weight is None else weight.lower() if weight.lower() in ["ren", "ursby"] else None
        # self.weight = "default" if weight is None else weight.lower() in ["ren", "ursby"]
        self.occupancy = occupancy

        # scaling
        if scale == "iso":
            self.f_model_light_scaled_obs = self.linear_scaling(self.f_model_light_scaled, self.f_model_dark_scaled)
        elif scale == "multi":
            F_dark_multiscaled = self.map_model_multiscaling(self.f_model_light_scaled_obs, self.f_model_dark_scaled)
            self.diff_map_obs = self.diff_map(self.f_model_dark_s, self.f_model_light_scaled_obs,
                                              F_dark_multiscaled)

        else:
            # print(self.indices)
            # self.f_model_light_scaled_obs = self.linear_scaling(self.f_model_light_scaled, self.f_model_dark_scaled)
            self.f_model_light_scaled_obs = self.aniso_scaling(self.f_model_light_scaled, self.f_model_dark_scaled)

        ### differnce map
        if scale != "multi":
            self.delta_F = self.f_model_light_scaled_obs.data() - self.f_model_dark_scaled.data()
            self.delta_F_model = self.F_light - self.F_calc
        else:
            self.delta_F = self.diff_map_obs.data()
            self.delta_F_model = self.F_light - self.F_calc

        ## weight calculation
        ## K weight
        self.weight_k = self.W_weight()
        self.weight_k = self.W_weight()
        self.q_weight = self.compute_q_weight().data()
        self.wDF = flex.double(self.weight_k * self.delta_F)
        self.rDF = flex.double(self.q_weight.as_numpy_array() * self.delta_F)
        self.mtz_dataset_object = self.mtz_dataset(f"%s/Difference_map_weighted_all.mtz" % self.path)
        self.mtz_extrapolated(self.delta_F, "Extrapolated_map_all_default.mtz")
        ## map coefficient
        self.FOFODF = self.map_coefficient(self.delta_F, coefficient="FO", write=self.write,
                                           name=f"%s/F_obs_minus_F_obs.ccp4" % self.path, map_type="diff")
        self.FOFODF_K = self.map_coefficient(self.wDF, coefficient="FO", write=self.write,
                                             name=f"%s/F_obs_minus_F_obs_ren_weight.ccp4" % self.path, map_type="diff")
        self.FOFODF_Q = self.map_coefficient(self.rDF, write=self.write, coefficient="FO",
                                             name=f"%s/F_obs_minus_F_obs_ursby_weight.ccp4" % self.path,
                                             map_type="diff")

        # ## extrapolated map
        self.F_cal_extra = defaultdict(int)
        if weight == "ren":
            self.cluster_negative_position_max, self.cluster_negative_position_all_n, \
            self.cluster_negative_position_all_p = self.cluster_position(
                self.FOFODF_K)
        elif weight == "ursby":
            self.cluster_negative_position_max, self.cluster_negative_position_all_n, self.cluster_negative_position_all_p \
                = self.cluster_position(
                self.FOFODF_Q)
        else:
            self.cluster_negative_position_max, self.cluster_negative_position_all_n, self.cluster_negative_position_all_p \
                = self.cluster_position(self.FOFODF)

        if occupancy:
            if weight == "ren":
                self.Extraplated_map_F_cal_occupancy(self.wDF, "2mFO-DFC", 1, "ren")
                self.mtz_extrapolated(self.wDF, "Extrapolated_map_all_ursby.mtz")
            elif weight == "ursby":
                self.Extraplated_map_F_cal_occupancy(self.rDF, "2mFO-DFC", 1, "ursby")
                self.mtz_extrapolated(self.rDF, "Extrapolated_map_all_ursby.mtz")
            else:
                self.Extraplated_map_F_cal_occupancy(self.delta_F, "2mFO-DFC", 1, "default")


    def check_sym(self):
        if self.F_dark_miller.space_group_number() == self.F_dark_miller.space_group_number():
            dark_vol = self.F_dark_miller.unit_cell().volume()
            light_vol = self.F_light_miller.unit_cell().volume()
            volume_change = (light_vol - dark_vol) / dark_vol
            volume_change_percent = volume_change * 100
            if volume_change_percent < 10:
                return False
            else:
                True
        else:
            True

    def write_extrapolate_map(self, weight, number):

        if weight == "ren":
            result = flex.double(self.F_dark + self.wDF)
            self.map_coefficient(result, coefficient="2mFO-DFC", write=True,
                                 name=f"%s/2mFO-DFC_%s_%s.ccp4" % (self.path, number, weight), map_type="Normal")

        elif weight == "ursby":
            result = flex.double(self.F_dark + self.rDF)
            self.map_coefficient(result, coefficient="2mFO-DFC", write=True,
                                 name=f"%s/2mFO-DFC_%s_%s.ccp4" % (self.path, number, weight), map_type="Normal")

        else:
            result = flex.double(self.F_dark + number * self.delta_F)
            self.map_coefficient(result, coefficient="2mFO-DFC", write=True,
                                 name=f"%s/2mFO-DFC_%s.ccp4" % (self.path, number), map_type="Normal")

    def cluster_position(self, map_coefficient):

        from residem.MapradiusSelection import MapPeakSelection_coff
        map_model = map_model_manager(map_manager=map_coefficient.as_map_manager(), model=self.model)
        self.map_coeff = MapPeakSelection_coff(map_model, self.high_resolution, self.sigma, write=True, path=self.path)

        all_out_n = reduce(operator.concat, self.map_coeff.cluster_n)
        all_array_n = flex.vec3_double(np.concatenate(tuple([x for x in all_out_n])))

        all_out_p = reduce(operator.concat, self.map_coeff.cluster_p)
        all_array_p = flex.vec3_double(np.concatenate(tuple([x for x in all_out_p])))

        cluster_num = self.map_coeff.cluster_num
        cluster_array = flex.vec3_double(self.map_coeff.cluster_n[cluster_num][0])
        fractional_array_max = self.map_coeff.uc.fractionalize(cluster_array)
        fractional_array_all_n = self.map_coeff.uc.fractionalize(all_array_n)
        fractional_array_all_p = self.map_coeff.uc.fractionalize(all_array_p)
        # peak_height = make_parallel(pk.real_map_unpadded.value_at_closest_grid_point)(pk.uc.fractionalize(flex.vec3_double(pk.cluster_n[60][0])))
        # peak_height_integral = np.array(peak_height).sum()
        return fractional_array_max, fractional_array_all_n, fractional_array_all_p

    def Extraplated_map_F_cal_occupancy(self, array_data_double, cof_type="Fo", step=1, weight="default"):
        list_range_occu = [(i, 200 / i) for i in range(1, 101, step) if i != 0]
        print("calculation of peak sum ")

        def cal_extrapole(x):
            map_data = self.f_model_dark_s.f_calc().as_amplitude_array().data() + x[
                1] * array_data_double  # Extrapolation
            self.F_cal_extra[f'F_cal_extra_%02d' % x[0]] = self.map_coefficient(map_data, cof_type, write=False,
                                                                                name=f'F_cal_extra_%02d.ccp4' % x[0])
            result = make_parallel(
                self.F_cal_extra['F_cal_extra_%02d' % x[0]].real_map_unpadded().value_at_closest_grid_point)(
                self.cluster_negative_position_max)
            result_sum = np.array(result, dtype=np.float32).sum()
            self.F_cal_extra[f'F_cal_extra_all%02d' % x[0]] = self.map_coefficient(map_data, cof_type)
            result_all_n = make_parallel(
                self.F_cal_extra['F_cal_extra_all%02d' % x[0]].real_map_unpadded().value_at_closest_grid_point)(
                self.cluster_negative_position_all_n)
            result_sum_all_n = np.array(result_all_n, dtype=np.float32).sum()

            result_all_p = make_parallel(
                self.F_cal_extra['F_cal_extra_all%02d' % x[0]].real_map_unpadded().value_at_closest_grid_point)(
                self.cluster_negative_position_all_p)
            result_sum_all_p = np.array(result_all_p, dtype=np.float32).sum()
            return [f'F_cal_extra_%02d' % x[0], result_sum, result_sum_all_n, result_sum_all_p]

        def obs_extrapole(x):
            map_data = self.f_model_dark_s.f_obs().as_amplitude_array().data() + x[
                1] * array_data_double  # Extrapolation
            self.F_cal_extra[f'F_obs_extra_%02d' % x[0]] = self.map_coefficient(map_data, cof_type, write=False,
                                                                                name=f'F_obs_extra_%02d.ccp4' % x[0])
            result = make_parallel(
                self.F_cal_extra['F_obs_extra_%02d' % x[0]].real_map_unpadded().value_at_closest_grid_point)(
                self.cluster_negative_position_max)
            result_sum = np.array(result).sum()
            self.F_cal_extra[f'F_obs_extra_all%02d' % x[0]] = self.map_coefficient(map_data, cof_type)
            result_all_n = make_parallel(
                self.F_cal_extra['F_obs_extra_all%02d' % x[0]].real_map_unpadded().value_at_closest_grid_point)(
                self.cluster_negative_position_all_n)
            result_sum_all_n = np.array(result_all_n).sum()

            result_all_p = make_parallel(
                self.F_cal_extra['F_obs_extra_all%02d' % x[0]].real_map_unpadded().value_at_closest_grid_point)(
                self.cluster_negative_position_all_p)
            result_sum_all_p = np.array(result_all_p).sum()
            return [f'F_obs_extra_%02d' % x[0], result_sum, result_sum_all_n, result_sum_all_p]

        warnings.filterwarnings("ignore")

        def logistic_growth(x, K, r, t0):
            return K / (1 + np.exp(-r * (x - t0)))

        def derivative(x, K, r, t0):
            return K * r * (1 - logistic_growth(x, K, r, t0) / K) * (np.exp(-r * (x - t0))) / (
                    1 + np.exp(-r * (x - t0))) ** 2

        result_cal = make_parallel(cal_extrapole)(tqdm.tqdm(list_range_occu))
        x_cal, y_cal, z_cal, i_cal = [int(x.split('_')[-1]) for x, y, z, i in result_cal], [y for x, y, z, i in
                                                                                            result_cal], [z for
                                                                                                          x, y, z, i
                                                                                                          in
                                                                                                          result_cal], [
                                         i for
                                         x, y, z, i
                                         in
                                         result_cal]
        result_obs = make_parallel(cal_extrapole)(tqdm.tqdm(list_range_occu))
        x_obs, y_obs, z_obs, i_obs = [int(x.split('_')[-1]) for x, y, z, i in result_obs], [y for x, y, z, i in
                                                                                            result_obs], [z for
                                                                                                          x, y, z, i
                                                                                                          in
                                                                                                          result_obs], [
                                         i for
                                         x, y, z, i
                                         in
                                         result_obs]

        df = pd.DataFrame()
        df["Percentage_cal"] = x_cal
        df["sum_negative_cal"] = y_cal
        df["sum_negative_all_cal"] = z_cal
        df["sum_positive_all_cal"] = i_cal
        x_axis = [x for x in range(1, 101)]
        df.sort_values(by=["Percentage_cal"], inplace=True)
        df1 = pd.DataFrame()
        df1["Percentage_obs"] = x_obs
        df1["sum_negative_obs"] = y_obs
        df1["sum_negative_all_obs"] = z_obs
        df1["sum_positive_all_obs"] = i_obs
        df1.sort_values(by=["Percentage_obs"], inplace=True)

        # initialParameters = np.array([1.0, 1.0, 1.0, 1.0])
        # if df["z"][0] > 0:
        #     data = df["z"]
        #     fittedParameters, pcov = curve_fit(func, x_axis, df["z"], initialParameters)
        # elif df["i"][0] > 0:
        #     data = df["i"]

        out_pdf = f'%s/estimated_occupancy_%s.pdf' % (self.path, weight)

        data = {r"sum negative max peak cal[$\Delta F(extra)+ F_C$]": df["sum_negative_cal"],
                r"sum negative all peaks cal[$\Delta F(extra)+ F_C$]": df["sum_negative_all_cal"],
                r"sum positive all peaks cal[$\Delta F(extra)+ F_C$]": df["sum_positive_all_cal"],
                r"sum negative max peak obs[$\Delta F(extra)+ F(obs)$]": df1["sum_negative_obs"],
                r"sum negative all peaks obs[$\Delta F(extra)+ F(obs)$]": df1["sum_negative_all_obs"],
                r"sum positive all peaks obs[$\Delta F(extra)+ F(obs)$]": df1["sum_positive_all_obs"]
                }
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        for x, y in data.items():
            if (y < 0).sum() > len(y) / 2:
                y = np.negative(y)

            if y.min() < 0:
                min_val = np.abs(np.min(y))
                y = np.add(y, min_val)
            fig, ax = plt.subplots(figsize=(10, 10))
            p0 = [y.max(), 0, 1]
            popt, _ = curve_fit(logistic_growth, x_axis, list(y), p0, maxfev=1000)
            model = logistic_growth(x_axis, *popt)
            K, r, t0 = popt
            if K < 0:
                popt, _ = curve_fit(logistic_growth, x_axis, list(np.array(y)), p0, maxfev=1000)
                model = logistic_growth(x_axis, *popt)
                K, r, t0 = popt
                ax.scatter(x_axis, y, label="Raw data")

            else:
                ax.scatter(x_axis, y, label="Raw data")
            growth_rate = derivative(x_axis, K, r, t0)
            p1 = [growth_rate.max(), 0, 1]
            popt1, _ = curve_fit(logistic_growth, x_axis, growth_rate, p1)
            K1, r1, t1 = popt1
            growth_rate_1 = derivative(x_axis, K1, r1, t1)
            if K1 < 0:
                res = minimize_scalar(lambda x: -derivative(x, K1, r1, t1))
                x_max = round(res.x)
            else:
                res = minimize_scalar(lambda x: derivative(x, K1, r1, t1))
                x_max = round(res.x)
            r2 = r2_score(y, model)
            ax.plot(x_axis, model, label=r'Logistic growth\nFitted model, $R^2=%.2f$' % r2, color='orange',
                    linewidth=2.5)

            ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
            ax.set_xlabel('Occupancy')
            ax.set_ylabel('Peak Height')
            if 100 > x_max > 0:
                ax.axvline(x=x_max, color='y', alpha=0.4, label=f"Estimated occupancy %s" % x_max)
                self.write_extrapolate_map(self.weight, x_max)

            axins = inset_axes(ax, width="30%", height="30%", loc=1)
            axins.scatter(x_axis, growth_rate_1, label="Logistic growth-\ncurve derivative", color='g', alpha=0.2)
            axins.legend(bbox_to_anchor=(-0.1, -0.25), loc='upper left', prop={'size': 9})
            plt.xticks(visible=True)
            plt.yticks(visible=True)
            ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
            ax.set_title("Occupancy Estimation(%s)" % x)
            plt.rcParams.update({'font.size': 8})
            pdf.savefig(fig, bbox_inches='tight')
        pdf.close()
        df.to_csv(f'%s/estimated_occupancy_%s_cal.csv' % (self.path, weight))
        df1.to_csv(f'%s/estimated_occupancy_%s_obs.csv' % (self.path, weight))
        return df, df1

    def map_model_scaling(self, f_model):
        obs = f_model.f_obs()
        f_obs_scale = 1.0 / f_model.k_anisotropic() / f_model.k_isotropic()
        obs = miller.array(miller_set=self.miller_set, data=obs.data() * f_obs_scale)
        obs.set_sigmas(f_model.f_obs().sigmas() * f_obs_scale)
        return obs

    def linear_scaling(self, f_light, f_dark):
        scale_k1 = 1
        den = flex.sum(flex.abs(f_light.data()) * flex.abs(f_light.data()))
        if den != 0:
            scale_k1 = flex.sum(flex.abs(f_dark.data()) * flex.abs(f_dark.data())) / den
        F_light_common_scaled_data = f_light.array(
            data=f_light.data() * scale_k1, sigmas=f_light.sigmas() * scale_k1)
        return F_light_common_scaled_data

    def aniso_scaling(self, f_light, f_dark):
        F_light_common_scaled_data = f_light.array(
            data=flex.double(self.exp_anisotropic_scaling(f_dark.data(), f_light.data(),
                                                          np.array(self.indices.as_vec3_double().as_numpy_array(),
                                                                   dtype=np.float64))),
            sigmas=flex.double(self.exp_anisotropic_scaling(f_dark.sigmas(),
                                                            f_light.sigmas(),
                                                            np.array(self.indices.as_vec3_double().as_numpy_array(),
                                                                     dtype=np.float64))))
        return F_light_common_scaled_data

    @staticmethod
    def aniso_exp_residuals(params, F_native, F_deriv, Miller_indices):
        C, B11, B22, B33, B12, B13, B23 = params
        h, k, l = Miller_indices.T
        scaling_factors = C * np.exp(
            -(h ** 2 * B11 + k ** 2 * B22 + l ** 2 * B33 + 2 * h * k * B12 + 2 * h * l * B13 + 2 * k * l * B23))
        F_deriv_scaled = F_deriv * scaling_factors
        residuals = np.abs(F_native - F_deriv_scaled)
        return residuals

    def exp_anisotropic_scaling(self, F_native, F_deriv, Miller_indices):
        initial_params = np.array([1.0] + [0.0] * 6)
        result = least_squares(self.aniso_exp_residuals, initial_params, args=(F_native, F_deriv, Miller_indices))
        C, B11, B22, B33, B12, B13, B23 = result.x
        h, k, l = Miller_indices.T
        scaling_factors = C * np.exp(
            -(h ** 2 * B11 + k ** 2 * B22 + l ** 2 * B33 + 2 * h * k * B12 + 2 * h * l * B13 + 2 * k * l * B23))
        F_deriv_scaled = F_deriv * scaling_factors
        return F_deriv_scaled

    def map_model_multiscaling(self, f_light, f_dark):
        fobs_dark = f_light.multiscale(other=f_dark, reflections_per_bin=250)
        return fobs_dark

    def phase_transfer(self, miller_array, phase_source):
        tmp = miller.array(miller_set=self.miller_set, data=flex.double(miller_array.indices().size(), 1)
                           ).phase_transfer(phase_source=phase_source)
        return miller.array(miller_set=self.miller_set, data=miller_array.data() * tmp.data())

    ## difference map map coefficint
    def diff_map(self, f_model, f_light, f_dark):
        from cctbx import miller
        diff = miller.array(
            miller_set=self.miller_set,
            data=f_light.data() - f_dark.data())
        map_coeff = self.phase_transfer(
            miller_array=diff,
            phase_source=f_model.f_model())
        return map_coeff

    ## weight calculation
    ### weight ##

    def W_weight(self):
        """ weight as suggested by Ren 2001 and hekstra 2016 """
        # SigDF = np.sqrt(self.SIGF_dark ** 2 + self.SIGF_light ** 2)
        SigDF = np.sqrt(self.f_model_dark_scaled.sigmas() ** 2 + self.f_model_light_scaled_obs.sigmas() ** 2)
        weight = self.compute_target_weight_k(self.delta_F.as_numpy_array(), SigDF)
        return weight

    @staticmethod
    def compute_target_weight_k(F, SIGF, alpha=1):
        ''' This is a simplified weighting scheme implemented by Ren 2001 and hekstra 2016.'''

        w = (1 + (SIGF ** 2 / (SIGF ** 2).mean()) + alpha * (F ** 2 / (F ** 2).mean()))
        return w ** -1

    def compute_q_weight(self):
        """ This is the weighting scheme implemented by ursby 1997 also implemented in xtrapol8 2022 """

        mult = self.f_model_dark_scaled.multiplicities().data().as_double()
        mult = miller.array(miller_set=self.miller_set.centric_flags(), data=mult)
        # stol = self.F_dark_miller.unit_cell().stol(self.indices)  # sin(theta)/lambda
        # res = 1 / (2 * stol)  # RESOLUTION
        # res = miller.array(miller_set=self.miller_set.centric_flags(), data=res)
        num = flex.double([2 if val == True else 1 for val in self.miller_set.centric_flags().data()])
        num = miller.array(miller_set=self.miller_set.centric_flags(), data=num)
        indices = flex.miller_index()
        weight = flex.double()
        self.f_model_dark_scaled.setup_binner(n_bins=20)
        self.f_model_light_scaled_obs.use_binning_of(self.f_model_dark_scaled)
        mult.use_binning_of(self.f_model_dark_scaled)
        # res.use_binning_of(self.f_model_dark_scaled)
        num.use_binning_of(self.f_model_dark_scaled)
        for bin in self.f_model_dark_scaled.binner().range_used():
            sel_miller_bin = self.f_model_dark_scaled.binner().selection(bin)
            sel_dark_data = self.f_model_dark_scaled.select(sel_miller_bin).deep_copy()
            sel_light_data = self.f_model_light_scaled_obs.select(sel_miller_bin).deep_copy()
            sel_mult_bin = mult.select(sel_miller_bin)
            sel_num_bin = num.select(sel_miller_bin)
            indices_bin = sel_dark_data.indices()
            if sel_dark_data.size() == 0: continue
            F_light_minus_F_dark = sel_light_data.data() - sel_dark_data.data()
            abs_F_light_minus_dark = flex.abs(F_light_minus_F_dark)
            sigma_sum_square = sel_light_data.sigmas() ** 2 + sel_dark_data.sigmas() ** 2
            sigma_square_ND = flex.mean((sel_num_bin.data() * F_light_minus_F_dark ** 2 / sel_mult_bin.data()))
            sigma_square_D = flex.mean(
                (sel_num_bin.data() * (F_light_minus_F_dark ** 2 - sigma_sum_square) / sel_mult_bin.data()))

            if sigma_square_D > sigma_square_ND / 2:
                numerator = sel_mult_bin.data() * sigma_square_D / sel_num_bin.data()
                w_bin = numerator / (sigma_sum_square + numerator)

            else:
                numerator = sel_mult_bin.data() * 0.5 * sigma_square_ND / sel_num_bin.data()
                w_bin = numerator / (sigma_sum_square + numerator)

            sel = w_bin * abs_F_light_minus_dark > 3 * np.sqrt(flex.mean(F_light_minus_F_dark ** 2))
            w_bin = w_bin.set_selected(sel, 0)
            indices = indices.concatenate(indices_bin)
            weight = weight.concatenate(w_bin)

        q_weight = miller.array(miller_set=miller.set(self.cs, indices, anomalous_flag=False),
                                data=weight)

        if q_weight.indices().size() != self.indices.size():
            raise ValueError("The size of miller not matching the q weight")

        else:
            return q_weight

        def compute_hekstra_weight(self,weight_term,alpha):
            mult = self.f_model_dark_scaled.multiplicities().data().as_double()
            mult = miller.array(miller_set=self.miller_set.centric_flags(), data=mult)
            # stol = self.F_dark_miller.unit_cell().stol(self.indices)  # sin(theta)/lambda
            # res = 1 / (2 * stol)  # RESOLUTION
            # res = miller.array(miller_set=self.miller_set.centric_flags(), data=res)
            num = flex.double([2 if val == True else 1 for val in self.miller_set.centric_flags().data()])
            num = miller.array(miller_set=self.miller_set.centric_flags(), data=num)
            indices = flex.miller_index()
            weight = flex.double()
            self.f_model_dark_scaled.setup_binner(n_bins=20)
            self.f_model_light_scaled_obs.use_binning_of(self.f_model_dark_scaled)
            mult.use_binning_of(self.f_model_dark_scaled)
            # res.use_binning_of(self.f_model_dark_scaled)
            num.use_binning_of(self.f_model_dark_scaled)
            SigDF = self.f_model_dark_scaled.sigmas() ** 2 + self.f_model_light_scaled_obs.sigmas() ** 2
            # sqSigDF = flex.double(np.sqrt(self.f_model_dark_scaled.sigmas() ** 2 + self.f_model_light_scaled_obs.sigmas() ** 2))


            DF_mean = flex.mean(self.delta_F ** 2)
            Sig_mean = flex.mean(SigDF)



            for bin in self.f_model_dark_scaled.binner().range_used():
                sel_miller_bin = self.f_model_dark_scaled.binner().selection(bin)
                sel_dark_data = self.f_model_dark_scaled.select(sel_miller_bin).deep_copy()
                sel_light_data = self.f_model_light_scaled_obs.select(sel_miller_bin).deep_copy()
                indices_bin = sel_dark_data.indices()
                if sel_dark_data.size() == 0: continue
                F_light_minus_F_dark = sel_light_data.data() - sel_dark_data.data()
                sigma_sum_square = sel_light_data.sigmas() ** 2 + sel_dark_data.sigmas() ** 2
                SigDF_op = sigma_sum_square / Sig_mean
                DF_term = alpha * (F_light_minus_F_dark **2 / DF_mean  )
                w_bin  = 1/ (1+SigDF_op + DF_term)
                indices = indices.concatenate(indices_bin)
                weight = weight.concatenate(w_bin)

            r_weight = miller.array(miller_set=miller.set(self.cs, indices, anomalous_flag=False),
                                    data=weight / flex.mean(flex.double(weight)))

            if r_weight.indices().size() != self.indices.size():
                raise ValueError("The size of miller not matching the q weight")

            else:
                return r_weight

    def mtz_extrapolated(self, data_array,name="Extrapolated_map_all.mtz"):

        mtz_dataset_extra = miller.array(miller_set=self.miller_set,
                                         data=self.f_model_dark_scaled.as_amplitude_array().data()).as_mtz_dataset(
            column_root_label="F_DARK", column_types="F")
        mtz_dataset_extra.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=self.f_model_dark_scaled.sigmas())
            , column_root_label="SIGF_DARK", column_types="Q")
        mtz_dataset_extra.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                                     data=self.f_model_light_scaled_obs.data()),
                                           column_root_label="F_LIGHT", column_types="F")

        mtz_dataset_extra.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=self.f_model_light_scaled_obs.sigmas())
            , column_root_label="SIGF_LIGHT", column_types="Q")

        mtz_dataset_extra.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set,
                                      data=flex.double(self.F_cal_common.phases(deg=True).data())),
            column_root_label="PHIC", column_types="P")

        list_range_occu = [(i, 200 / i) for i in range(1, 101, 1) if i != 0]

        

        def add_extra(x):
            mtz_dataset_extra.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                                         data=self.f_model_dark_s.f_calc().as_amplitude_array().data() +
                                                                              x[1] * data_array),
                                               column_root_label=f"F_extrapolated_%s" % x[0], column_types="F")
        make_parallel(add_extra)(tqdm.tqdm(list_range_occu))
        mtz_dataset_extra.mtz_object().write(name)

    def mtz_dataset(self, name="Difference_map_weighted_all.mtz"):
        mtz_dataset = miller.array(miller_set=self.miller_set,
                                   data=self.f_model_dark_scaled.as_amplitude_array().data()).as_mtz_dataset(
            column_root_label="F_DARK", column_types="F")
        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=self.f_model_dark_scaled.sigmas())
            , column_root_label="SIGF_DARK", column_types="Q")
        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                               data=self.f_model_light_scaled_obs.data()),
                                     column_root_label="F_LIGHT", column_types="F")

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=self.f_model_light_scaled_obs.sigmas())
            , column_root_label="SIGF_LIGHT", column_types="Q")

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set,
                                      data=flex.double(self.F_cal_common.phases(deg=True).data())),
            column_root_label="PHIC", column_types="P")

        mtz_dataset.add_miller_array(miller_array=self.f_model_dark.f_calc().as_amplitude_array(),
                                     column_root_label="F_CAL")
        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                               data=flex.double(self.delta_F)),
                                     column_root_label="DF", column_types="F")

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.weight_k)),
            column_root_label='Ren_weight', column_types="W")

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.q_weight)),
            column_root_label='Ursby_weight', column_types="W")
        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.wDF)),
                                     column_root_label='kDF', column_types='F')
        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.rDF)),
                                     column_root_label='qDF', column_types='F')

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.delta_F_model)),
            column_root_label='DF_cal', column_types='F')
        mtz_dataset.mtz_object().write(name)

        return mtz_dataset

    ### crystal gridding
    def crystal_gridding(self):
        """crystal grid parameter"""
        resolution_factor = self.grid_step / self.high_resolution
        mp = mmtbx.masks.mask_master_params.extract()
        mp.grid_step_factor = 1. / resolution_factor
        mmtbx_masks_asu_mask_obj = mmtbx.masks.asu_mask(
            xray_structure=self.xrs,
            d_min=self.high_resolution,
            mask_params=mp)

        bulk_solvent_mask = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()

        sel = bulk_solvent_mask > 0

        bulk_solvent_mask = bulk_solvent_mask.set_selected(sel, 1)

        cr_gr = maptbx.crystal_gridding(
            unit_cell=self.xrs.unit_cell(),
            space_group_info=self.xrs.space_group_info(),
            pre_determined_n_real=bulk_solvent_mask.focus())

        return cr_gr

    def map_coefficient(self, data, coefficient="FO", write=False, name="map.ccp4", map_type=None):
        # dm = DataManager()
        # dm.set_overwrite(True)
        # miller_da = miller.array(miller_set=self.miller_set, data=flex.double(data))
        # miller_da = miller_da.phase_transfer(self.f_model_dark_s.f_model().phases(deg=True))
        # miller_da = miller_da.map_to_asu().average_bijvoet_mates()
        f_model = mmtbx.f_model.manager(f_obs=miller.array(miller_set=self.miller_set,
                                                           data=data), xray_structure=self.xrs)

        map_coefficients = map_tools.electron_density_map(fmodel=f_model).map_coefficients(
            map_type=coefficient, isotropize=False, fill_missing=False, sharp=False)

        fft_map = miller.fft_map(crystal_gridding=self.crystal_gridding(),
                                 fourier_coefficients=map_coefficients.as_non_anomalous_array())
        if map_type == "diff":
            miller_da = miller.array(miller_set=self.miller_set, data=flex.double(data))
            miller_da = miller_da.phase_transfer(f_model.f_model())
            miller_da = miller_da.map_to_asu().average_bijvoet_mates()
            fft_map = miller_da.fft_map(
                resolution_factor=self.resolution_factor)  # crystal_gridding=self.crystal_gridding(),
            # resolution_factor=self.resolution_factor,
            # symmetry_flags=maptbx.use_space_group_symmetry).apply_sigma_scaling()

            # if self.inverse_scaling:
            #     fft_map = fft_map.apply_scaling(-1)
            # fft_map = fft_map.apply_sigma_scaling()
        else:
            fft_map = fft_map.apply_sigma_scaling()

        if write:
            # miller_scaled = miller_da.fft_map(crystal_gridding=self.crystal_gridding(),
            #                                   resolution_factor=self.resolution_factor,
            #                                   symmetry_flags=maptbx.use_space_group_symmetry).apply_sigma_scaling()
            # miller_scaled.as_ccp4_map(name)
            fft_map.as_ccp4_map(f"{name}")
            # mmm = map_model_manager(  # a new map_model_manager
            #     model=self.model,  # initializing with a model
            #     map_manager=miller_scaled.as_map_manager())
            # box_mmm = mmm.extract_all_maps_around_model(  # extract a box around model
            #     selection_string="all")
            # dm.write_real_map_file(box_mmm.map_manager(),
            #                        filename="Data_folder/" + pathlib.Path(name).stem + "_box.ccp4")
        return fft_map


def run(ref, mtz, tmtz, refl, tril, high, low, sigma, grid, multiscaling, weight, phil_occu, inverse, scale):
    path = os.getcwd()
    sc = ScaleModel(reference_pdb=ref, reference_mtz=mtz, Triggered_mtz=tmtz, reference_label=refl,
                    triggered_label=tril, high_resolution=high, low_resolution=low, sigma=sigma, grid_step=grid,
                    multiscaling=multiscaling, weight=weight, occupancy=phil_occu, path=path, inverse_scaling=inverse,
                    scale=scale)


def run_all():
    import argparse
    import sys
    from datetime import datetime

    start = datetime.now()
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--ref", type=str, help="reference pdb file")
    parser.add_argument("-m", "--mtz", type=str, help="reference hkl file in mtz format")
    parser.add_argument("-t", "--tmtz", type=str, help="triggered hkl file in mtz format")
    parser.add_argument("-hi", "--high", type=float, help="high resolution")
    parser.add_argument("-l", "--low", type=float, help="low resolution")
    parser.add_argument("-g", "--grid", type=float, help="grid step")
    parser.add_argument("-rl", "--refl", action='store', help="label of the reference file ")
    parser.add_argument("-tl", "--tril", action='store', help="label of the triggred file")
    parser.add_argument("-s", "--sigma", type=float, help="sigma cutoff for map processing")
    parser.add_argument("-w", "--weight", type=str, help="weight for choosing the map calculation")
    parser.add_argument("-x", "--mul", type=bool, help="if multiscaling is true or false")
    parser.add_argument("-o", "--occu", type=bool, help="occupancy estimation")
    parser.add_argument("-i", "--inverse", type=bool, help="inverse scaling")
    parser.add_argument("--scale", type=str,
                        help="scaling if it is 'iso'(isotropic scaling) or 'aniso'(anisotropic scalig) or multi scaling")

    # parser.add_argument("-c", "--method", help="method to choose the occupancy estimation, cluster,map_atom,all")

    args = parser.parse_args()

    if len(sys.argv) == 3:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        refl = args.refl if args.refl is not None else None
        tril = args.tril if args.tril is not None else None
        high = args.high if args.high is not None else None
        low = args.low if args.low is not None else None
        sigma = args.sigma if args.sigma is not None else 3.0
        grid = args.sgrid if args.grid is not None else 0.5
        max_distance = args.mul if args.mul is not None else False
        phil_occu = args.occu if args.occu is not None else False
        inverse = args.inverse if args.inverse is not None else False
        scale = args.scale if args.scale in ["iso", "aniso"] else "aniso"
        # cluster = "cluster" if args.method is None else args.method if args.method.lower() in ["map_atom", "all",
        #                                                                                        "cluster"] else "cluster "
        weight = "default" if args.weight is None else args.weight.lower() if args.weight.lower() in ["ren",
                                                                                                      "ursby"] else "default"
        if args.ref and args.mtz and args.tmtz:
            print(scale)
            run(args.ref, args.mtz, args.tmtz, refl, tril, high, low, sigma, grid, max_distance, weight, phil_occu,
                inverse, scale)
        else:
            print("Please provide the reference pdb, reference hkl file in mtz format and triggered file in mtz format")

    print(datetime.now() - start)


if __name__ == '__main__':
    run_all()
