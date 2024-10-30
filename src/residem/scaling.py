import os
from residem.master import DataInput
from mmtbx.scaling import make_param
import iotbx.phil
import mmtbx.model
from mmtbx.scaling import matthews
from mmtbx.scaling import absolute_scaling, relative_scaling
from mmtbx.scaling import pair_analyses
from mmtbx.scaling import fa_estimation
from cctbx import miller, maptbx
import mmtbx.f_model
from cctbx.array_family import flex
from scipy.optimize import minimize, least_squares
import numpy as np
from mmtbx import map_tools
from iotbx.map_model_manager import map_model_manager
from libtbx.utils import null_out
import sys
import iotbx.map_tools
from cctbx.xray import structure
import matplotlib.pyplot as plt
from residem.logger_set import setup_logger
import sympy as sp
import pandas as pd

plt.rcParams['image.cmap'] = 'RdBu_r'
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

matplotlib.use('Agg')
font = {'weight': 'bold', 'size': 18, 'family': 'sans-serif', 'sans-serif': ['Helvetica']}

matplotlib.rc('font', **font)
import seaborn as sns


class ScaleNativeDerivative(DataInput):
    def __init__(self, reference_pdb, reference_mtz, Triggered_mtz, reference_label=None, triggered_label=None,
                 fcalclabel=None,
                 high_resolution=None, low_resolution=None, weight="default", grid_step=0.5, scale="iso", write=True,
                 path=None, alpha=1, logger="output.log"):
        super().__init__(reference_pdb, reference_mtz, Triggered_mtz, reference_label, triggered_label, fcalclabel,
                         high_resolution, low_resolution, logger)
        self.logger = setup_logger(logger)
        self.alpha = alpha

        self.model = mmtbx.model.manager(model_input=self.pdb_in)
        self.grid_step = grid_step
        self.resolution_factor = self.grid_step / self.high_resolution

        self.F_dark_miller = self.I_2_F(self.dark_dict[self.reference_label])
        self.F_light_miller = self.I_2_F(self.light_dict[self.triggered_label])
        self.F_cal = self.dark_dict[self.map_obs_label(reference_label, triggered_label, fcalclabel)[2]]

        info0 = self.F_dark_miller.info()
        info1 = self.F_light_miller.info()
        self.assertion = self.check_sym()
        self.path = path

        self.F_dark_miller = self.F_dark_miller.resolution_filter(d_min=self.high_resolution,
                                                                  d_max=self.low_resolution).set_info(info0)
        self.F_light_miller = self.F_light_miller.resolution_filter(d_min=self.high_resolution,
                                                                    d_max=self.low_resolution).set_info(info1)
        self.logger.info(f" ")
        self.logger.info(
            f"""Symmetry of reference structure is {self.F_dark_miller.crystal_symmetry()}\nSymmetry of Triggered structure is {self.F_light_miller.crystal_symmetry()}""")

        if not self.check_sym():
            self.F_light_miller = self.F_light_miller.customized_copy(
                crystal_symmetry=self.F_dark_miller.crystal_symmetry())

        self.dark_miller_native = self.Native_miller(self.F_dark_miller)
        self.light_miller_native = self.Native_miller(self.F_light_miller)
        self.dark_before_common = self.F_dark_miller.common_set(self.F_light_miller)
        self.light_before_common = self.F_light_miller.common_set(self.F_dark_miller)
        self.logger.info(f"""Selecting common reflection between reference and triggered state""")

        if not self.check_sym():
            self.light_miller_native = self.light_miller_native.customized_copy(
                crystal_symmetry=self.dark_miller_native.crystal_symmetry())

        self.xrs = self.pdb_in.xray_structure_simple()
        # self.light_xrs = structure(scatterers=self.xrs.scatterers(),
        #                            crystal_symmetry=self.F_light_miller.crystal_symmetry())

        self.native_pre_scale = Prescaling(reference_pdb, self.dark_miller_native).miller_array_scaled

        self.derivative_pre_scale = Prescaling(reference_pdb, self.light_miller_native).miller_array_scaled

        self.params = self.param_4_sf_correction()
        print("Scaling the Data")
        self.logger.info(f"Scaling the Data")
        sys.stdout = open(os.devnull, "w")
        self.miller_array_native, self.miller_array_derivative = self.fa_estimate(self.native_pre_scale,
                                                                                  self.derivative_pre_scale)

        sys.stdout = sys.__stdout__

        self.miller_array_native = self.miller_array_native.common_set(self.miller_array_derivative)
        self.miller_array_derivative = self.miller_array_derivative.common_set(self.miller_array_native)

        self.f_model_dark = mmtbx.f_model.manager(f_obs=self.miller_array_native, xray_structure=self.xrs)
        self.f_model_dark.update_all_scales()
        self.f_model_light = mmtbx.f_model.manager(f_obs=self.miller_array_derivative, xray_structure=self.xrs)
        self.f_model_light.update_all_scales()

        self.f_obs_dark_common = self.f_model_dark.f_obs().common_set(self.f_model_light.f_obs(),
                                                                      assert_is_similar_symmetry=self.assertion)
        self.f_obs_light_common = self.f_model_light.f_obs().common_set(self.f_model_dark.f_obs(),
                                                                        assert_is_similar_symmetry=self.assertion)
        self.f_model_dark_common = mmtbx.f_model.manager(f_obs=self.f_obs_dark_common, xray_structure=self.xrs)
        self.f_model_dark_common.update_all_scales()
        self.f_model_light_common = mmtbx.f_model.manager(f_obs=self.f_obs_light_common, xray_structure=self.xrs)
        self.f_model_light_common.update_all_scales()
        self.F_cal_common = self.F_cal.common_set(self.f_model_dark_common.f_obs())

        self.indices = self.f_obs_dark_common.indices()
        self.cs = self.f_model_dark.f_model().crystal_symmetry()
        self.uc = self.cs.unit_cell()
        self.miller_set = miller.set(self.cs, self.indices, anomalous_flag=False)

        self.f_model_dark_scaled = self.map_model_scaling(self.f_model_dark_common, scale)
        self.f_model_light_scaled = self.map_model_scaling(self.f_model_light_common, scale)
        self.F_dark = self.f_model_dark_scaled.data().as_numpy_array()
        self.SIGF_dark = self.f_model_dark_scaled.sigmas().as_numpy_array()
        self.F_light = self.f_model_light_scaled.data().as_numpy_array()
        self.SIGF_light = self.f_model_light_scaled.sigmas().as_numpy_array()
        self.F_calc = self.F_cal_common.data().as_numpy_array()
        self.R_Free_array = self.f_model_dark_scaled.generate_r_free_flags(
            fraction=0.05, use_lattice_symmetry=True, format="shelx").data().as_numpy_array()
        self.R_Free_array[self.R_Free_array == -1] = 0
        self.R_Free = miller.array(miller_set=self.miller_set,
                                   data=flex.int(self.R_Free_array))

        # scaling
        if scale == "iso":
            self.logger.info(f"Isotropic scaling ")
            self.f_model_light_scaled_obs = self.linear_scaling(self.f_model_light_scaled, self.f_model_dark_scaled)

        else:
            self.logger.info(f"Anisotropic scaling ")
            self.f_model_light_scaled_obs = self.aniso_scaling(self.f_model_light_scaled, self.f_model_dark_scaled)

        ### differnce map
        if scale == "no":
            self.logger.info(f"No scaling ")
            self.f_model_dark_scaled = self.F_dark_miller.common_set(self.f_obs_dark_common)
            self.f_model_light_scaled_obs = self.F_light_miller.common_set(self.f_obs_light_common)

            self.delta_F = self.f_model_light_scaled_obs.data() - self.f_model_dark_scaled.data()

        else:
            self.delta_F = self.f_model_light_scaled_obs.data() - self.f_model_dark_scaled.data()
        # self.delta_F_model = self.F_light - self.F_calc
        self.phase = self.phase_flip(self.delta_F)

        # ## weight calculation
        # ## K weightSC.F
        self.logger.info(f"Calculating different weights in isomorphous difference density map")
        self.weight_ren = self.W_weight("ren")
        self.ursby_weight = self.compute_ursby_weight().data()

        self.sigma_df = self.sigma_cal()

        self.ren_factor = self.weight_ren / flex.mean(flex.double(self.weight_ren))
        self.ursby_factor = self.ursby_weight / flex.mean(flex.double(self.ursby_weight))
        self.weight_hekstra = self.W_weight("hekstra", alpha=alpha)  # .data()

        self.hekstra_factor = self.weight_hekstra / flex.mean(flex.double(self.weight_hekstra))
        self.hDF = flex.double(self.hekstra_factor * self.delta_F)
        self.wDF = flex.double(self.ren_factor.as_numpy_array() * self.delta_F)
        self.rDF = flex.double(self.ursby_factor.as_numpy_array() * self.delta_F)
        self.sigma_ren_df = self.sigma_df * self.ren_factor
        self.sigma_ursby_df = self.sigma_df * self.ursby_factor
        self.sigma_hekstra_df = self.sigma_df * self.hekstra_factor
        if self.alpha != 1:
            self.weight_hekstra_1 = self.W_weight("hekstra", alpha=1)
            self.hekstra_factor_1 = self.weight_hekstra_1 / flex.mean(flex.double(self.weight_hekstra_1))
            self.hDF_1 = flex.double(self.hekstra_factor_1 * self.delta_F)
            self.sigma_hekstra_df_1 = self.sigma_df * self.hekstra_factor_1

    def I_2_F(self, miller_array):
        if miller_array.is_xray_intensity_array():
            miller_array_copy = miller_array.deep_copy()
            miller_array_copy_I = miller_array_copy.french_wilson()
        elif miller_array.is_xray_amplitude_array():
            miller_array_copy_I = miller_array.deep_copy()

        return miller_array_copy_I

    def sigma_cal(self):
        Sigma_DF = np.sqrt(self.f_model_dark_scaled.sigmas() ** 2 + self.f_model_light_scaled_obs.sigmas() ** 2)
        return Sigma_DF

    @staticmethod
    def Native_miller(miller_in):
        miller_out = miller_in.map_to_asu().select(miller_in.indices() != (0, 0, 0)).deep_copy()
        if miller_out.is_xray_intensity_array():
            miller_out = miller_out.select(
                miller_out.data() > 0)
            miller_out = miller_out.f_sq_as_f()

        if miller_out.anomalous_flag():
            miller_out = miller_out.average_bijvoet_mates() \
                .set_observation_type(miller_out)
        return miller_out

    def param_4_sf_correction(self):
        """:parameter of default single isomorphus replacement is used"""
        params_generator = make_param.phil_lego()
        master_params = iotbx.phil.parse(params_generator.default_sir())
        effective_params = master_params.fetch()
        params = effective_params.extract()
        return params

    def fa_estimate(self, miller_array_native, miller_array_derivative):
        """Estimation of the structure factor correction factor amplitude using maximum liklihood method"""

        scalar = fa_estimation.combined_scaling(miller_array_native, miller_array_derivative,
                                                self.params.scaling.input.scaling_strategy.iso_protocol, out=null_out())
        miller_array_native = scalar.x1.deep_copy()
        miller_array_derivative = scalar.x2.deep_copy()
        del scalar
        return miller_array_native, miller_array_derivative

    def map_model_scaling(self, f_model, scale="no"):

        obs = f_model.f_obs()
        f_obs_scale = 1.0 / f_model.k_anisotropic() / f_model.k_isotropic()
        if scale == "no":
            f_obs_scale = 1
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

    def map_coefficient(self, data, coefficient="FO", write=False, name="map", map_type=None):

        f_model = mmtbx.f_model.manager(f_obs=miller.array(miller_set=self.miller_set,
                                                           data=data), xray_structure=self.xrs)

        map_coefficients = map_tools.electron_density_map(fmodel=f_model).map_coefficients(
            map_type=coefficient, isotropize=False, fill_missing=False, sharp=False, acentrics_scale=2.0,
            centrics_pre_scale=1.0)

        fft_map = miller.fft_map(crystal_gridding=self.crystal_gridding(),
                                 fourier_coefficients=map_coefficients.as_non_anomalous_array())
        if map_type == "diff":
            miller_da = miller.array(miller_set=self.miller_set, data=flex.double(data))
            miller_da = miller_da.phase_transfer(f_model.f_model())
            miller_da = miller_da.map_to_asu().average_bijvoet_mates()

            fft_map = miller_da.fft_map(
                resolution_factor=self.resolution_factor)
            fft_map = fft_map.apply_sigma_scaling()

        else:
            fft_map = fft_map.apply_sigma_scaling()

        if write:
            mmm = map_model_manager(model=self.model, map_manager=fft_map.as_map_manager())
            box_mmm = mmm.extract_all_maps_around_model(selection_string="all")
            box_mmm.write_map(f"{name}")

        return fft_map

    def map_coefficient_cal(self, data, sigma=None, coefficient="FO", write=False, name="map", map_type=None,
                            mtz_name="name"):

        f_model = mmtbx.f_model.manager(f_obs=miller.array(miller_set=self.miller_set,
                                                           data=data, sigmas=sigma), xray_structure=self.xrs)

        map_coefficients = map_tools.electron_density_map(fmodel=f_model).map_coefficients(
            map_type=coefficient, isotropize=False, fill_missing=False, sharp=False, acentrics_scale=2.0,
            centrics_pre_scale=1.0)
        # map_coefficients = map_coefficients.phase_transfer(flex.double(self.phase_flip(data)))

        fft_map = miller.fft_map(crystal_gridding=self.crystal_gridding(),
                                 fourier_coefficients=map_coefficients.as_non_anomalous_array())
        if map_type == "diff":
            miller_da = miller.array(miller_set=self.miller_set, data=flex.double(data), sigmas=sigma)
            miller_da = miller_da.phase_transfer(f_model.f_model())
            # miller_da = miller_da.map_to_asu().average_bijvoet_mates()

            fft_map = miller_da.fft_map(
                resolution_factor=self.resolution_factor)
            fft_map = fft_map.apply_sigma_scaling()

        else:

            fft_map = fft_map.apply_sigma_scaling()

        if write:
            if map_type == "diff":
                mmm = map_model_manager(model=self.model, map_manager=fft_map.as_map_manager())
                box_mmm = mmm.extract_all_maps_around_model(selection_string="all")
                box_mmm.write_map(f"{name}")

            else:
                fft_4_map_save = map_coefficients.fft_map(crystal_gridding=self.crystal_gridding(),
                                                          resolution_factor=0.25).apply_sigma_scaling()
                iotbx.map_tools.write_ccp4_map(sites_cart=self.xrs.sites_cart(),
                                               unit_cell=fft_4_map_save.unit_cell(),
                                               map_data=fft_4_map_save.real_map(),
                                               n_real=fft_4_map_save.n_real(), buffer=5.0,
                                               file_name=f"{name}")

            if coefficient == "2mFO-DFC" or coefficient == "mFO-DFC":
                column_name_F = "2FEXTRFCWT" if coefficient == "2mFO-DFC" else "FEXTRFCWT" if coefficient == "mFO-DFC" else None
                column_name_P = "PHI2FEXTRFCWT" if coefficient == "2mFO-DFC" else "PHIFEXTRFCWT" if coefficient == "mFO-DFC" else None
                mtz_dataset = miller.array(miller_set=self.miller_set,
                                           data=map_coefficients.as_amplitude_array().data()).as_mtz_dataset(
                    column_root_label=column_name_F, column_types="F")
                mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                                       data=flex.double(
                                                                           map_coefficients.phases(deg=True).data())),
                                             column_root_label=column_name_P, column_types="P")
                mtz_dataset.add_miller_array(
                    miller_array=miller.array(miller_set=self.miller_set,
                                              data=flex.int(self.R_Free.data())), column_root_label="FREE",
                    column_types="I")

                mtz_dataset.mtz_object().write(f"{mtz_name}")

        return fft_map

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
            self.logger.info(
                f"""\nBoth the symmetry does not match so I am negating  the exact Symmetry match option by fixing the symmetry of reference to triggered state.\nThis is done to generate isomorphous difference density map without strict isomorphism.\nIf you have a better option please provide your mtz with no scaling option with Residem.""")

            True

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

    def structure_factor_2_intensity(self, miller_array):
        miller_array_copy = miller_array.deep_copy()
        miller_array_sf = miller_array_copy.data()
        miller_array_sigma = miller_array_copy.sigmas()
        I_miller_array = np.array(miller_array_sf ** 2)
        mask = np.array(miller_array_copy.data() < 0)
        I_miller_array[mask] *= -1
        I_sigmas = np.array(miller_array_sigma * miller_array_copy.data())
        I_from_sf = miller.array(miller_set=self.miller_set, data=flex.double(I_miller_array),
                                 sigmas=flex.double(I_sigmas))
        return I_from_sf

    def Plot_wilson(self, path):
        out_pdf = f'%s/Wilson_plot_comparison.pdf' % path
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)

        resolution_bin = flex.double()
        wilson_plot_bin_dark_b = flex.double()
        wilson_plot_bin_dark_a = flex.double()
        wilson_plot_bin_light_b = flex.double()
        wilson_plot_bin_light_a = flex.double()
        r_iso_4_bin = flex.double()
        cc_iso_4_bin = flex.double()
        bin_iso_4_bin = []

        dark_after = self.structure_factor_2_intensity(self.f_model_dark_scaled)
        light_after = self.structure_factor_2_intensity(self.f_model_light_scaled_obs)
        dark_before = self.structure_factor_2_intensity(self.f_obs_dark_common)
        light_before = self.structure_factor_2_intensity(self.f_obs_light_common)

        dark_before.setup_binner(n_bins=20)
        dark_after.use_binner_of(dark_before)
        light_before.use_binner_of(dark_before)
        light_after.use_binner_of(dark_before)

        for i_bin in dark_before.binner().range_all():
            sel_bin = dark_before.binner().selection(i_bin)
            dark_reso_bin_before = dark_before.select(sel_bin)
            dark_reso_bin_after = dark_after.select(sel_bin)
            light_reso_bin_before = light_before.select(sel_bin)
            light_reso_bin_after = light_after.select(sel_bin)
            if light_reso_bin_before.size() == 0: continue
            bin_res_cent = np.median(dark_before.binner().bin_d_range(i_bin))
            reso_bin = f"%s-%s" % (
                round(dark_before.binner().bin_d_range(i_bin)[0], 2),
                round(dark_before.binner().bin_d_range(i_bin)[1], 2))
            bin_iso_4_bin.append(reso_bin)
            resolution_bin.append(bin_res_cent)
            wilson_plot_bin_dark_b.append(dark_reso_bin_before.wilson_plot())
            wilson_plot_bin_dark_a.append(dark_reso_bin_after.wilson_plot())
            wilson_plot_bin_light_b.append(light_reso_bin_before.wilson_plot())
            wilson_plot_bin_light_a.append(light_reso_bin_after.wilson_plot())

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
        # Plotting on the first subplot

        normalized_wilson_plot_bin_dark_b = (wilson_plot_bin_dark_b[1:] - np.min(wilson_plot_bin_dark_b[1:])) / (
                    np.max(wilson_plot_bin_dark_b[1:]) - np.min(wilson_plot_bin_dark_b[1:]))
        normalized_wilson_plot_bin_light_b = (wilson_plot_bin_light_b[1:] - np.min(wilson_plot_bin_light_b[1:])) / (
                    np.max(wilson_plot_bin_light_b[1:]) - np.min(wilson_plot_bin_light_b[1:]))
        normalized_wilson_plot_bin_dark_a = (wilson_plot_bin_dark_a[1:] - np.min(wilson_plot_bin_dark_a[1:])) / (
                    np.max(wilson_plot_bin_dark_a[1:]) - np.min(wilson_plot_bin_dark_a[1:]))
        normalized_wilson_plot_bin_light_a = (wilson_plot_bin_light_a[1:] - np.min(wilson_plot_bin_light_a[1:])) / (
                    np.max(wilson_plot_bin_light_a[1:]) - np.min(wilson_plot_bin_light_a[1:]))
        ax1.scatter(resolution_bin.as_numpy_array()[1:], normalized_wilson_plot_bin_dark_b, label="Dark-before")
        ax1.scatter(resolution_bin.as_numpy_array()[1:], normalized_wilson_plot_bin_light_b, label="Light-before")
        ax1.set_title('Wilson intensity plot before scaling')
        ax1.set_xlabel('Resolution')
        ax1.set_ylabel('Normalized <I>')

        # # Plotting on the second subplot
        ax2.scatter(resolution_bin.as_numpy_array()[1:], normalized_wilson_plot_bin_dark_a, label="Dark-after")
        ax2.scatter(resolution_bin.as_numpy_array()[1:], normalized_wilson_plot_bin_light_a, label="Light-after")
        ax2.set_title('Wilson intensity plot after scaling')
        ax2.set_xlabel('Resolution')
        ax2.set_ylabel('Normalized <I>')

        xtick = resolution_bin.as_numpy_array()[1:]
        ax1.set_xlim(np.max(xtick) + 0.2, np.min(xtick) - 0.2)
        ax2.set_xlim(np.max(xtick) + 0.2, np.min(xtick) - 0.2)
        ax1.legend()
        ax2.legend()
        plt.tight_layout()
        pdf.savefig(fig)
        pdf.close()

    def Plot_R_iso_CCiso(self, path):
        R_iso_overall, CC_iso_overall = self.overall_R_iso_Cciso()
        self.logger.info(f"{20 * '*'}")
        self.logger.info(f"""Over all R-iso between the two data is {R_iso_overall}""")
        self.logger.info(f"{20 * '*'}")
        self.logger.info(f"""Over all CC-iso between the two data is {CC_iso_overall}""")

        dark = self.f_obs_dark_common.deep_copy()
        light = self.f_obs_light_common.deep_copy()
        resolution_bin = flex.double()
        r_iso_4_bin = flex.double()
        cc_iso_4_bin = flex.double()
        bin_iso_4_bin = []
        dark.setup_binner(n_bins=20)
        light.use_binner_of(dark)
        for i_bin in dark.binner().range_all():
            sel_bin = dark.binner().selection(i_bin)
            dark_reso_bin = dark.select(sel_bin)
            light_reso_bin = light.select(sel_bin)
            if light_reso_bin.size() == 0: continue
            bin_res_cent = np.median(dark.binner().bin_d_range(i_bin))
            reso_bin = f"%s-%s" % (
            round(dark.binner().bin_d_range(i_bin)[0], 2), round(dark.binner().bin_d_range(i_bin)[1], 2))
            bin_iso_4_bin.append(reso_bin)
            resolution_bin.append(bin_res_cent)
            r_work_bin = self.R_iso(light_reso_bin, dark_reso_bin)
            r_iso_4_bin.append(r_work_bin)
            cc_work_bin = dark_reso_bin.correlation(light_reso_bin).coefficient()
            cc_iso_4_bin.append(cc_work_bin)

        df_r_iso_cc_iso = pd.DataFrame()
        df_r_iso_cc_iso["Resolution bin"] = np.array(bin_iso_4_bin)
        df_r_iso_cc_iso["R-iso"] = np.array(r_iso_4_bin)
        df_r_iso_cc_iso["CC-iso"] = np.array(cc_iso_4_bin)
        self.logger.info(df_r_iso_cc_iso)
        out_pdf = f'%s/R_iso_CC_iso.pdf' % path
        self.logger.info(f"{20 * '*'}")
        self.logger.info(f"""The plot of R-iso and CC-iso is saved in {out_pdf}""")
        self.logger.info(f"{20 * '*'}")
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        x = resolution_bin.as_numpy_array()
        fig, ax1 = plt.subplots(figsize=(20, 10))
        sns.lineplot(x=x, y=r_iso_4_bin.as_numpy_array(), ax=ax1, color='b',
                     label=f'R-iso by resolution bin, overall R-iso (%.4f)' % R_iso_overall)
        ax1.set_ylabel('R-iso', color='b', fontsize=24, fontweight='bold')
        ax1.set_xlabel(r'Resolution ($\AA$)', fontsize=24, fontweight='bold')
        ax2 = ax1.twinx()
        sns.lineplot(x=x, y=cc_iso_4_bin.as_numpy_array(), ax=ax2, color='r',
                     label=f'CC-iso by resolution bin, overall CC-iso (%.4f)' % CC_iso_overall)
        ax2.set_ylabel('CC-iso', color='r', fontsize=24, fontweight='bold')
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        all_lines = lines1 + lines2
        all_labels = labels1 + labels2
        ax1.legend(all_lines, all_labels, loc='upper right')
        plt.title('R-iso and CC-iso by high resolution bins', fontsize=24, fontweight='bold')
        ax1.set_xlim(np.max(x[1:]), np.min(x[1:]))
        ax2.get_legend().remove()
        pdf.savefig(fig)
        pdf.close()

    def overall_R_iso_Cciso(self):
        dark = self.f_obs_dark_common.deep_copy()
        light = self.f_obs_light_common.deep_copy()
        R_iso = self.R_iso(dark, light)
        CC_iso = dark.correlation(light).coefficient()

        if R_iso < 0.25:
            pass
        else:
            print("The overall R_iso is above 0.25, Please review the results carefully")
        return R_iso, CC_iso

    @staticmethod
    def R_iso(light, dark):
        Num = np.sum(np.abs(light.data() - dark.data()))
        den = np.sum(np.abs(light.data() + dark.data())) / 2
        R_iso = Num / den if den != 0 else None
        return R_iso

    def phase_flip(self, array):
        arr1 = array.deep_copy()
        arr1 = arr1.as_numpy_array()
        arr2 = self.F_cal_common.phases(deg=True).data().as_numpy_array()
        negative_mask = arr1 < 0
        arr2_flipped = np.copy(arr2)
        arr2_flipped[negative_mask] = arr2_flipped[negative_mask] + 180
        arr2_flipped = np.where(arr2_flipped > 180, arr2_flipped - 360, arr2_flipped)
        arr2_flipped = np.where(arr2_flipped < -180, arr2_flipped + 360, arr2_flipped)
        return arr2_flipped

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
        mtz_dataset.add_miller_array(miller_array=self.f_model_dark.f_calc().as_amplitude_array(),
                                     column_root_label="FC", column_types="F")

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set,
                                      data=flex.double(self.F_cal_common.phases(deg=True).data())),
            column_root_label="PHIC", column_types="P")

        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                               data=flex.double(np.abs(self.delta_F))),
                                     column_root_label="DF_P", column_types="F")
        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                               data=flex.double(np.full(len(self.delta_F), 10))),
                                     column_root_label="DF_10", column_types="F")

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set,
                                      data=flex.double(self.phase)),
            column_root_label="PHIC_Flip", column_types="P")

        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                               data=flex.double(self.delta_F)),
                                     column_root_label="DF", column_types="F")
        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                               data=flex.double(self.sigma_df)),
                                     column_root_label="SIGF_DF", column_types="Q")
        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.wDF)),
            column_root_label='DF_Ren', column_types="F")
        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                               data=flex.double(self.sigma_ren_df)),
                                     column_root_label="SIGF_Ren", column_types="Q")

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.rDF)),
            column_root_label='DF_Ursby', column_types="F")
        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                               data=flex.double(self.sigma_ursby_df)),
                                     column_root_label="SIGF_Ursby", column_types="Q")

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.hDF)),
            column_root_label='DF_hekstra', column_types="F")
        mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                               data=flex.double(self.sigma_hekstra_df)),
                                     column_root_label="SIGF_hekstra", column_types="Q")

        if self.alpha != 1:
            mtz_dataset.add_miller_array(
                miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.hDF_1)),
                column_root_label='DF_hekstra_A_1', column_types="F")
            mtz_dataset.add_miller_array(miller_array=miller.array(miller_set=self.miller_set,
                                                                   data=flex.double(self.sigma_hekstra_df_1)),
                                         column_root_label="SIGF_hekstra_A_1", column_types="Q")

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.ursby_factor)),
            column_root_label='ursby_w', column_types='W')
        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.ren_factor)),
            column_root_label='ren_w', column_types='W')
        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.hekstra_factor)),
            column_root_label='hekstra_w', column_types='W')
        if self.alpha != 1:
            mtz_dataset.add_miller_array(
                miller_array=miller.array(miller_set=self.miller_set, data=flex.double(self.hekstra_factor_1)),
                column_root_label='schmidt_w', column_types='W')

        mtz_dataset.add_miller_array(
            miller_array=miller.array(miller_set=self.miller_set, data=flex.int(self.R_Free.data())),
            column_root_label='FREE', column_types='I')
        mtz_dataset.mtz_object().write(name)

        return mtz_dataset

    def Plot_weight(self, path):

        out_pdf = os.path.join(path, "Difference_structure_factor_distribution.pdf")
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)

        # List of dictionaries containing data to plot
        list_to_plot = [
            {"x": self.delta_F.as_numpy_array(), "y": np.abs(self.delta_F.as_numpy_array()) / np.array(self.sigma_df),
             "z": np.array(self.sigma_df), "label": "No weighting scheme"},
            {"x": self.wDF.as_numpy_array(), "y": np.abs(self.wDF.as_numpy_array()) / np.array(self.sigma_ren_df),
             "z": np.array(self.sigma_ren_df),
             "label": "Weighting scheme by Ren"},
            {"x": self.rDF.as_numpy_array(), "y": np.abs(self.rDF.as_numpy_array()) / np.array(self.sigma_ursby_df),
             "z": np.array(self.sigma_ursby_df),
             "label": "Weighting scheme by Ursby"},
            {"x": self.hDF.as_numpy_array(), "y": np.abs(self.hDF.as_numpy_array()) / np.array(self.sigma_hekstra_df),
             "z": np.array(self.sigma_hekstra_df),
             "label": f"Weighting scheme by Hekstra, alpha=%s" % self.alpha}
        ]
        if self.alpha != 1:
            list_to_plot.append({"x": self.hDF_1.as_numpy_array(),
                                 "y": np.abs(self.hDF_1.as_numpy_array()) / np.array(self.sigma_hekstra_df_1),
                                 "z": np.array(self.sigma_hekstra_df_1),
                                 "label": "Weighting scheme by Schmidt, alpha=1"})

        for item in list_to_plot:
            fig, ax = plt.subplots(figsize=(10, 7))
            percentage_above_threshold = np.sum(item["y"] > 10) / len(item["y"]) * 100

            Y_Tick_label = r"$\frac{\left| \Delta F \right|}{\sigma_{|\Delta F|}}$"

            filtered_indices = (item["y"] > -5) & (item["y"] < 5)
            x_filtered = item["x"][filtered_indices]
            y_filtered = item["y"][filtered_indices]

            sns.kdeplot(x=x_filtered, y=y_filtered, cmap="Blues", fill=True, bw_adjust=0.5, ax=ax, label=item["label"])
            ax.set_xlabel(r"$\left| \Delta F \right|$", size=20)
            ax.set_ylabel(Y_Tick_label, size=24)
            ax.set_title('Distribution of difference Fourier coefficients')
            ax.legend()

            ax_inset = ax.inset_axes([0.7, 0.7, 0.25, 0.25])  # Create inset at a specific location and size

            sns.kdeplot(x=item["x"], y=item["z"], cmap="Reds", fill=True, bw_adjust=0.5, ax=ax_inset)

            ax_inset.set_xlabel(r"$\left| \Delta F \right|$", fontsize=10)
            ax_inset.set_ylabel(r'$\sigma_{|\Delta F |}$', fontsize=10)
            label_text = f"Percentage of {Y_Tick_label} Points above {10} : {percentage_above_threshold:.2f}%"
            ax.legend([label_text, item["label"]], loc='upper left')
            ax.set_ylim(-2, 5)

            pdf.savefig(fig)
            plt.close()

        # Close the PDF
        pdf.close()

    def W_weight(self, weight="ren", alpha=1):
        """ weight as suggested by Ren 2001 and hekstra 2016 """
        SigDF = np.sqrt(self.f_model_dark_scaled.sigmas() ** 2 + self.f_model_light_scaled_obs.sigmas() ** 2)
        if weight == "ren":
            weight_re = flex.double(self.compute_target_weight_k(self.delta_F.as_numpy_array(), SigDF))
        elif weight == "hekstra":
            # weight_re = flex.double(self.compute_hekstra_weight(alpha).data())
            weight_re = flex.double(self.compute_target_weight_h(self.delta_F.as_numpy_array(), SigDF, alpha))

        return weight_re

    @staticmethod
    def compute_target_weight_k(F, SIGF):
        ''' This is a simplified weighting scheme implemented by Ren 2001 '''

        w = np.power((1 + (SIGF ** 2 / np.mean(np.array(SIGF)) ** 2) + (F ** 2 / (F).mean() ** 2)), -1)
        return w  # /np.mean(w)

    @staticmethod
    def compute_target_weight_h(F, SIGF, alpha):
        ''' This is a simplified weighting scheme implemented by Ren 2001 '''

        w = np.power((1 + (SIGF ** 2 / np.mean(np.array(SIGF ** 2))) + alpha * (F ** 2 / (F ** 2).mean())), -1)
        return w  # /np.mean(w)

    def compute_ursby_weight(self):
        """ This is the weighting scheme implemented by ursby 1997 also implemented in xtrapol8 2022 
        In developing this function, we referred to the following work, licensed under the MIT License:
        * De Zitter, E., Coquelle, N., Oeser, P., Barends, T. & Colletier, J.-P. 
        Xtrapol8 enables automatic elucidation of low-occupancy intermediate states in crystallographic studies. 
        Communications Biology 5 (2022). DOI: 10.1038/s42003-022-03575-7. 
        * Xtrapol8, Copyright (c) 2021 Elke De Zitter, Nicolas Coquelle, Thomas Barends and Jacques-Philippe Colletier 
        https://github.com/ElkeDeZitter/Xtrapol8/
        The original license text can be found in the LICENSE file.

        """

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
                                data=weight / flex.mean(flex.double(weight)))

        if q_weight.indices().size() != self.indices.size():
            raise ValueError("The size of miller not matching the q weight")

        else:
            return q_weight


class Prescaling:
    def __init__(self, pdb, miller_arrays):
        self.pdb_in = iotbx.pdb.input(file_name=pdb)
        self.miller_array = miller_arrays
        self.pdb_hierarchy = mmtbx.model.manager(model_input=self.pdb_in).get_hierarchy()
        self.n_residues = self.pdb_hierarchy.overall_counts().n_residues
        self.Anios_corrected = self.AnisoCorrect()
        self.u_star_correct_nat = self.Anios_corrected.u_star
        self.miller_array_scaled = absolute_scaling.anisotropic_correction(self.miller_array,
                                                                           self.Anios_corrected.p_scale,
                                                                           self.u_star_correct_nat)

    def matthews_analyses(self):
        matthews_analyses = matthews.matthews_rupp(
            crystal_symmetry=self.miller_array,
            n_residues=self.n_residues)
        return matthews_analyses

    def AnisoCorrect(self):
        aniso_correct = absolute_scaling.ml_aniso_absolute_scaling(
            miller_array=self.miller_array,
            n_residues=self.n_residues * self.miller_array.crystal_symmetry().space_group().order_z() * self.matthews_analyses().n_copies)
        return aniso_correct
