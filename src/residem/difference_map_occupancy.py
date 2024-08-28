from collections import OrderedDict
from cctbx import miller, maptbx
from iotbx import ccp4_map, mtz
import mmtbx.f_model
import iotbx.pdb
import numpy as np
from iotbx import reflection_file_reader
import math
from cctbx.array_family import flex
import mmtbx.utils
from residem.structUtils import ModelSelection
import pandas as pd
from collections import defaultdict
from residem.multiprocess_func import make_parallel
import re
import json
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import math
from residem.scan_all_points import Peaks

import cctbx
from mmtbx import map_tools

matplotlib.use('Agg')
font = {'weight': 'bold', 'size': 8}

matplotlib.rc('font', **font)




class Difference_map:

    def __init__(self, pdb, Light_dark_scaled, atom_list=None, grid_step=0.5, msd=0.35, cluster=None, weight = None):
        # input parameters
        self.pdb = iotbx.pdb.input(file_name=pdb)
        self.xrs = self.pdb.xray_structure_simple()
        self.msd = msd
        self.grid_step = grid_step
        self.Light_dark_scaled = Light_dark_scaled
        self.cluster = "cluster" if cluster is None else cluster if cluster.lower() in ["map_atom", "all",
                                                                                        "cluster"] else "cluster "
        self.atom_list = ModelSelection(pdb).atom_df if atom_list is None else self.atom_list_df(atom_list)
        self.weight = None if weight is None else weight.lower() if weight.lower() in ["ren","ursby"] else None
        # miller process
        self.reflection_file = reflection_file_reader.any_reflection_file(self.Light_dark_scaled)
        self.miller_arrays = self.reflection_file.as_miller_arrays()
        self.miller_dict = {ma.info().label_string(): ma for ma in self.miller_arrays}
        self.high_resolution = round(
            max(self.miller_dict["F_DARK,SIGF_DARK"].d_min(), self.miller_dict['F_LIGHT,SIGF_LIGHT'].d_min()), 2)
        self.cr_gd = self.crystal_gridding()

        # common miller indices
        self.F_dark_miller = self.miller_dict['F_DARK,SIGF_DARK'].common_set(self.miller_dict['F_LIGHT,SIGF_LIGHT']).deep_copy()
        self.F_light_miller = self.miller_dict['F_LIGHT,SIGF_LIGHT'].common_set(self.miller_dict['F_DARK,SIGF_DARK']).deep_copy()
        self.F_cal_common = self.miller_dict['F_CAL'].common_set(self.F_dark_miller).as_amplitude_array().deep_copy()
        self.PHIC_common = self.miller_dict['DF,PHIC'].common_set(self.F_dark_miller).phases(deg=True).deep_copy()
        self.ren_F = self.miller_dict['kDF'].common_set(self.miller_dict['F_DARK,SIGF_DARK']).as_amplitude_array().data()
        self.ursby_F = self.miller_dict['qDF'].common_set(self.miller_dict['F_DARK,SIGF_DARK']).as_amplitude_array().data()
        self.DF = self.miller_dict['DF,PHIC'].common_set(self.miller_dict['F_DARK,SIGF_DARK']).as_amplitude_array().data()
        self.DF_cal = self.miller_dict['DF_cal'].common_set(self.miller_dict['F_DARK,SIGF_DARK']).as_amplitude_array().data()

        ## cell symmetry
        self.indices = self.F_dark_miller.indices()
        self.cs = self.F_dark_miller.crystal_symmetry()
        self.uc = self.cs.unit_cell()
        self.miller_set = miller.set(self.cs, self.indices, anomalous_flag=False)


        ### miller to map_coefficient :





        # f_model dark and light
        f_model_dark = self.f_model(self.F_dark_miller)
        f_model_light = self.f_model(self.F_light_miller)


        # #map_coefficient calculation for occupancy estimation.
        ##There are nine type of map coefficient which has to be estimated based on user inputs.
        ##creating a map coefficient is not time and memory intensive only occupancy estimation is so, dictionary is
        ##created for map coefficient

        #self.map_coefficient_dict_extrapolated ={"", }

        ##occupancy estimation
        ## lowest negative peak sum in differnce elecrton density map for that position in extrapolated map.
        ## first position of voxel in the
        ###difference electron density map is calculated using map coefficient

        # plotting
        # self.d_sigma = defaultdict(int)
        # self.d_vol = defaultdict(int)

    ######################  methods ###################################

    def atom_list_df(self, atom_list):
        """Reads the clustered data as inputfile for choosing the residues"""
        atom_list_df = pd.read_csv(atom_list)
        if set(["x", "y", "z", "resn"]) <= set([x.lower() for x in list(atom_list_df.columns)]):
            atom_list_df = atom_list_df[["resn", 'resi', 'atom name', 'X', 'Y', 'Z']]
            atom_list_df.set_axis(['resnum', 'resname', 'atomname', 'X', 'Y', 'Z'], axis=1, inplace=True)

        else:
            raise ImportError(f"The given dataframe is not vaild, it should contain X,Y,Z co-ordinates \n "
                              "use csv from cluster analysis or map analysis ")
        return atom_list_df



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

    def f_model(self, miller):
        """ creates f_model with given miller arrays"""
        f_model = mmtbx.f_model.manager(f_obs=miller, xray_structure=self.xrs)

        # f_model_dark.update_all_scales(fast=True, remove_outliers=False)
        k = f_model.k_isotropic() * f_model.k_anisotropic()
        fd = f_model.f_obs().customized_copy(data=f_model.f_obs().data() / k)
        fd = fd.phase_transfer(phase_source=f_model.f_model())
        return fd

    ### weight ##



    ##### map coefficient ######
    def map_coefficient(self, data, name):
        f_model = mmtbx.f_model.manager(f_obs=miller.array(miller_set=self.miller_set,
                                                           data=data), xray_structure=self.xrs)
        map_coefficients = map_tools.electron_density_map(fmodel=f_model).map_coefficients(
            map_type="Fo", isotropize=True, fill_missing=False)
        fft_map = miller.fft_map(crystal_gridding=self.crystal_gridding(), fourier_coefficients=map_coefficients)
        fft_map = fft_map.apply_sigma_scaling()
        fft_map.as_ccp4_map(name)
        return fft_map


    def fft_map(self, miller_object_array, input_type):
        f_model = mmtbx.f_model.manager(
            f_obs=miller_object_array,
            xray_structure=self.xrs)

        # map_coefficients = map_tools.electron_density_map(fmodel=f_model).map_coefficients(
        #     map_type="Fo", isotropize=True, fill_missing=False
        # )
        # fft_map = cctbx.miller.fft_map(
        #     crystal_gridding=map_model.crystal_gridding(), fourier_coefficients=map_coefficients
        # )
        # fft_map = fft_map.apply_sigma_scaling()
        # fft_map.as_ccp4_map("test_Fo.ccp4")

        if input_type != "diff":
            # f_model.update_all_scales(fast=True, remove_outliers=False)
            # scale the model to observed
            k = f_model.k_isotropic() * f_model.k_anisotropic()
        elif input_type == "diff":
            k = 1
        fo = f_model.f_obs().customized_copy(data=f_model.f_obs().data() / k)
        fo = fo.phase_transfer(phase_source=f_model.f_model())
        fc = f_model.f_calc().customized_copy(data=f_model.f_calc().data())

        if input_type == "diff" or input_type == "extra":
            miller_array_diff_extra = miller.array(miller_set=self.miller_set,
                                                   data=fo.data())
            fft_diff_extra = self.fft_map_c(miller_array_diff_extra)

            return fft_diff_extra

        elif input_type == "extra_diff_cal":
            miller_array_diff = miller.array(miller_set=self.miller_set,
                                             data=fo.data() - fc.data())
            fft_diff = self.fft_map_c(miller_array_diff)
            return fft_diff

        elif input_type == "extra_diff_dark":
            miller_array_dark = miller.array(miller_set=self.miller_set,
                                             data=fo.data() - self.f_model_dark().data())

            fft_dark = self.fft_map_c(miller_array_dark)

            return fft_dark

    def fft_map_c(self, miller_array_data):
        fft_map = miller.fft_map(
            crystal_gridding=self.cr_gd,
            fourier_coefficients=miller_array_data.average_bijvoet_mates())
        # f_000=mmtbx.utils.f_000(self.xrs, self.msd).f_000)

        fft_map = fft_map.apply_sigma_scaling()

        return fft_map

    ### plot ####
    @staticmethod
    def plot_weight(DF, SIGDF, weight):
        fig, ax = plt.subplots(figsize=(10, 7))
        pts = ax.scatter(DF, abs(DF) / SIGDF, c=weight)
        ax.set_xlabel(r"$\Delta F$")
        ax.set_ylabel(r"$\frac{\left| \Delta F \right|}{\sigma_{\Delta F}}$")

        # Inset
        axins = ax.inset_axes([0.6, 0.6, 0.37, 0.37])
        axins.scatter(DF, abs(DF) / SIGDF, c=weight)
        x1, x2, y1, y2 = -25, 25, 0, 10
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        ax.indicate_inset_zoom(axins)

        fig.colorbar(pts, label="Weight")
        plt.savefig("weight.png", dpi=600)
