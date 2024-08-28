from iotbx.reflection_file_reader import any_reflection_file
import iotbx.map_tools
import os
from cctbx import miller, maptbx
from cctbx.maptbx import crystal_gridding
from cctbx import maptbx


class mtz2msp:
    def __init__(self, mtz, high_resolution, low_resolution, grid_step=0.5):
        self.mtz = any_reflection_file(os.path.realpath(mtz))
        self.hkl = self.mtz.as_miller_arrays()
        self.hkl_dict = {ma.info().label_string(): ma for ma in self.hkl}
        self.phic = self.hkl_dict["PHIC"]
        self.DF = self.hkl_dict["DF"]
        self.ren = self.hkl_dict['kDF']
        self.ursby = self.hkl_dict['qDF']
        self.low_resolution = low_resolution
        self.high_resolution = high_resolution
        self.grid_step = grid_step
        self.resolution_factor = self.grid_step / self.high_resolution
        # self.map_coff_to_map(self.DF, name="Data_folder/F_obs_minus_F_obs.ccp4")
        # self.map_coff_to_map(self.ren, name="Data_folder/F_obs_minus_F_obs_ren_weight.ccp4")
        # self.map_coff_to_map(self.ursby, name="Data_folder/F_obs_minus_F_obs_ursby_weight.ccp4")

    def map_coff_to_map(self, miller_array, name):
        map_coeffs = iotbx.map_tools.combine_f_phi_and_fom(f=miller_array, phi=self.phic, fom=None)
        map_coeffs = map_coeffs.map_to_asu().average_bijvoet_mates()
        map_coeffs = map_coeffs.resolution_filter(d_min=self.high_resolution,
                                                  d_max=self.low_resolution)

        # map_coeffs.fft_map(resolution_factor=self.resolution_factor).apply_sigma_scaling().as_ccp4_map(name)

        cr_gr = maptbx.crystal_gridding(
            unit_cell=map_coeffs.unit_cell(),
            space_group_info=map_coeffs.space_group_info(), pre_determined_n_real=map_coeffs.data().all())

        map_coeffs.fft_map(resolution_factor=self.resolution_factor,
                           symmetry_flags=maptbx.use_space_group_symmetry,
                           crystal_gridding=cr_gr).apply_sigma_scaling().as_ccp4_map(name)
