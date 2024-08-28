from iotbx.reflection_file_reader import any_reflection_file
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
from cctbx import miller, maptbx
from cctbx.array_family import flex
from mmtbx import map_tools


class MTZ2MAP:
    def __init__(self, pdb, mtz, label=None, map_type="2mFo-DFc", name="output"):
        self.pdb_in = iotbx.pdb.input(pdb)
        self.model = mmtbx.model.manager(model_input=self.pdb_in)

        self.map_type = map_type
        self.mtz = any_reflection_file(mtz)
        self.hkl = self.mtz.as_miller_arrays()
        self.miller_dict = {ma.info().label_string(): ma for ma in self.hkl}
        self.keys = [k for k, v in self.miller_dict.items()]
        self.low_resolution, self.high_resolution = self.miller_dict[self.keys[0]].d_max_min()
        self.indices = self.miller_dict[self.keys[0]].indices()
        self.crystal_symmetry = self.miller_dict[self.keys[0]].crystal_symmetry()
        self.miller_set = miller.set(self.crystal_symmetry, self.indices, anomalous_flag=False)
        self.label = label
        self.xrs = self.model.get_xray_structure()
        self.resolution_factor = 0.5 / self.high_resolution
        self.array = self.miller_array_to_convert()
        self.fft_map = self.write_ccp4_map(name)

    def miller_array_to_convert(self):
        try:
            array = self.miller_dict[self.label]

        except KeyError:
            print(f"Please choose from the available label for the file \n {self.mtz}:{list(self.miller_dict.keys())}")
            exit()
        return array


    def write_ccp4_map(self, name):
        # f_model = mmtbx.f_model.manager(f_obs=miller.array(miller_set=self.miller_set,
        #                                                    data=self.array.data(), sigmas=self.array.sigmas()),
        #                                 xray_structure=self.xrs)
        f_model = mmtbx.f_model.manager(f_obs=self.array,
                                        xray_structure=self.xrs)
        f_model.update_all_scales()


        map_coefficients = map_tools.electron_density_map(fmodel=f_model).map_coefficients(
            map_type=self.map_type, isotropize= True, fill_missing=False, sharp=False,exclude_free_r_reflections=True,acentrics_scale=2.0, centrics_pre_scale=1.0)

        map_coefficients = map_coefficients.average_bijvoet_mates()

        fft_map = map_coefficients.fft_map(d_min=self.high_resolution)
        # fft_map = miller.fft_map(crystal_gridding=self.crystal_gridding(),
        #                          fourier_coefficients=map_coefficients.as_non_anomalous_array())

        fft_map.as_ccp4_map(f"{name}_{self.map_type}.ccp4")

        return fft_map

    def crystal_gridding(self):
        """crystal grid parameter"""
        resolution_factor = 0.5 / self.high_resolution
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


def run(pdb, mtz, label, map_type, name):
    MTZ2MAP(pdb, mtz, label, map_type, name)


def run_all():
    import argparse
    from datetime import datetime

    start = datetime.now()
    parser = argparse.ArgumentParser(prog='residem_mtz2map', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f""" This script generates ccp4 map for visualization in pymol. 
        Few example to run the code: \n\n
        1. residem_mtz2map -r 5b6v.pdb -m 5b6v.mtz -l 'FP,SIGFP' -t 2mFo-Dfc/mFo-DFc/2Fo-Fc/mFo-Fc/mFo/Fo -o output \n\n
        2. residem_mtz2map -r 5b6w.pdb -m 5b6w.mtz -l 'FP,SIGFP' -t 2mFo-DFc -o 5b6w
    """)

    parser.add_argument("-r", "--ref", type=str, help="reference pdb file")
    parser.add_argument("-m", "--mtz", type=str, help="mtz file for which map need to be transformed")
    parser.add_argument("-l", "--label", type=str, help="label")

    parser.add_argument("-t", "--map_type", type=str, help="map_cofficient for which map has to be generated")
    parser.add_argument("-o", "--output", type=str, help="output file name")

    args = parser.parse_args()
    map_type = args.map_type if args.map_type is not None else "2mFo-DFc"
    name = args.output if args.output is not None else "output"
    run(args.ref, args.mtz, args.label, map_type, name)

    print(datetime.now() - start)


if __name__ == '__main__':
    run_all()
