from residem.master import DataInput
import mmtbx
import mmtbx.model
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager

import mmtbx.utils
import iotbx.pdb
import os


class path:
    def __init__(self, file):
        self.file = file
        self.stem = self.get_extension()[0]
        self.suffix = self.get_extension()[1]

    def get_extension(self):
        basename = os.path.basename(self.file)
        name_list = basename.split(".")
        stem = name_list[0]
        suffix = "." + name_list[1]
        return stem, suffix


class CctbxFunc(DataInput):
    def __init__(self, pdb, map_file, high_resolution, sigma, msd=0.35):
        self.pdb = pdb
        self.pdb_in = iotbx.pdb.input(pdb)
        self.map = map_manager(map_file)
        self.msd = msd
        self.high_resolution = high_resolution
        self.model = mmtbx.model.manager(model_input=self.pdb_in)
        self.map_model = map_model_manager(map_manager=self.map, model=self.model)
        self.map_model.box_all_maps_around_density_and_shift_origin()
        self.out_put_map = path(map_file).stem
        self.f000 = self.f000_model()
        self.map_coefficient = self.map_model.map_as_fourier_coefficients(d_min=self.high_resolution)
        self.real_map = self.map_coefficient.fft_map()
        self.sigma = sigma

    def map_write(self):
        self.map_model.box_all_maps_around_model_and_shift_origin()
        self.map_model.write_map(self.out_put_map + '_around_protein.ccp4')

    def f000_model(self):
        """calculate f000 from the model"""
        xrs = self.pdb_in.xray_structure_simple()
        f_000 = mmtbx.utils.f_000(xrs, self.msd).f_000
        return f_000


def run(pdb, map_file, high_resolution, sigma):
    map_box_file = CctbxFunc(pdb, map_file, high_resolution, sigma)
    map_box_file.map_write()


def main():
    import argparse
    from datetime import datetime
    start = datetime.now()
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file", help="pdb_file")
    parser.add_argument("map_file", help="Fo-Fo_map_file in CCP4 format")
    parser.add_argument("d_min", help="resolution cutoff")
    parser.add_argument("sigma", help="sigma cutoff")
    args = parser.parse_args()
    run(args.pdb_file, args.map_file, float(args.d_min), float(args.sigma))
    print(datetime.now() - start)


if __name__ == "__main__":
    main()
