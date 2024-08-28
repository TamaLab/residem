import re
from scitbx.array_family import flex
import itertools
import pathlib
import numpy as np
import iotbx.pdb
import mmtbx
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager
from residem.multiprocess_func import make_parallel
import tqdm
import os


class MapCofSet:
    """ set map value to zero based on non selection, selection value as similar to phenix selection syntax:
    https://phenix-online.org/documentation/reference/atom_selections.html
    resseq 182:214, all, name CA
    """

    def __init__(self, pdb, Map, selection="all", radius=5.0, box_cushion=5.0, write=True,output="map"):
        self.pdb_in = iotbx.pdb.input(pdb)
        self.model = mmtbx.model.manager(model_input=self.pdb_in)
        self.mapp = map_manager(Map)
        self.write = write
        self.output = output
        self.path = os.getcwd()
        self.map_name = pathlib.PurePath(Map).stem
        self.selection = "all" if selection is "all" else selection
        self.map_model = map_model_manager(map_manager=self.mapp, model=self.model)
        self.map_model_grid_size = self.map_model.map_manager().data.all()
        self.map_cof_2 = self.map_model.as_map_model_manager().deep_copy()
        self.map_cof_2.box_all_maps_around_model_and_shift_origin(selection_string=self.selection,
                                                                  soft_mask_radius=radius, box_cushion=box_cushion)
        self.map_cof_2_orgin_shift_in_grid_units = self.map_cof_2.map_manager().origin_shift_grid_units

        self.map_grid_shift = np.array(self.map_model_grid_size,dtype=np.float32) + np.array(self.map_cof_2_orgin_shift_in_grid_units,dtype=np.float32)

        self.map_model_grid_size = self.map_model.map_manager().data.all()
        self.map_cof_2_grid_size = self.map_cof_2.map_manager().data.all()
        self.map_grid_shift = np.array(self.map_model_grid_size,dtype=np.float32) + np.array(self.map_cof_2_orgin_shift_in_grid_units,dtype=np.float32)
        self.diff_map_length = self.map_grid_out()
        # self.diff_map_list_2_zero = self.map_model_grid_2_zero()

        if self.write:
            # self.diff_map_list_2_zero.write_map(f"%s/%s_Zeroed.ccp4" % (self.path, self.map_name))
            self.map_cof_2.write_map(f"%s/%s_4_SVD.ccp4" % (self.path, self.output))

        ### set map_data to zero for regions except specified

    def map_grid_out(self):
        map_cof_2_grid = list(itertools.product(range(self.map_cof_2_grid_size[0]),
                                                range(self.map_cof_2_grid_size[1]),
                                                range(self.map_cof_2_grid_size[2])))

        assert len(map_cof_2_grid) == self.map_cof_2.map_data().size()

        map_model_2_grid = list(itertools.product(range(self.map_model_grid_size[0]),
                                                  range(self.map_model_grid_size[1]),
                                                  range(self.map_model_grid_size[2])))

        assert len(map_model_2_grid) == self.map_model.map_data().size()

        self.map_cof_2_grid = map_cof_2_grid
        self.map_model_2_grid = map_model_2_grid
        return len(map_model_2_grid) - len(map_cof_2_grid)

    def map_model_not(self, list_item):
        i = list_item
        x = i[0] + self.map_grid_shift[0] \
            if (i[0] + self.map_grid_shift[0]) < self.map_model_grid_size[0] \
            else i[0] + self.map_grid_shift[0] - self.map_model_grid_size[0]

        y = i[1] + self.map_grid_shift[1] \
            if (i[1] + self.map_grid_shift[1]) < self.map_model_grid_size[1] \
            else i[1] + self.map_grid_shift[1] - self.map_model_grid_size[1]

        z = i[2] + self.map_grid_shift[2] \
            if (i[2] + self.map_grid_shift[2]) < self.map_model_grid_size[2] \
            else i[2] + self.map_grid_shift[2] - self.map_model_grid_size[2]

        return (x, y, z)

    def differnce_list(self):
        print("Difference list calculation")
        result = make_parallel(self.map_model_not)(tqdm.tqdm(self.map_cof_2_grid))
        return result

    def map_model_grid_2_zero(self):
        print(
            f"Setting the value of map data to zero for the regions other than selected residues(%s) " % self.selection)
        map_object = self.map_model.deep_copy()

        def func(i):
            map_object.map_data()[i] = 0

        difference_grid = list(set(self.map_model_2_grid).difference(set(self.differnce_list())))
        make_parallel(func)(tqdm.tqdm(difference_grid))
        return map_object


def run(pdb, Map, residue, radius, box,output):
    sc = MapCofSet(pdb=pdb, Map=Map, selection=residue, radius=radius, box_cussion=box,output=output)


def run_2_all():
    import argparse

    parser = argparse.ArgumentParser(prog='residem_extract_around_model', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f""" This is script set extract the map around the specified region and sets zero to other region which can be futher used for SVD analysis.
    selection value as similar to phenix selection syntax:
    https://phenix-online.org/documentation/reference/atom_selections.html
    resseq 182:214, all, name CA
    Example: \n\n
    $ residem_extract_around_model -r 5b6v -m F_obs_minus_F_obs.ccp4 -s "resseq 182:214" 
    $ residem_extract_around_model -r 5b6v -m F_obs_minus_F_obs.ccp4 -s "resseq 182:214" -d 5.0 -b 5.0
    
    """)
    parser.add_argument("-r","--pdb", help=" reference pdb_file")
    parser.add_argument("-m","--map", help="Fo-Fo_map_file in CCP4 format")
    parser.add_argument("-s","--selection",nargs='?', action='store', type=str,help="residue selection")
    parser.add_argument("-d","--radius", nargs='?', action='store', type=float, help="radius")
    parser.add_argument("-b","--box", nargs='?', action='store', type=float, help="box cushion")
    parser.add_argument("-o", "--output", nargs='?', action='store', type=str, help="output name")
    args = parser.parse_args()
    residue = "all" if args.radius is None else args.selection
    radius = 5.0 if args.radius is None else args.radius
    box = 5.0 if args.box is None else args.box
    output = "map" if args.output is None else args.output
    run(args.pdb, args.map, residue, radius, box,output)

if __name__ == "__main__":
    run_2_all()
