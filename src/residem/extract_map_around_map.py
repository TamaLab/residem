import mmtbx
import mmtbx.model
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager

import mmtbx.utils
import iotbx.pdb
import os
from iotbx.data_manager import DataManager



class MapOutput:
    def __init__(self,pdb,mapin,output="output"):
        self.pdb = pdb 
        self.mapin = mapin
        self.output = output

        self.pdb_in = iotbx.pdb.input(self.pdb) 
        self.map_manager =  map_manager(self.mapin)
        self.model = mmtbx.model.manager(model_input=self.pdb_in)
        self.map_model = map_model_manager(map_manager=self.map_manager, model=self.model)

        self.dm = DataManager()
        self.dm.set_overwrite(True)

        self.box_mmm = self.map_model.extract_all_maps_around_model(  # extract a box around model
            selection_string="all")

        self.dm.write_real_map_file(self.box_mmm.map_manager(),
                            filename= f"{self.output}.ccp4")



def run(pdb, mapin,output ):
    sc = MapOutput(pdb=pdb, mapin=mapin, output=output)


def run_2_all():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--pdb", help="pdb_file")
    parser.add_argument("-m","--map", help="ccp4 file name")
    parser.add_argument("-o","--output",help="output_file_name")

    args = parser.parse_args()
    output = "output" if args.output is not None else args.output

    run(args.pdb, args.map,output)


if __name__ == "__main__":
    run_2_all()