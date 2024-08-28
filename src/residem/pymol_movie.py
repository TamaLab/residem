import os

import pymol
from pymol import cmd
import numpy as np
from glob import glob
from residem.pymol_image import PymolImage
from residem.multiprocess_func import make_parallel



class Pymolmovie(PymolImage):
    """A program to create a moview of extrapolated map and extrapolated differnce map superimposed to isomorphous difference map"""
    def __init__(self,pdb,ccp4_map,ccp4_map_dir,map_type,selection="all",sigma=3.0,map_sigma=1.5,carve=2.0,positve=True,negative=True,all_map=True,name="output",rot_x=0,rot_y=0,rot_z=0,path=os.getcwd()):
        super().__init__(pdb,ccp4_map,resilist=selection,sigma=sigma,carve=carve,positve=positve,negative=negative,all_map=all_map,name=name,rot_x=rot_x,rot_y=rot_y,rot_z=rot_z)
        self.pdb = cmd.load(pdb,"pdb")
        self.difference_map = cmd.load(ccp4_map,"ccp4_diff_map")

        self.ccp4_map = glob(ccp4_map_dir) # a directory containing ccp4 map files
        self.map_type = map_type if map_type in ["extrapolated","extrapolated_diff"] else "extrapolated_diff"
        self.sorted_ccp4_map = self.sort_ccp4_map()
        self.name = name
        self.path = path
        self.map_sigma = map_sigma
        self.generate_image_default()
        self.generate_difference_map()

        # Three things to do here, one generate the image of the difference map which is singel, and extrapolated map, and the difference map are many
        # generate the image of the difference map what has to be done ther


    def sort_ccp4_map(self):
        occupancy = [float(file.split("_")[-1].split(".")[0]) for file in self.ccp4_map]
        ccp4_map = [file for _, file in sorted(zip(occupancy, self.ccp4_map))]
        return ccp4_map

    def generate_difference_map(self):
        os.makedirs(f"{self.path}/difference_map",exist_ok=True)
        if self.all_map :
            self.map_all()

        elif not self.all_map and self.positive :
            self.map_p()

        elif not self.all_map and self.negative :
            self.map_n()
        cmd.center("pdb and site", 0, 1)
        cmd.zoom("pdb and site", "5.0")
        # cmd.clip("slab",10,"pdb and site")
        cmd.bg_color("white")
        cmd.show("nb_spheres", "pdb and site")
        cmd.show("cartoon","all")
        cmd.ray("1024", "1024")
        cmd.png(f"%s/%s_%s.png" % (f"{self.path}/difference_map",self.name, "cartoon"))


    def map_all(self):
        cmd.do(f"isomesh map_p_diff,ccp4_diff_map,%s,pdb and {self.resilist},carve=%s" % (self.sigma, self.carve))
        cmd.do(f"isomesh map_n,ccp4_diff_map,-%s,pdb and {self.resilist},carve=%s" % (self.sigma, self.carve))
        cmd.show("mesh", "map_n_diff")
        cmd.show("mesh", "map_p_diff")
        cmd.color("green", "map_p_diff")
        cmd.color("red", "map_n_diff")

    def map_p(self):
        cmd.do(f"isomesh map_p_diff,ccp4_diff_map,%s, pdb and {self.resilist},carve=%s" % (self.sigma, self.carve))
        cmd.show("mesh", "map_p_diff")
        cmd.color("green", "map_p_diff")

    def map_n(self):
        cmd.do(f"isomesh map_n_diff,ccp4_diff_map,-%s,pdb and {self.resilist},carve=%s" % (self.sigma, self.carve))
        cmd.show("mesh", "map_n_diff")
        cmd.color("red", "map_n_diff")


    def generate_image_map_xtra(self):
        os.makedirs(f"{self.path}/extrapolated_map",exist_ok=True)

        def generate_image_extra(x):
            cmd.load(x,"ccp4_map")
            self.name = f"%s_extra.png"%x.split(".")[0]
            cmd.do(f"isomesh map_extra,ccp4_map,%s, pdb and {self.resilist},carve=" % (self.map_sigma, self.carve))
            cmd.show("mesh", "map_extra")
            cmd.color("grey", "map_extra")
            cmd.center("pdb and site", 0, 1)
            cmd.zoom("pdb and site", "5.0")
            self.map_all()
            cmd.bg_color("white")
            cmd.show("nb_spheres", "pdb and site")
            cmd.show("cartoon", "all")
            cmd.ray("1024", "1024")
            cmd.png(f"%s/%s_%s.png" % (f"{self.path}/extrapolated_map",self.name, "cartoon"))
            cmd.delete("ccp4_map")
            cmd.delete("map_extra")





    def generate_image_default(self):
        cmd.select("site", self.resilist)
        if cmd.count_atoms("site") == 0:
            print("Selection 'site' has no atoms with valid coordinates.")
            return
        cmd.show("stick", "pdb and site")
        # cmd.center("pdb and site")
        cmd.hide("all")
        cmd.set("valence", "0")
        cmd.remove("hydrogen")

        cmd.turn("x", self.rot_x)
        cmd.turn("y", self.rot_y)
        cmd.turn("z", self.rot_z)
        cmd.set("ray_trace_fog", "0")
        cmd.set("depth_cue", "0")
        cmd.set("ray_shadows", "off")
        cmd.set("cartoon_transparency", "0.8", "all")





    def generate_image_diff_xtra(self):
        os.makedirs(f"{self.path}/extrapolated_diff_map",exist_ok=True)

        def generate_image_p(x):
            cmd.load(x, "ccp4_map")
            self.name = f"%s.png"%x.split(".")[0]
            cmd.do(f"isomesh map_extra_p,ccp4_map,%s, pdb and {self.resilist},carve=" % (self.map_sigma, self.carve))
            cmd.show("mesh", "map_extra_p")
            cmd.color("cyan", "map_extra_p")
            cmd.do(f"isomesh map_extra_n,ccp4_map,%s, pdb and {self.resilist},carve=" % (self.map_sigma, self.carve))
            cmd.show("mesh", "map_extra_n")
            cmd.color("orange", "map_extra_n")
            cmd.center("pdb and site", 0, 1)
            cmd.zoom("pdb and site", "5.0")
            self.map_all()
            cmd.bg_color("white")
            cmd.show("nb_spheres", "pdb and site")
            cmd.show("cartoon", "all")
            cmd.ray("1024", "1024")
            cmd.png(f"%s/%s_%s_diff.png" % (f"{self.path}/extrapolated_diff_map",self.name, "cartoon"))
            cmd.delete("ccp4_map")
            cmd.delete("map_extra_p")
            cmd.delete("map_extra_n")

        make_parallel(add_extra)(tqdm.tqdm(self.sorted_ccp4_map))






