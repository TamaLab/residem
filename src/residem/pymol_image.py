import pymol
from pymol import cmd
import numpy as np


class PymolImage:
    def __init__(self, pdb, ccp4_map, resilist="all", sigma=3.0, carve=1.6, positve=True, negative=True, all_map=True,
                 name="output", rot_x=0, rot_y=0, rot_z=0):
        self.pdb = cmd.load(pdb, "pdb")
        if "ccp4" in ccp4_map:
            self.map_to = True
            self.map = cmd.load(ccp4_map, "ccp4_map")
        else:
            self.map_to = False
        self.resilist = str(resilist)
        self.sigma = sigma
        self.carve = carve
        self.positive = positve
        self.negative = negative
        self.all_map = all_map
        self.name = name
        self.rot_x = rot_x
        self.rot_y = rot_y
        self.rot_z = rot_z
        print(self.resilist)

    def generate_image(self):
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
        if self.all_map and self.map_to:
            self.map_all()

        elif not self.all_map and self.positive and self.map_to:
            self.map_p()

        elif not self.all_map and self.negative and self.map_to:
            self.map_n()

        cmd.set("ray_trace_fog", "0")
        cmd.set("depth_cue", "0")
        cmd.set("ray_shadows", "off")
        cmd.set("cartoon_transparency", "0.8", "all")
        for x in ["stick,pdb and site", "cartoon,all"]:
            cmd.center("pdb and site", 0, 1)
            cmd.zoom("pdb and site", "5.0")
            # cmd.clip("slab",10,"pdb and site")
            cmd.bg_color("white")
            cmd.show("nb_spheres", "pdb and site")
            cmd.show(x.split(",")[0], x.split(",")[1])
            cmd.ray("1024", "1024")
            cmd.png(f"%s_%s.png" % (self.name, x.split(",")[0]))

    def map_all(self):
        cmd.do(f"isomesh map_p,ccp4_map,%s,pdb and {self.resilist},carve=%s" % (self.sigma, self.carve))
        cmd.do(f"isomesh map_n,ccp4_map,-%s,pdb and {self.resilist},carve=%s" % (self.sigma, self.carve))
        cmd.show("mesh", "map_n")
        cmd.show("mesh", "map_p")
        cmd.color("green", "map_p")
        cmd.color("red", "map_n")

    def map_p(self):
        cmd.do(f"isomesh map_p,ccp4_map,%s, pdb and {self.resilist},carve=%s" % (self.sigma, self.carve))
        cmd.show("mesh", "map_p")
        cmd.color("green", "map_p")

    def map_n(self):
        cmd.do(f"isomesh map_n,ccp4_map,-%s,pdb and {self.resilist},carve=%s" % (self.sigma, self.carve))
        cmd.show("mesh", "map_n")
        cmd.color("red", "map_n")




def run(pdb, ccp4_map, resilist, sigma, carve, positve, negative, all_map, name, rot_x, rot_y, rot_z):
    sc = PymolImage(pdb=pdb, ccp4_map=ccp4_map, resilist=resilist, sigma=sigma, carve=carve, positve=positve,
                    negative=negative, all_map=all_map, name=name,rot_x=rot_x, rot_y=rot_y,rot_z= rot_z)
    sc.generate_image()


def run_all():
    import argparse
    from datetime import datetime

    start = datetime.now()

    parser = argparse.ArgumentParser(prog='residem_pymol_image', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f"""It creates images of the residue with map cutoff using pymol \n
    Few example two run the code: \n\n
    1. residem_pymol_image -r 5b6v.pdb -m FO-FO.ccp4 -resi 'resi 182+216+300' -all True -o output_all -s 3.0 -c 1.6 \n\n
    2. residem_pymol_image -r 5b6v.pdb -m FO-FO.ccp4 -resi 'resi 182+216+300' -p True -o output_p \n\n
    3. residem_pymol_image -r 5b6v.pdb -m FO-FO.ccp4 -resi 'resi 182+216+300' -n True -o output_n \n\n
    4. residem_pymol_image -r 5b6v.pdb -m FO-FO.ccp4 -resi 'resi 182+216+300' -n True -o output_n  -x 90 -y 90 -z 90\n\n

    The resi syntax has to be followed as similar to the pymol selection syntax for reside selection. 
    The following link can be further noted for selection algebra.
    "https://pymolwiki.org/index.php/Selection_Algebra"
""")
    parser.add_argument("-r", "--pdb", type=str, help="pdb_file")
    parser.add_argument("-m", "--map", type=str, help="Fo-Fo_map_file in CCP4 format")
    parser.add_argument("-resi", "--resi", type=str, help="list of residue for which map and image has to be generated")
    parser.add_argument("-all", "--all", type=bool, help="If positive and negative map has to be printed")
    parser.add_argument("-p", "--positive", type=bool, help="If positive map alone has to be printed")
    parser.add_argument("-n", "--negative", type=bool, help="If negative map alone has to be printed")
    parser.add_argument("-s", "--sigma", type=float, help="sigma cutoff of the map")
    parser.add_argument("-c", "--carve", type=float, help="carve radius around the atom")
    parser.add_argument("-o", "--output", type=str, help="output name of the png")
    parser.add_argument("-x", "--rot_x", type=int, help="rotate x, angle")
    parser.add_argument("-y", "--rot_y", type=int, help="rotate y, angle")
    parser.add_argument("-z", "--rot_z", type=int, help="rotate y, angle")

    args = parser.parse_args()
    resi = args.resi if args.resi is not None else "all"
    sigma = args.sigma if args.sigma is not None else 3.0
    carve = args.carve if args.carve is not None else 1.6
    output = args.output if args.output is not None else "output"
    all_map = args.all if args.all is not None else True
    positive = args.positive if args.positive is not None else False
    negative = args.negative if args.negative is not None else False
    rot_x = args.rot_x if args.rot_x in list(np.arange(361)) else 0
    rot_y = args.rot_y if args.rot_y in list(np.arange(361)) else 0
    rot_z = args.rot_z if args.rot_z in list(np.arange(361)) else 0
    if positive:
        all_map = False
    elif negative:
        all_map = False

    run(args.pdb, args.map, resi, sigma, carve, positive, negative, all_map, output,rot_x,rot_y,rot_z)
    print(datetime.now() - start)


if __name__ == "__main__":
    run_all()
