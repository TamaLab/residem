import os

import pymol

from pymol import cmd
import subprocess


class PymolImage_ALT:
    def __init__(self, pdb1, pdb2, ccp4_map, diff_map, resilist="all", sigma=3.0, carve=1.6, name="output"):
        self.pdb1 = cmd.load(pdb1, "pdb1")
        self.pdb2 = cmd.load(pdb2, "pdb2")
        self.ccp4_map = cmd.load(ccp4_map, "ccp4_map")
        self.diff_map = cmd.load(diff_map, "diff_map")
        self.resilist = str(resilist)
        self.sigma = sigma
        self.carve = carve
        self.name = name
        print(self.resilist)

    def generate_image(self):
        cmd.select("site", self.resilist)
        if cmd.count_atoms("site") == 0:
            print("Selection 'site' has no atoms with valid coordinates.")
            return
        cmd.show("stick", "pdb1 and site and alt A")
        cmd.color("green", "pdb1 and site and alt A and elem C")
        cmd.show("stick", "pdb2 and site and alt A")
        cmd.color("green", "pdb2 and site and alt A and elem C")
        cmd.show("stick", "pdb1 and site and alt B")
        cmd.color("yellow", "pdb1 and site and alt B and elem C")
        cmd.show("stick", "pdb2 and site and alt B")
        cmd.set_color("coral", [1.0, 0.7, 0.2])  # Set custom RGB color
        cmd.color("coral", "pdb2 and site and alt B and elem C")
        cmd.hide("all")
        cmd.set("valence", "0")
        cmd.remove("hydrogen")

        cmd.do(f"isomesh map, ccp4_map, 1.0, pdb1 and {self.resilist}, carve=%s" % self.carve)
        cmd.show("mesh", "map")
        cmd.color("grey50", "map")
        cmd.set("transparency", 0.5, "map")  # Set transparency to 0.5 (50% transparency) for "map_p"
        cmd.do(f"isomesh map_p, diff_map, 3.0, pdb1 and {self.resilist}, carve=%s,buffer=4" % 3)
        cmd.set("transparency", 0.8, "map_p")
        cmd.show("mesh", "map_p")
        cmd.color("green", "map_p")
        cmd.do(f"isomesh map_n, diff_map, -3.0, pdb1 and {self.resilist}, carve=%s,buffer=4" % 3)
        cmd.show("mesh", "map_n")
        cmd.color("red", "map_n")
        cmd.set("transparency", 0.8, "map_n")
        cmd.set("ray_trace_fog", "1")
        cmd.set("depth_cue", "0")
        cmd.set("ray_shadows", "on")
        cmd.set("cartoon_transparency", "0.8", "all")

        for x in ["stick,pdb1 and site or pdb2 and site ", "cartoon,all"]:
            cmd.orient("pdb1 and site")
            cmd.center("pdb1 and site")
            cmd.zoom("pdb1 and site",5,complete=1)
            cmd.show(x.split(",")[0], x.split(",")[1])
            cmd.bg_color("white")
            cmd.ray(1024, 1024)
            cmd.png(f"{self.name}_{x.split(',')[0]}.png")

        cmd.show("stick", "pdb1 and site")
        cmd.show("stick", "pdb2 and site")
        cmd.zoom("pdb1 and site", complete=1)
        cmd.hide("cartoon", "all")
        cmd.mset("1 x30")
        cmd.mview("store")
        cmd.viewport(800, 600)
        os.makedirs(f"movie_{self.name}", exist_ok=True)
        for i in range(1, 21):
            cmd.mview("clear")
            cmd.frame(i)
            cmd.turn("y", 18)
            cmd.png(f"movie_%s/%s_movie_%02d.png" % (self.name,self.name, i))

        for i in range(21, 41):
            cmd.mview("clear")
            cmd.frame(i)
            cmd.turn("x", 18)
            cmd.png(f"movie_%s/%s_movie_%02d.png" % (self.name,self.name, i))


# PymolImage_ALT("pdb1.pdb", "pdb2.pdb", "2fo-fc.map", "mfo-fc.map", resilist="alt A or alt B").generate_image()


def run(pdb1, pdb2, ccp4_map, diff_map, resilist, sigma, carve, name,rot_x,rot_y,rot_z):
    sc = PymolImage_ALT(pdb1=pdb1, pdb2=pdb2, ccp4_map=ccp4_map, diff_map=diff_map, resilist=resilist, sigma=sigma,
                        carve=carve, name=name,)
    sc.generate_image()
    command_1 = f'ffmpeg -i movie_{name}/{name}_movie_%02d.png -vf "setpts=10.0*PTS" -s 960x960 -c:v mjpeg -q:v 0 movie_{name}/{name}.avi'
    command_2 = f'ffmpeg -i movie_{name}/{name}.avi -c:v libx264 -vf "scale=1280:720" -crf 18 -b:v 4000k -pix_fmt yuv420p movie_{name}/{name}.mp4'

    # command_1 = f'ffmpeg -i movie_{name}/{name}_movie_%02d.png -vcodec mjpeg -vf "setpts=2.0*PTS" -s 960x960 movie_{name}/{name}.avi'
    # command_2 = f"ffmpeg -i movie_{name}/{name}.avi -b:v 4000k -vcodec libx264 -pix_fmt yuv420p movie_{name}/{name}.mp4"
    subprocess.run(command_1, shell=True)
    subprocess.run(command_2, shell=True)


def run_all():
    import argparse
    from datetime import datetime

    start = datetime.now()

    parser = argparse.ArgumentParser(prog='residem_pymol_image', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f"""It creates images of the residue with map cutoff using pymol \n
    Few example two run the code: \n\n
    1. residem_pymol_alt_viz -r1 5b6x.pdb  -r2 MD1.pdb -m1 2FO-Fc.ccp4 -m2 FO-FO.ccp4 -resi 'resi 182+216+300' -o output_all -s 3.0 -c 1.6 \n\n


    The resi syntax has to be followed as similar to the pymol selection syntax for reside selection. 
    The following link can be further noted for selection algebra.
    "https://pymolwiki.org/index.php/Selection_Algebra"
""")
    parser.add_argument("-r1", "--pdb1", type=str, help="pdb_file 1")
    parser.add_argument("-r2", "--pdb2", type=str, help="pdb_file 2")
    parser.add_argument("-m1", "--map1", type=str, help="2Fo-Fc_map_file in CCP4 format")
    parser.add_argument("-m2", "--map2", type=str, help="Fo-Fo_map_file in CCP4 format")
    parser.add_argument("-resi", "--resi", type=str, help="list of residue for which map and image has to be generated")
    parser.add_argument("-s", "--sigma", type=float, help="sigma cutoff of the map")
    parser.add_argument("-c", "--carve", type=float, help="carve radius around the atom")
    parser.add_argument("-o", "--output", type=str, help="output name of the png")
    parser.add_argument("-x", "--rot_x", type=int, help="rotate x, angle")
    parser.add_argument("-y", "--rot_y", type=int, help="rotate y, angle")
    parser.add_argument("-z", "--rot_z", type=int, help="rotate y, angle")
    args = parser.parse_args()
    rot_x = args.rot_x if args.rot_x is not None and args.rot_x in list(np.arange(361)) else 0
    rot_y = args.rot_y if args.rot_y is not None and args.rot_y in list(np.arange(361)) else 0
    rot_z = args.rot_z if args.rot_z is not None and args.rot_z in list(np.arange(361)) else 0
    resi = args.resi if args.resi is not None else "all"
    sigma = args.sigma if args.sigma is not None else 3.0
    carve = args.carve if args.carve is not None else 1.6
    output = args.output if args.output is not None else "output"
    run(args.pdb1, args.pdb2, args.map1, args.map2, resi, sigma, carve, output,rot_x,rot_y,rot_z)
    print(datetime.now() - start)

if __name__ == "__main__":
    run_all()
