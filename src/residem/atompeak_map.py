import residem.atompeakheight
from residem.MapradiusSelection import MapPeakSelection
import os


def run_atom_map(refernce_pdb, map, high_resolution, sigma, max_dist, min_dist, name, path, weight, map_do=True):
    if map_do:
        MapPeakSelection(pdb_file=refernce_pdb, map_file=map, resolution=high_resolution, sigma=sigma,
                         max_dist=max_dist, min_dist=min_dist, name=name, path=path)
    residem.atompeakheight.AtomPeaks(pdb_file=refernce_pdb, map_file=map, resolution=high_resolution, sigma=sigma,
                                    name=name, path=path, weight=weight)


def run_all():
    import argparse

    parser = argparse.ArgumentParser(prog='residem_map_identify', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f"""
    It identifies the residues corresponding to the isomorphous difference map. 

    1. residem_map_identify -r 5b6v.pdb -m FO-FO.ccp4 -h 1.4 -s 3.0 -d 2.0 -o map_output\n\n
    2. residem_map_identify -r 5b6v.pdb -m FO-FO.ccp4 -h 1.4 \n\n
""")
    parser.add_argument("-r", "--pdb", type=str, help="pdb_file")
    parser.add_argument("-m", "--map", type=str, help="Fo-Fo_map_file in CCP4 format")
    parser.add_argument("-h", "--high", type=str, help="high resolution")
    parser.add_argument("-s", "--sigma", type=float, help="sigma cutoff")
    parser.add_argument("-d", "--max", type=float, help="maximum distance to identify atom")
    parser.add_argument("-o", "--output", type=str, help="output name of the png")

    args = parser.parse_args()
    sigma = args.sigma if args.sigma is not None else 3.0
    max = args.max if args.max is not None else 2.0
    output = args.output if args.output is not None else "output"

    name = "map_Data_folder"
    exists = True
    while exists:
        if not os.path.exists(name + "_%s" % I):
            os.mkdir(name + "_%s" % I)
            exists = False
        else:
            I = I + 1

    name1 = name + "_%s" % I
    name_path = os.getcwd() + f"/%s" % name1

    run_atom_map(args.pdb, args.map, args.high, sigma, max, 0, output, name_path, "default", map_do=True)


if __name__ == "__main__":
    run_all()
