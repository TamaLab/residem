from pymol import cmd


class SplitConfOcc:
    """This class creates alternative conformer based on the residue given and adds occupancy to the file"""

    def __init__(self, pdb, residue_number, chain=None, occupancy_b=0, output="merged"):
        self.pdb = pdb
        self.name = output

        self.residue_number = residue_number
        self.occupancy_b = float(occupancy_b) if occupancy_b < 1 else round(occupancy_b / 100, 2)
        self.residue_to_select = "+".join(self.residue_to_list())
        self.chain = " " if chain is None else chain
        self.save_pdb()

    def residue_to_list(self):
        residue_list = []
        for res in self.residue_number.split(","):
            if "_" in res:
                residue_list.append(str(res.replace("-", "-")))
            else:
                residue_list.append(str(res))

        residue_list.remove("") if "" in residue_list else ""
        return residue_list

    def save_pdb(self):
        cmd.load(self.pdb, "pdb")
        cmd.select(f"alt_residue", f" pdb and chain {self.chain} and resi {self.residue_to_select}")
        cmd.create("alt_b", "alt_residue")
        cmd.alter("alt_residue", "alt='A'")
        cmd.alter("alt_b", "alt='B'")
        cmd.alter("alt_residue", f"occupancy={1 - self.occupancy_b}")
        cmd.alter("alt_b", "occupancy=0.5")
        cmd.create("merged", "alt_b or pdb")
        cmd.save(f"{self.name}.pdb", "merged")


def main():
    import argparse
    parser = argparse.ArgumentParser(prog='residem_pdb', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f""" It creates pdb with alternative conformer for the selected 
                                     residue for refinement \n
        
        example:\n
       $ residem -r 5b6v.pdb -n 182,216-218,300 -o 30/0.30 -f merged\n
       $ residem -r 5b6v.pdb -n 182,216-218,300 -c A -o 30/0.30 -f merged\n
       $ residem -r 5b6v.pdb -t residue_list.txt -o 30/0.30 -f merged\n
       $ residem -r 5b6v.pdb -t residue_list.txt -c A -o 30/0.30 -f merged\n

""")
    parser.add_argument("-r", "--pdb", action="store", type=str,
                        help="pdb file")

    parser.add_argument("-n", "--residue", action='store', type=str,
                        help="list of residue list")

    parser.add_argument("-o", "--occu", action="store", type=float,
                        help="occupancy of the B conformer to be set")

    parser.add_argument("-c", "--chain", action="store", type=float,
                        help="chain name to be specified if any")

    parser.add_argument("-f", "--output", action="store", type=str,
                        help="output pdb name")
    parser.add_argument("-t", "--residue_list", action="store", type=str,
                        help="residue list file from coot")

    args = parser.parse_args()

    chain = args.chain if args.chain is not None else " "
    occupancy = args.occu if args.occu is not None else 0
    output = args.output if args.output is not None else "output"
    residue = args.residue if args.residue is not None and args.residue_list is None else open(args.residue_list).readline() if args.residue_list is not None else exit()
    SplitConfOcc(args.pdb, residue, chain, occupancy,output)


if __name__ == '__main__':
    main()
