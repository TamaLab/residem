from libtbx import easy_run
import os
import re
import subprocess



import residem

path_dir = os.path.dirname(os.path.abspath(residem.__file__))
path = path_dir +"/"+ "coot_view_peak.py"


class ShellRun:

    @staticmethod
    def easy_run_in(command, verbose=False):
        """Input the bash commandline argument for run"""
        if verbose:
            easy_run.fully_buffered(command).raise_if_errors().show_stdout()
        else:
            easy_run.fully_buffered(command)


class Coot(ShellRun):
    def __init__(self,pdb,map):
        self.pdb = pdb
        self.map = map
        self.write_image_path()


    def write_image_path(self):
        with open(path, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if re.search("pathdirc = os.getcwd()", line):
                lines[i] = f"pathdirc = \"%s\""%str(path_dir)
        
        with open(path, "w") as f:
            f.writelines(lines)



    def coot(self):

        template_f = f"""coot --pdb {self.pdb} --map {self.map} --script {path}"""
        # self.easy_run_in(template_f, False)
        try:
            output = subprocess.check_output(template_f, shell=True)
            print(output.decode())  # Print the command output
        except subprocess.CalledProcessError as e:
            print("Command execution failed:", e)


def run(pdb, map):
    sc = Coot(pdb=pdb, map=map)
    sc.coot()




def run_all():
    import argparse

    parser = argparse.ArgumentParser(prog='residem_coot', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f""" This is script opens coot and assumes coot to be installed and coot path is set.
    Example: \n\n
    $ residem_coot -r 5b6v -m F_obs_minus_F_obs.ccp4
    
    """)
    parser.add_argument("-r","--pdb", help="pdb_file")
    parser.add_argument("-m","--map", help="Fo-Fo_map_file in CCP4 format")
    args = parser.parse_args()
    run(args.pdb, args.map)

if __name__ == "__main__":
    run_all()
