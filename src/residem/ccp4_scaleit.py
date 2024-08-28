import subprocess
import os
from subprocess import PIPE, Popen, run
from residem.master import DataInput
import time
from residem.logger_set import setup_logger


class Run_scaleit(DataInput):
    def __init__(self, reference_pdb, reference_mtz, Triggered_mtz, reference_label=None, triggered_label=None,fcalclabel=None,
                 high_resolution=None,
                 low_resolution=None,scaling="linear",path=os.getcwd(),logger=None):
        super().__init__(reference_pdb, reference_mtz, Triggered_mtz, reference_label, triggered_label,fcalclabel,
                 high_resolution,
                 low_resolution)
        self.logger = logger
        self.high_resolution = high_resolution if high_resolution is not None else self.high_resolution
        self.low_resolution = low_resolution if low_resolution is not None else self.low_resolution

        # self.add_sigma()
        self.path = path
        self.run_cad_r = self.run_cad()
        # if scaling == "linear":
        self.run_scaleit_r = self.run_scaleit()
        # elif scaling == "multi":
        # self.run_scaleit_r = self.run_scaleit_advanced()



    def add_sigma(self):
        """Adding sigma column to the mtz file"""
        template = f"""mtz2various HKLIN {self.reference_mtz} HKLOUT {self.path}/model_phs.hkl<<EOF
LABIN FP={self.label_phase[0]} PHIC={self.label_phase[1]}
OUTPUT USER '(3I5,F12.3,'  1.00  ',F12.3)'
EOF"""
        symm = self.pdb_in.crystal_symmetry().unit_cell().parameters()

        template_1 =f"""f2mtz HKLIN model_phs.hkl HKLOUT {self.path}/FC_dark.mtz << f2m_phs

CELL {symm[0]}   {symm[1]}   {symm[2]} {symm[3]} {symm[4]} {symm[5]}  # angles default to 90
SYMM {self.space_group_number}
LABOUT H   K  L   FC_D SIG_FC_D PHI_D
CTYPE  H   H  H   F     Q        P
f2m_phs"""


        return self.run_result(template),self.run_result(template_1)

    def run_cad(self):
        template = f"""cad hklin1 {self.reference_mtz} hklin2 {self.Triggered_mtz}  \
hklout {self.path}/cad_light_dark.mtz <<EOF
LABIN FILE_NUMBER 1 E1={self.reference_label.split(",")[0]} E2={self.reference_label.split(",")[1]} E4={self.label_phase[1]} E3={self.label_phase[0]}
LABIN FILE_NUMBER 2 E1={self.triggered_label.split(",")[0]} E2={self.triggered_label.split(",")[1]}
LABOUT FILE_NUMBER 1 E1=F_DARK E2=SIGF_DARK E4={self.label_phase[1]} E3={self.label_phase[0]}
LABOUT FILE_NUMBER 2 E1=F_LIGHT E2=SIGF_LIGHT 
XNAME FILE_NUMBER 1 ALL=DARK
XNAME FILE_NUMBER 2 ALL=LIGHT
resolution file 1 {self.low_resolution} {self.high_resolution}
resolution file 2 {self.low_resolution} {self.high_resolution}
EOF"""
        self.logger.info(f"Running CCP4 CAD to combine the refernce and triggered data with the following command\n {template}")
        return self.run_result(template)

    def run_scaleit(self):

        template = f"""scaleit hklin {self.path}/cad_light_dark.mtz hklout {self.path}/Light-dark-scaled.mtz <<EOF | tee {self.path}/Light-dark_scaleit.log
LABIN FP=F_DARK SIGFP=SIGF_DARK FPH1=F_LIGHT SIGFPH1=SIGF_LIGHT
EOF"""
        self.logger.info(f"Running CCP4 scaleit to scale the refernce and triggered data with the following command\n {template}")
        return self.run_result(template)

    def run_scaleit_advanced(self):
        template = f"""scaleit hklin {self.path}/cad_light_dark.mtz hklout {self.path}/Light-dark-scaled_model.mtz <<EOF | tee {self.path}/Light-dark_scaleit_model.log
        TITLE FPHs scaled to FP
        reso {self.low_resolution} {self.high_resolution}  
        EXCLUDE FP SIG 4 FMAX 10000000
        REFINE ANISOTROPIC
        LABIN FP={self.label_phase[0]} SIGFP=SIGF_DARK FPH1=F_DARK SIGFPH1=SIGF_DARK FPH2=F_LIGHT SIGFPH2=SIGF_LIGHT
        CONV ABS 0.0001 TOLR  0.000000001 NCYC 150
        END
        EOF
        """
        self.run_result(template)

        template_1 = f"""scaleit hklin {self.path}/Light-dark-scaled_model.mtz hklout {self.path}/Light-dark-scaled.mtz <<EOF | tee {self.path}/Light-dark_scaleit.log
        TITLE FPHs scaled to FP
        reso {self.low_resolution} {self.high_resolution}  
        EXCLUDE FP SIG 4 FMAX 10000000
        REFINE ANISOTROPIC
        LABIN FP=F_DARK SIGFP=SIGF_DARK FPH1=F_LIGHT SIGFPH1=SIGF_LIGHT
        CONV ABS 0.0001 TOLR  0.000000001 NCYC 150
        END
        EOF
        """
        self.run_result(template_1)


    def run_result(self, template):
        try:
            result = subprocess.run(template, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, shell=True, check=True)
            return result.returncode
        except subprocess.CalledProcessError as e:
            print(f"Background error: {e.stderr}")
            return e.returncode


def run_all():
    import argparse
    parser = argparse.ArgumentParser(description='Run scaleit.')
    parser.add_argument('-r', '--reference_pdb', type=str,
                        required=True, help='Reference PDB file')
    parser.add_argument('-m', '--reference_mtz', type=str,
                        required=True, help='Reference MTZ file')
    parser.add_argument('-t', '--Triggered_mtz', type=str,
                        required=True, help='Triggered MTZ file')
    parser.add_argument("-rl", "--refl", action='store',
                        help="label of the reference file ")
    parser.add_argument("-tl", "--tril", action='store',
                        help="label of the triggred file")
    parser.add_argument("-hi", "--high", type=float, help="high resolution")
    parser.add_argument("-l", "--low", type=float, help="low resolution")

    args = parser.parse_args()
    refl = args.refl if args.refl is not None else None
    tril = args.tril if args.tril is not None else None
    high = args.high if args.high is not None else None
    low = args.low if args.low is not None else None


    run_scaleit = Run_scaleit(reference_pdb=args.reference_pdb,
                              reference_mtz=args.reference_mtz, Triggered_mtz=args.Triggered_mtz, reference_label=refl, triggered_label=tril,high_resolution=high,low_resolution=low)


if __name__ == "__main__":
    run_all()
