import os
import shutil

import iotbx.phil
import iotbx.pdb
from iotbx.reflection_file_reader import any_reflection_file
import warnings
import sys
import math
from iotbx.map_model_manager import map_model_manager
import mmtbx.f_model
import logging
from residem.logger_set import setup_logger
import sympy as sp
from datetime import datetime


fo_minus_fo_master_params_st = """\
reference = None
    .type = path
    .help = reference state pdb file
reference_mtz = None
  .type = path
  .help = Reference MTZ file
  .short_caption = Reference Reflections
reference_label = "FP,SIGFP"
  .type = str
  .short_caption = Reference Label
  .input_size = 160
triggered_mtz = None
  .type = path
  .help = Triggered MTZ file
  .short_caption = Reflections of triggered state 
  .style = bold file_type:hkl process_hkl child:fobs:f_obs_2_label force_data
triggered_label = "FP,SIGFP"
  .type = str
  .short_caption = Triggered  Label
  .input_size = 160
fcalclabel = None
  .type = str
  .help = label of the refined f_calc file example FC,PHIC
high_resolution = None
  .type = float
  .help = High resolution data cutoff
low_resolution = None
  .type = float
  .help = Low resolution data cutoff
sigma = 3.0
  .type = float
  .help = sigma cutoff

weight = default
  .type = str
  .help = Weight used to calculate the difference map, the available option are 'ursby','ren','hekstra' and 'default'.
  
alpha = 0.05
  .type = float
  .help = Value of alpha in hekstra Weight used to calculate the difference map, It ccan be set to 0.05 or 1.0 .
  
  
scale = aniso
  .type = str
  .help = anisotropic or isotropic scaling of the native and derivative data, available option, aniso or iso or no
                 
atom_profile = mean
  .type = str
  .help = plots mean or max profile of the atom.
  
scaling_method = cctbx
  .type = str
  .help = scaling method is done either using cctbx or ccp4
  
ccp4_scaling_method = linear
  .type = str
  .help = CCP4 scaling method  if the scaling is done using ccp4, the avalible option is linear 

    
  """


############Functions##############
def deg_convert(theta):
    """ convert negative thetha to positive radians and converts as theta to radians
    """
    theta = math.fmod(math.fmod(theta, 360) + 360, 360) if theta < 0 else math.radians(theta)
    return theta


###################


class DataInput:
    def __init__(self, reference_pdb, reference_mtz, Triggered_mtz, reference_label=None, triggered_label=None,fcalclabel=None,
                 high_resolution=None,
                 low_resolution=None,logger=None):
        self.reference_mtz = reference_mtz
        self.Triggered_mtz = Triggered_mtz
        self.pdb = reference_pdb
        self.logger = setup_logger(logger)
        self.Dark_mtz = any_reflection_file(reference_mtz)
        self.Light_mtz = any_reflection_file(Triggered_mtz)

        self.dark_dict = self.map_model_manager()[0]
        self.dark_col = list(self.dark_dict.keys())
        self.light_dict = self.map_model_manager()[1]
        self.light_col = list(self.light_dict.keys())
        self.reference_label = reference_label if reference_label is not None else \
        self.map_obs_label(reference_label, triggered_label,fcalclabel)[0]
        self.triggered_label = triggered_label if triggered_label is not None else \
        self.map_obs_label(reference_label, triggered_label,fcalclabel)[1]
        self.label_phase = self.map_obs_label(reference_label, triggered_label,fcalclabel)[2].split(",")


        self.high_resolution = high_resolution if high_resolution is not None else self.map_resolution()[1]
        self.low_resolution = low_resolution if low_resolution is not None else self.map_resolution()[0]
        self.logger.info(self.map_resolution()[2])

        self.pdb_in = iotbx.pdb.input(file_name=self.pdb)
        self.map_iso = self.map_isomorphous()
        self.space_group_number = self.pdb_in.crystal_symmetry().space_group_number()
        self.logger.info(f"""The crystal space group number from the pdb is {self.pdb_in.crystal_symmetry_from_cryst1().space_group_number()},\nthe space group number of the reference and triggered file are {self.Dark_mtz.file_content().space_group_number()},{self.Light_mtz.file_content().space_group_number()} respectively.""")
        self.logger.info(f"""The column label used for reference data set is {self.reference_label}\nthe column label used for triggered data set is {self.triggered_label}\nthe column label used for F-calculated from reference data set is {self.label_phase}""")



    def __repr__(self):
        return f'DataInput(refernce pdb={self.pdb}, reference mtz={self.reference_mtz} , Triggered mtz - {self.Triggered_mtz},' \
               f'label {self.reference_label},{self.triggered_label}, high resolution {self.high_resolution}' \
               f', low resolution {self.low_resolution}, phase_label={self.label_phase} )'

    def map_isomorphous(self):
        try:
            space_group_number_pdb = self.pdb_in.crystal_symmetry_from_cryst1().space_group_number()
            space_group_number_dark_hkl = self.Dark_mtz.file_content().space_group_number()
            space_group_number_light_hkl = self.Light_mtz.file_content().space_group_number()

            if space_group_number_pdb == space_group_number_dark_hkl == space_group_number_light_hkl:
                return True
            else:
                print("The given model and map are not isomorphous")

        except AttributeError:
            print("The space group is not defined, The program will quit")

        return False

    def map_model_manager(self):
        dark_hkl = self.Dark_mtz.as_miller_arrays()
        light_hkl = self.Light_mtz.as_miller_arrays()
        dark_dict = {ma.info().label_string(): ma for ma in dark_hkl}
        light_dict = {ma.info().label_string(): ma for ma in light_hkl}

        return dark_dict, light_dict

    def map_resolution(self):
        if self.reference_label is not None and self.triggered_label is not None:
            dark_dmax = round(self.dark_dict[self.reference_label].d_max_min()[0], 1)
            dark_dmin = round(self.dark_dict[self.reference_label].d_max_min()[1], 1)
            light_dmax = round(self.light_dict[self.triggered_label].d_max_min()[0], 1)
            light_dmin = round(self.light_dict[self.triggered_label].d_max_min()[1], 1)
            message = f"""The Resolution cutoff of the Dark and light are not matching\nThe high resolution are ({dark_dmax},{light_dmax}) and low resolution are \n({dark_dmin},{light_dmin}), taking the lowest resolution cutoff"""

            if dark_dmin != light_dmin or dark_dmax != light_dmax:

                # self.logger.info(message)
                warnings.warn(message)
                # self.logger.info(message)
                return min(light_dmax, dark_dmax), max(dark_dmin, light_dmin),message
            else:
                return dark_dmax, dark_dmin,message

    def map_obs_label(self, A, B,C):
        label_dark = [X for X in list(self.dark_col) if
                      X == 'FP,SIGFP' or X == 'F,SIGF' or X == "F-obs-filtered,SIGF-obs-filtered" or X=="I-obs,SIGI-obs" or X=="I,SIGI" or X == A][0]
        label_light = [X for X in list(self.light_col) if
                       X == 'FP,SIGFP' or X == 'F,SIGF' or X == "F-obs-filtered,SIGF-obs-filtered" or X=="I-obs,SIGI-obs" or X=="I,SIGI" or X == B][0]

        label_phase = [X for X in list(self.dark_col) if X == 'FC,PHIC' or X == 'FC_ALL,PHIC_ALL' or X == 'FC,PHIFC' or X == "F-model,PHIF-model" or X == C][0]


        return label_dark, label_light, label_phase



def write_phil(output_path,refernce_pdb, refernce_mtz, triggered_mtz, reference_label, triggred_label, high_resolution, low_resolution,
        sigma, weight, scale,atom_profile,alpha,scaling_method,ccp4_scaling_method,fcalclabel):

    params = iotbx.phil.parse(fo_minus_fo_master_params_st, process_includes=True)
    working_phil = params.extract()
    working_phil.reference = refernce_pdb
    working_phil.reference_mtz = refernce_mtz
    working_phil.reference_label = reference_label
    working_phil.triggered_mtz = triggered_mtz
    working_phil.triggered_label = triggred_label
    working_phil.high_resolution = high_resolution
    working_phil.low_resolution = low_resolution
    working_phil.sigma = sigma
    working_phil.weight = weight
    working_phil.alpha = alpha
    working_phil.fcalclabel = fcalclabel
    working_phil.scale = scale
    working_phil.atom_profile = atom_profile
    working_phil.scaling_method = scaling_method
    working_phil.ccp4_scaling_method = ccp4_scaling_method
    modified_phil = params.format(python_object=working_phil)
    modified_phil.show(out=open(f"{output_path}/input.phil", "w"), attributes_level=2)





def run(refernce_pdb, refernce_mtz, triggered_mtz, reference_label, triggred_label, high_resolution, low_resolution,
        sigma, weight, scale,atom_profile,alpha,scaling_method,ccp4_scaling_method,fcalclabel):
    start_calculation = datetime.now()
    from residem.Extrapolate import ExtrapoledMap

    from residem.MapradiusSelection import MapPeakSelection_coff
    import residem.atompeak_map
    from residem.ccp4_scaleit import Run_scaleit
    refernce_pdb_path = os.path.realpath(refernce_pdb)
    refernce_mtz_path = os.path.realpath(refernce_mtz)
    triggered_mtz_path = os.path.realpath(triggered_mtz)
    I = 0
    name = "Data_folder"
    exists = True
    while exists:
        if not os.path.exists(name + "_%s" % I):
            os.mkdir(name + "_%s" % I)
            exists = False
        else:
            I = I + 1

    name1 = name + "_%s" % I
    name_path = os.getcwd() + f"/%s" % name1

    write_phil(name_path, refernce_pdb_path, refernce_mtz_path, triggered_mtz_path, reference_label, triggred_label, high_resolution,
               low_resolution,sigma, weight, scale,  atom_profile,
               alpha, scaling_method, ccp4_scaling_method,fcalclabel)

    log_file = "Residem.log"
    log_path = os.path.join(name_path, log_file)
    logger = setup_logger(log_path)
    logger.info(f"""{"#"*20} 
Running Residem calculation for the identification of the residue in difference density map.
Reference data set consists of pdb file: {refernce_pdb} and {refernce_mtz}.
Triggered data set with mtz {triggered_mtz} is taken for calculation. 
The input file for the current calculation is written as {f"{name_path}/input.phil"} which can be used to run to recalculation.""")



    if scaling_method == "cctbx":

        Map_file = ExtrapoledMap(refernce_pdb_path, refernce_mtz_path, triggered_mtz_path, reference_label=reference_label,
                                 triggered_label=triggred_label, high_resolution=high_resolution,fcalclabel=fcalclabel,
                                 low_resolution=low_resolution, sigma=sigma, grid_step=0.5,
                                 weight=weight, write=True, path=name_path, scale=scale,
                                 alpha=alpha,logger=log_path)
    elif scaling_method == "ccp4":
        run_scaleit = Run_scaleit(refernce_pdb_path, refernce_mtz_path, triggered_mtz_path,reference_label=reference_label,triggered_label=triggred_label,fcalclabel=fcalclabel, high_resolution=high_resolution,
                                 low_resolution=low_resolution,scaling=ccp4_scaling_method,path=name_path,logger=log_path)
        result_mtz_path = os.path.realpath(f"%s/Light-dark-scaled.mtz"% name_path)

        Map_file = ExtrapoledMap(refernce_pdb_path, result_mtz_path, result_mtz_path, reference_label="F_DARK,SIGF_DARK",
                                 triggered_label="F_LIGHT,SIGF_LIGHT", high_resolution=high_resolution,
                                 low_resolution=low_resolution, sigma=sigma, grid_step=0.5,
                                 weight=weight, write=True, path=name_path, scale="no",
                                  alpha=alpha, logger=log_path)



    dict_map = {"ren": "F_obs_minus_F_obs_ren_weight.ccp4", "ursby": "F_obs_minus_F_obs_ursby_weight.ccp4",
                "default": "F_obs_minus_F_obs.ccp4"}



    if weight == "ren":
        map_model = map_model_manager(map_manager=Map_file.FOFODF_K.as_map_manager(), model=Map_file.model)


        to_identify = MapPeakSelection_coff(map_model, Map_file.high_resolution, sigma, name=weight, path=name_path,write=True)

        residem.atompeakheight.AtomPeaks(refernce_pdb_path, map_model,
                                        Map_file.high_resolution, sigma,
                                        path=name_path,
                                        residue=None, weight="ren",atom_profile=atom_profile)
    elif weight == "ursby":
        map_model = map_model_manager(map_manager=Map_file.FOFODF_Q.as_map_manager(), model=Map_file.model)



        to_identify = MapPeakSelection_coff(map_model, Map_file.high_resolution, sigma, name=weight, path=name_path,write=True)

        residem.atompeakheight.AtomPeaks(refernce_pdb_path, map_model,
                                        Map_file.high_resolution, sigma,
                                        path=name_path,
                                        residue=None, weight="ursby",atom_profile=atom_profile)
    elif weight == "hekstra":
        map_model = map_model_manager(map_manager=Map_file.FOFODF_H.as_map_manager(), model=Map_file.model)


        to_identify = MapPeakSelection_coff(map_model, Map_file.high_resolution, sigma, name=weight, path=name_path,write=True)

        residem.atompeakheight.AtomPeaks(refernce_pdb_path, map_model,
                                        Map_file.high_resolution, sigma,
                                        path=name_path,
                                        residue=None, weight="hekstra",atom_profile=atom_profile)


    else:
        map_model = map_model_manager(map_manager=Map_file.FOFODF.as_map_manager(), model=Map_file.model)

        to_identify = MapPeakSelection_coff(map_model, Map_file.high_resolution, sigma, name=weight, path=name_path,write=True)

        residem.atompeakheight.AtomPeaks(refernce_pdb_path, map_model,
                                        Map_file.high_resolution, sigma,
                                        path=name_path,
                                        residue=None, weight="default",atom_profile=atom_profile)

    logger.info(f"""{"#" * 20}""" )
    logger.info(f"""Clustering Done.
All the result of the clustering are inside the folder {f"{name_path}/map_dump_default"} which can be used for further calculations.
The atom peak profile is also available in the folder with chain identifier 
example:{f"{name_path}/map_dump_default/chain_A_U/Atom_peak_height_chain_A_U.pdf"} or
{f"{name_path}/map_dump_default/chain_A_U/Residual_peak_height_max_chain_A_U.pdf"}""")
    logger.info(f"""{"#" * 20}""" )
    logger.info(f"""Time taken for over all calculation is {datetime.now() - start_calculation}""")




def run_all():
    import argparse
    from datetime import datetime



    start = datetime.now()
    parser = argparse.ArgumentParser(prog='residem', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f""" This is the main script which calculates the difference map and occupancy of reference and triggred data. 
        Few example to run the code: \n\n
        1. residem -r 5b6v.pdb -m 5b6v.mtz -t 5b6w.mtz  \n\n
        2. residem -r 5b6v.pdb -m 5b6v.mtz -t 5b6w.mtz  --scale aniso/iso -w ursby/ren/hekstra/default -sm cctbx/ccp4 
        3. residem -r 2vja.pdb -m 2vja.mtz -t 2vjb.mtz -w ursby 
    """)

    parser.add_argument("-p", "--phil", type=str, help="process parameter file")
    parser.add_argument("-d", "--default", action='store_true', help="write default parameter file")
    parser.add_argument("-r", "--ref", type=str, help="reference pdb file")
    parser.add_argument("-m", "--mtz", type=str, help="reference hkl file in mtz format")
    parser.add_argument("-t", "--tmtz", type=str, help="triggered hkl file in mtz format")
    parser.add_argument("-hi", "--high", type=float, help="high resolution")
    parser.add_argument("-l", "--low", type=float, help="low resolution")
    parser.add_argument("-rl", "--refl", action='store', help="label of the reference file ")
    parser.add_argument("-tl", "--tril", action='store', help="label of the triggred file")
    parser.add_argument("-fl", "--fcalclabel", action='store', help="label of the refined f_calc file example FC,PHIC")

    parser.add_argument("-s", "--sigma", type=float, help="sigma cutoff for map processing")
    parser.add_argument("-w", "--weight", type=str, help="weight for choosing the map calculation")
    parser.add_argument("-ws","--alpha",type=float,help = "value of alpha in Hekstra weight, default 0.05 ")
    parser.add_argument("--scale", type=str,
                        help="scaling if it is 'iso'(isotropic scaling) or 'aniso'(anisotropic scalig) or (no) scaling")
    parser.add_argument("-v", "--atom_profile", type=str,
                        help="If difference density profile if max or mean")
    parser.add_argument("-sm", "--scaling_method", type=str,
                        help="scaling method is done either using cctbx or ccp4")

    parser.add_argument("-cs", "--ccp4_scaling_method", type=str,
                        help="scaling method is done using linear")

    args = parser.parse_args()




    if args.default:
        global fo_minus_fo_master_params_st
        with open("Residem.phil", 'w') as f:
            f.write(fo_minus_fo_master_params_st)
            f.close()
        print("Please provide a modified parameter file,similar to Residem.phil")
        parser.print_help()

    elif args.phil:
        fo_minus_fo_master_params_st = " ".join(open(args.phil, "r").readlines())
        params = iotbx.phil.parse(fo_minus_fo_master_params_st, process_includes=True).extract()
        reference_pdb = params.reference if params.reference is not None else args.ref
        reference_mtz = params.reference_mtz if params.reference_mtz is not None else args.mtz
        Triggered_mtz = params.triggered_mtz if params.triggered_mtz is not None else args.tmtz
        reference_label = params.reference_label if params.reference_label is not None else None
        triggered_label = params.triggered_label if params.triggered_label is not None else None
        fcalclabel = params.fcalclabel if params.fcalclabel is not None else None
        high_resolution = params.high_resolution if params.high_resolution is not None else None
        low_resolution = params.low_resolution  if params.low_resolution is not None else None
        phil_sigma = params.sigma if params.sigma is not None else args.sigma if args.sigma is not None else 3.0
        phil_grid = params.grid_step
        phil_weight = params.weight  if params.weight is not None else "default"
        phil_alpha = params.alpha if params.alpha is not None else 0.05
        scale = params.scale if params.scale in ["iso", "aniso","no"] else "aniso"
        atom_profile = params.atom_profile
        scaling_method = params.scaling_method

        if scaling_method == "ccp4":
            scale = "no"
        ccp4_scaling_method = params.ccp4_scaling_method

        run(reference_pdb, reference_mtz, Triggered_mtz, reference_label, triggered_label, high_resolution, low_resolution, phil_sigma, phil_weight,scale, atom_profile, phil_alpha,scaling_method,ccp4_scaling_method,fcalclabel)


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif not args.phil:
        refl = args.refl if args.refl is not None else None
        tril = args.tril if args.tril is not None else None
        fcalclabel = args.fcalclabel if args.fcalclabel is not None else None
        high = args.high if args.high is not None else None
        low = args.low if args.low is not None else None
        sigma = args.sigma if args.sigma is not None else 3.0
        scale = args.scale if args.scale in ["iso", "aniso","no"] else "aniso"
        atom_profile = args.atom_profile if args.atom_profile is not None else "mean"
        weight = "default" if args.weight is None else args.weight.lower() if args.weight.lower() in ["ren",
                                                                                                      "ursby","hekstra"] else "default"
        alpha = args.alpha if args.alpha is not None else 0.05
        scaling_method = args.scaling_method if args.scaling_method in ["cctbx","ccp4"] else "cctbx"
        ccp4_scaling_method = args.ccp4_scaling_method if args.scaling_method == "ccp4" and args.ccp4_scaling_method in ["linear"] else "linear"
        if scaling_method == "ccp4":
            scale = "no"
        if args.ref and args.mtz and args.tmtz:
            run(refernce_pdb=args.ref, refernce_mtz=args.mtz, triggered_mtz=args.tmtz, reference_label=refl, triggred_label=tril,
                high_resolution=high, low_resolution=low, sigma=sigma, weight=weight, scale=scale,atom_profile=atom_profile,
                alpha=alpha,scaling_method=scaling_method,ccp4_scaling_method=ccp4_scaling_method,fcalclabel=fcalclabel)
        else:
            print("Please provide the reference pdb, reference hkl file in mtz format and triggered file in mtz format")

    print(datetime.now() - start)


if __name__ == '__main__':
    run_all()
