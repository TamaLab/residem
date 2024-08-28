import argparse

residem_overview = (""" `residem` is a command line utility that computes isomorphous difference density (DED) maps. 
You will need to provide one reference PDB file along with an MTZ file of the reference and triggered system. 
The input can be intensity profiles or structure factor amplitudes, which the program will try to detect labels for. 
If labels for reflection files are not provided, they must be supplemented. The ResiDEM performs scaling of the 
reference and triggered structure factors, followed by the calculation of the difference structure factor and DED maps.
""")

residem_pdb_mtz_overview = (""" `residem_pdb_mtz` is a command line utility used to download PDB and MTZ files from 
https://www.rcsb.org/.
""")

residem_coot_overview = (""" `residem_coot` is a command line utility used to visualize the DED maps. It uses Coot with 
an additional script to view residues associated with DED maps. It calls Coot and provides an input script for residual 
tracking associated with DED maps.
""")

residem_extract_around_model_overview = (""" `residem_extract_around_model` is a command line utility used to extract 
maps around a specified residue in the map.
""")

residem_svd_csv_overview = (""" `residem_svd_csv` is a command line utility used to analyze the results of 
difference density output. The utility can generate one-dimensional residue-wise plots, two-dimensional heat map plots, 
Pearson correlation plots, and network plots for residues with difference density in DED maps.
""")

residem_pymol_image_overview = (""" `residem_pymol_image` is a command line utility used to generate images of residues 
in DED maps using PyMOL.
""")

residem_args = (
    (("--phil", "-p"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """The PHIL file is a Python Hierarchical Input Language file, containing all inputs required to run the 
        ResiDEM command line utility. Provide it with one argument."""
    }),
    (("-d", "--default"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": None,
        "help": """Prints a default PHIL file which can be modified and provided as input."""
    }),
    (("--ref", "-r"), {
        "nargs": 1,
        "metavar": None,
        "required": True,
        "type": str,
        "default": None,
        "help": """Reference PDB file used to generate phases and perform other calculations such as scaling, figure of merit, 
        and phase error."""
    }),
    (("--mtz", "-m"), {
        "nargs": 1,
        "metavar": None,
        "required": True,
        "type": str,
        "default": None,
        "help": """Refined reference MTZ file with experimental observed intensities or structure factors. It should also have 
        FC and PHIC columns."""
    }),
    (("--tmtz", "-t"), {
        "nargs": 1,
        "metavar": None,
        "required": True,
        "type": str,
        "default": None,
        "help": """Triggered MTZ file with experimental observed intensities or structure factors."""
    }),
    (("--high", "-hi"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": None,
        "help": """High resolution cutoff. By default, the highest common resolution between the reference and triggered 
        structure factors is used."""
    }),
    (("--low", "-l"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": None,
        "help": """Low resolution cutoff. By default, the highest common resolution between the reference and triggered 
        structure factors is used."""
    }),
    (("--refl", "-rl"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": None,
        "help": """Label of the structure factor or intensity of the reference MTZ to be used for processing. By default, it 
        identifies from labels such as 'FP,SIGFP', 'F,SIGF', 'F-obs-filtered,SIGF-obs-filtered', 'I-obs,SIGI-obs', or 
        'I,SIGI'. If no label is identified, it must be provided."""
    }),
    (("--tril", "-tl"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": None,
        "help": """Label of the structure factor or intensity of the triggered MTZ to be used for processing. By default, it 
        identifies from labels such as 'FP,SIGFP', 'F,SIGF', 'F-obs-filtered,SIGF-obs-filtered', 'I-obs,SIGI-obs', or 
        'I,SIGI'. If no label is identified, it must be provided."""
    }),
    (("--fcalclabel", "-fl"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": None,
        "help": """Label of the calculated structure factor and phase of the reference MTZ to be used for processing. By 
        default, it identifies from labels such as 'FC,PHIC', 'FC_ALL,PHIC_ALL', 'FC,PHIFC', or 'F-model,PHIF-model'. 
        If no label is identified, it must be provided."""
    }),
    (("--sigma", "-s"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": float,
        "default": 3.0,
        "help": """Sigma cutoff for data processing. By default, it is 3.0 sigma."""
    }),
    (("--weight", "-w"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": None,
        "help": """Weighting scheme for the treatment of noise in DED maps."""
    }),
    (("--alpha", "-ws"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": float,
        "default": 0.05,
        "help": """Alpha for modified Ren weight. Use alpha=0.05 for Hekstra or alpha=1 for Schmidt weighting scheme."""
    }),
    (("--scale",), {
        "help": (
            """Scaling of structure factors between the triggered and reference systems (derivative to native scaling) 
        before generating difference structure factors and DED maps. It can be isotropic(`iso`), anisotropic(`aniso`), or no (`no`)scaling if 
        the structure factors are already scaled."""
        ),
        "nargs": 1,
        "required": False,
        "type": str,
        "default": "aniso",
    }),
    (("--atom_profile", "-v"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": "mean",
        "help": """One-dimensional atom-wise representation of difference density, based on `max` or `mean` peaks around the 
        atom profile. By default, it is mean."""
    }),
    (("--scaling_method", "-sm"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": "cctbx",
        "help": """Scaling method using methods implemented in CCTBX or CCP4-based programs."""
    }),
    (("--ccp4_scaling_method", "-cs"), {
        "nargs": 1,
        "metavar": None,
        "required": False,
        "type": str,
        "default": "linear",
        "help": """If the scaling method selected is based on CCP4 programs, linear scaling can be chosen by default."""
    }),
)

residem_pdb_mtz_args = (
    (("--ref", "-r"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """PDB ID to be downloaded in lowercase. Example: 5b6v"""
    }),
)

residem_coot_args = (
    (("--ref", "-r"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """Reference PDB for model visualization"""
    }),
    (("--map", "-m"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """Isomorphous difference density (DED) map for visualization"""
    }),
)

residem_extract_around_model_args = (
    (("--pdb", "-r"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """Reference PDB for mask calculation"""
    }),
    (("--map", "-m"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """Isomorphous difference density (DED) map for masking in CCP4 format"""
    }),
    (("--selection", "-s"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """Residue selection for masking the DED map. Phenix/CCTBX atom naming should be given as input.
        Example: 'resseq 182:214'"""
    }),
    (("--radius", "-d"), {
        "nargs": 1,
        "required": False,
        "type": float,
        "default": 5.0,
        "help": """Radius around the residue where the map should be truncated. Default is 5 Å."""
    }),
    (("--box", "-b"), {
        "nargs": 1,
        "required": False,
        "type": float,
        "default": 5.0,
        "help": """Box cushion around the residue where the map should be truncated. Default is 5 Å."""
    }),
    (("--output", "-o"), {
        "nargs": 1,
        "required": False,
        "type": str,
        "default": "output",
        "help": """Output name of the truncated map."""
    }),
)

residem_svd_csv_args = (
    (("--file", "-f"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """Files for analysis. List of files can be given, or if all files are in the same folder, "*.csv" 
        can be provided as input."""
    }),
    (("--cmap", "-c"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": "Spectral",
        "help": """Matplotlib colormap for plotting. Default is 'Spectral'. More can be found 
        at https://matplotlib.org/stable/users/explain/colors/colormaps.html"""
    }),
    (("--ncmap", "-s"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": "Blues",
        "help": """Network colormap for plotting. Default is 'Blues' but can be changed to any sequential colors 
        from https://matplotlib.org/stable/users/explain/colors/colormaps.html."""
    }),
    (("--closeness", "-k"), {
        "nargs": 1,
        "required": False,
        "type": float,
        "default": 0.05,
        "help": """Subplot node closeness. Default is 0.05, but can be adjusted depending on the layout."""
    }),
    (("--cutoff", "-i"), {
        "nargs": 1,
        "required": False,
        "type": float,
        "default": 0,
        "help": """Cutoff score for the data frame. The script either uses 'electron_sum_around_the_site' or 'peak_sum' 
        from the input CSV. A cutoff score can be provided to mask values above it. Default is zero."""
    }),
    (("--residue", "-n"), {
        "nargs": 1,
        "required": False,
        "type": str,
        "default": None,
        "help": """Residues for analysis and plotting network. Example: '182,216-218,300'"""
    }),
    (("--pdb", "-r"), {
        "nargs": 1,
        "required": False,
        "type": str,
        "default": None,
        "help": """Input reference PDB from which distance is calculated for network representation"""
    }),
    (("--dist", "-d"), {
        "nargs": 1,
        "required": False,
        "type": str,
        "default": None,
        "help": """Distance calculation type. Available options are CA (C-alpha), CM (center of mass for non-hydrogen atoms), 
        or SC (center of mass for non-hydrogen side chain atoms)."""
    }),
    (("--low", "-l"), {
        "nargs": 1,
        "required": False,
        "type": float,
        "default": 0.05,
        "help": """Low electron density sum cutoff threshold for plotting node size. Default is 0.05."""
    }),
    (("--medium", "-m"), {
        "nargs": 1,
        "required": False,
        "type": float,
        "default": 0.1,
        "help": """Medium electron density sum cutoff threshold for plotting node size."""
    }),
    (("--map_type", "-t"), {
        "nargs": 1,
        "required": False,
        "type": str,
        "default": "denisty",
        "help": """Type of map input. Options are 'peak' or 'density'. By default, 'density' is chosen."""
    }),
    (("--norm", "-z"), {
        "nargs": 1,
        "required": False,
        "type": bool,
        "default": True,
        "help": """Whether to normalize the input values. Must be either True or False."""
    }),
    (("--ew",), {
        "nargs": 1,
        "required": False,
        "type": float,
        "default": 0.150,
        "help": """Edge width between nodes. Default is 0.150."""
    }),
)

residem_pymol_image_args = (
    (("--pdb", "-r"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """PDB file for analysis. List of files can be given, or if all files are in the same folder, "*.csv" 
        can be provided as input."""
    }),
    (("--map", "-m"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """Isomorphous difference density map in CCP4 format."""
    }),
    (("--resi", "-resi"), {
        "nargs": 1,
        "required": True,
        "type": str,
        "default": None,
        "help": """List of residues for which map and image should be generated, using PyMOL atom selection syntax.
        Example: 'resi 182+216+300'. More details about [atom selection syntax](https://pymolwiki.org/index.php/Selection_Algebra)"""
    }),
    (("--all", "-all"), {
        "nargs": 1,
        "required": False,
        "type": bool,
        "default": True,
        "help": """If true, images of both positive and negative DED maps will be generated."""
    }),
    (("--positive", "-p"), {
        "nargs": 1,
        "required": False,
        "type": bool,
        "default": False,
        "help": """If true, only an image of the positive DED map will be generated."""
    }),
    (("--negative", "-n"), {
        "nargs": 1,
        "required": False,
        "type": bool,
        "default": False,
        "help": """If true, only an image of the negative DED map will be generated."""
    }),
    (("--sigma", "-s"), {
        "nargs": 1,
        "required": False,
        "type": float,
        "default": 3.0,
        "help": """Sigma cutoff of the DED map. Default is +/- 3.0 sigma."""
    }),
    (("--carve", "-c"), {
        "nargs": 1,
        "required": False,
        "type": float,
        "default": 1.2,
        "help": """Carve radius around the atom. Default is 1.2 Å."""
    }),
    (("--output", "-o"), {
        "nargs": 1,
        "required": False,
        "type": str,
        "default": "output",
        "help": """Output name of the PNG file."""
    }),
    (("--rot_x", "-x"), {
        "nargs": 1,
        "required": False,
        "type": int,
        "default": 0,
        "help": """Rotation around the x-axis in degrees."""
    }),
    (("--rot_y", "-y"), {
        "nargs": 1,
        "required": False,
        "type": int,
        "default": 0,
        "help": """Rotation around the y-axis in degrees."""
    }),
    (("--rot_z", "-z"), {
        "nargs": 1,
        "required": False,
        "type": int,
        "default": 0,
        "help": """Rotation around the z-axis in degrees."""
    }),
)

residem_parser = argparse.ArgumentParser(description=residem_overview)
residem_pdb_mtz_parser = argparse.ArgumentParser(description=residem_pdb_mtz_overview)
residem_coot_parser = argparse.ArgumentParser(description=residem_coot_overview)
residem_extract_around_model_parser = argparse.ArgumentParser(description=residem_extract_around_model_overview)
residem_svd_csv_parser = argparse.ArgumentParser(description=residem_svd_csv_overview)
residem_pymol_image_parser = argparse.ArgumentParser(description=residem_pymol_image_overview)

for args, kwargs in residem_args:
    residem_parser.add_argument(*args, **kwargs)

for args, kwargs in residem_pdb_mtz_args:
    residem_pdb_mtz_parser.add_argument(*args, **kwargs)

for args, kwargs in residem_coot_args:
    residem_coot_parser.add_argument(*args, **kwargs)

for args, kwargs in residem_extract_around_model_args:
    residem_extract_around_model_parser.add_argument(*args, **kwargs)

for args, kwargs in residem_svd_csv_args:
    residem_svd_csv_parser.add_argument(*args, **kwargs)

for args, kwargs in residem_pymol_image_args:
    residem_pymol_image_parser.add_argument(*args, **kwargs)
