
[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3100/)
[![CI](https://github.com/TamaLab/residem/actions/workflows/ci.yml/badge.svg)](https://github.com/TamaLab/residem/actions/workflows/ci.yml)


Platform | CI Status
---------|:---------
OSX      | ![macOS Build Status](https://github.com/TamaLab/residem/actions/workflows/ci.yml/badge.svg?branch=main&event=push&job=test_macos)
Linux    | ![Ubuntu Build Status](https://github.com/TamaLab/residem/actions/workflows/ci.yml/badge.svg?branch=main&event=push&job=test_ubuntu)



# Table of Contents

1. [About](#About)
2. [Installation](#Installation)
3. [License](#License)
4. [Contact](#Contact)

---

# ResiDEM.

Resi-DEM is a tool that identifies residues in isomorphous difference density map captured using Time-resolved serial femtosecond crystallographic data. It uses DBSCAN based clustering algorithm to cluster the co-ordinates of the difference peak in the difference density map and assigns the cluster to the neighbouring atoms. It can also perform various calculations and plot their outputs. A breif documentation of ResiDEM can be found [here](https://tamalab.github.io/residem/).

If ResiDEM was useful for your analysis, please cite:
Raghavan, S. S. & Miyashita, O. ResiDEM: Analytical Tool for Isomorphous Difference Electron Density Maps Utilizing Dynamic Residue Identification via Density Clustering. Journal of Chemical Information and Modeling (2024). [DOI: 10.1021/acs.jcim.4c00858](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00858).

## Installation
Quick installation of ResiDEM can be done using the following:

1. Creating a conda virtual environemt:

```bash
# you can change the name if you want
conda create -n residem python=3.10
```

2. Activation of conda environment

```bash
conda activate residem
```

3. Updating the environment.
Environment file can be found [here](environment.yml).

```bash
# After activation of residem environment, we can check the location with the command `which python`
# do the following to install in particular residem environemnt

full_path=$(which python)
base_path=$(dirname $(dirname "$full_path"))
echo $base_path

conda env update --prefix $base_path --file environment.yml --prune

```

4. To uninstall the package :

```
pip uninstall residem
```




# License
Residem is an open-source tool available under the MIT License, providing users with the freedom to use, modify, and distribute the software. Additionally, Residem incorporates adapted modules from the Computational Crystallography Toolbox (cctbx), and users are advised to also refer to the cctbx license to comply with its terms and conditions. For our implementation of weighting schemes, we referred to Xtrapol8, a tool licensed under the MIT License. The original license text can be found in the LICENSE file.

# Acknowledgement

In the development of the weighting schemes in this software, we referred to the following works and recommend citing these references when relevant:
* De Zitter, E., Coquelle, N., Oeser, P., Barends, T. & Colletier, J.-P. Xtrapol8 enables automatic elucidation of low-occupancy intermediate states in crystallographic studies. Communications Biology 5 (2022). DOI: 10.1038/s42003-022-03575-7. 
* Ursby, T. & Bourgeois, D. Improved estimation of structure-factor difference amplitudes from poorly accurate data. Acta Crystallographica Section A: Foundations of Crystallography 53, 564–575 (1997). DOI: 10.1107/S0108767397004522.
* Ren, Z. et al. A molecular movie at 1.8 Å resolution displays the photocycle of photoactive yellow protein, a eubacterial blue-light receptor, from nanoseconds to seconds. Biochemistry 40, 13788–13801 (2001). DOI: 10.1021/bi0107142.
* Hekstra, D. R. et al. Electric-field-stimulated protein mechanics. Nature 540, 400–405 (2016). DOI:  10.1038/nature20571.
* Pandey, S. et al. Time-resolved serial femtosecond crystallography at the European XFEL. Nature Methods 17, 73–78 (2019). DOI: 10.1038/s41592-019-0628-z.

Additionally, this software uses CCTBX. Please refer to the CCTBX license documentation for further details on their licensing terms and conditions: https://cctbx.github.io.




# Contact
* Sriram Srinivasa Raghavan : hypowergravity@gmail.com, sriram.srinivasaraghavan@riken.jp
* Osamu Miyashita : osamu.miyashita@riken.jp
