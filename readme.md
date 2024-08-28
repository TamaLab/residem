
[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3100/)
[![CI](https://github.com/hypowergravity/test_residem/actions/workflows/ci.yml/badge.svg)](https://github.com/hypowergravity/test_residem/actions/workflows/ci.yml)


Platform | CI Status
---------|:---------
OSX      | ![macOS Build Status](https://github.com/hypowergravity/test_residem/actions/workflows/ci.yml/badge.svg?branch=main&event=push&job=test_macos)
Linux    | ![Ubuntu Build Status](https://github.com/hypowergravity/test_residem/actions/workflows/ci.yml/badge.svg?branch=main&event=push&job=test_ubuntu)




# Table of Contents

1. [About](#About)
2. [Installation](#Installation)
3. [License](#License)
4. [Contact](#Contact)

---

# Resi-DEM.

Resi-DEM is a tool that identifies residues in isomorphous difference density map captured using Time-resolved serial femtosecond crystallographic data. It uses DBSCAN based clustering algorithm to cluster the co-ordinates of the difference peak in the difference density map and assigns the cluster to the neighbouring atoms. It can also perform various calculations and plot their outputs. A breif documentation of ResiDEM can be found [here](https://hypowergravity.github.io/test_residem/) 

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
Residem is an open-source tool available under the MIT License, providing users with the freedom to use, modify, and distribute the software. Additionally, Residem incorporates adapted modules from the Computational Crystallography Toolbox (cctbx), and users are advised to also refer to the cctbx license to comply with its terms and conditions. 

# Contact
* Sriram Srinivasa Raghavan : hypowergravity@gmail.com, sriram.srinivasaraghavan@riken.jp
* Osamu Miyashita : osamu.miyashita@riken.jp
