# Installation of ResiDEM

```{toctree}
    :titleonly:
    :maxdepth: 1
    :caption: "Contents:"
    :hidden:
        self

```

Installation of ResiDEM requires few dependency such as [Computational Crystallography Toolbox (CCTBX)](https://cctbx.github.io/).
It also some requires some basic python packages for data processing such as numpy, pandas, scipy etc...

## Installation

ResiDEM is a python tool written in python 3, and it is recommended that a new environment is created for better package
management and to avoid any dependency conflicts.

1. Creating a conda virtual environment and installation :

```bash
# you can change the name if you want
conda create -n residem python=3.10
conda activate residem

```

2. Updating the environment

```bash
# After activation of residem environment, we can check the location with command `which python`
# do the following to install in particular residem environment

full_path=$(which python)
base_path=$(dirname $(dirname "$full_path"))
echo $base_path

conda env update --prefix $base_path --file environment.yml --prune
```

2. Install ResiDEM package

```bash
pip install dist/residem-0.1.0.tar.gz

```

```{admonition} Additional Information
:class: note

The ResiDEM package might also need [CCP4 scaleit](https://www.ccp4.ac.uk/html/scaleit.html). It
will be used for scaling between reference and triggered state if that particular user option in given.
It also has inbuilt scaling method, but works well with
isomorphous data sets.
```

To check if the CCP4 package or scaleit in the path the following command can used.

```bash
if command -v ccp4 &> /dev/null
then
    echo "CCP4 in command line is Present"
else
    echo "CCP4 in command line is Not Present please install and give the path"
fi
```

## Unit test

Testing the tool with Bacteriorhodopsin data.
There are two python scripts which can be used to run to reproduce images as in the [published article](https://doi.org/10.1021/acs.jcim.4c00858)
.
The name of the scripts are [residem_unit_test.py](residem_unit_test.py) and [SVD_unit_test.py](SVD_unit_test.py).
The main script [residem_unit_test.py](residem_unit_test.py) can be used to test the tool and reproduce the certain figures as in [published article](https://doi.org/10.1021/acs.jcim.4c00858)
.
This script computes the isomorphous difference density for Bacteriorhodopsin(bR) for 13 time delays as published by [Nango et al.](https://www.science.org/doi/10.1126/science.aah3497)
.
Single computation approximately takes around 3-5 minutes in personal laptop. The script [residem_unit_test.py](residem_unit_test.py)
may take some time (~1 hour) to compute for 13 datasets.

The testing can be done using the following commands.

```bash
conda activate residem
python reside_unit_test.py
# after completeion of the above the follwoing can be submitted.
python SVD_unit_test.py

```

After running the script, image corresponding to Figure 3, Figure 4 and Figure 6(a) as in [published article](https://doi.org/10.1021/acs.jcim.4c00858)
can be reproduced.
Ater successful completion, the image can be found in the following places.

- Figure 3(a) can be seen in the file `SVD/760ns/one_dimensional_linear_plot.pdf` .
- Figure 3(b) can be seen in the file `SVD/negative/one_dimensional_linear_plot.pdf`.
- Figure 3(c) can be seen in the file `SVD/negative/SVD_original.pdf`.
- Figure 4 can be seen in the file `SVD/760ns/Difference_density_network_CA_nolabel.pdf`. In the generated image,
  the edge and nodes will be same, but the layout might be different as random layout is generated every time.
- Figure 6(a) can be seen in the file `SVD/all/correlation_plot.pdf`.

<script>
document.querySelectorAll('a[href^="http"]').forEach(link => {
    link.setAttribute('target', '_blank');
});
</script>
