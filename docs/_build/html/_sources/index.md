---
title: ResiDEM, How to start, without being mystified 🧙
subject: Tutorial
subtitle: Analyze the isomorphous difference density map.
short_title: How to ResiDEM
authors:
  - name: Sriram Srinivasa Raghavan
    affiliations:
      - RIKEN Center for Computational Science, Kobe, Hyogo 650-0047, Japan

    orcid: 0000-0003-3781-3921
    email: hypowergravity@gmail.com
  - name: Osamu Miyashita
    affiliations:
      - RIKEN Center for Computational Science, Kobe, Hyogo 650-0047, Japan

    orcid: 000-0002-2530-5674
    email: osamu.miyashita@riken.jp
license: MIT
keywords: Isomorphous difference map analysis, Network representation, Pearson correlation
---

# Welcome to ResiDEM Documentation.
```{toctree}
    :titleonly:
    :maxdepth: 1
    :caption: "Contents:"
    :hidden:
        self

Installation
Quickstart
Algorithm
Result
command
Limitation
About
```
ResiDEM is a tool to analyze isomorphous difference density map and
estimate residues having difference density feature associated with it.

This documentation provides the Command line functionality of running the ResiDEM in Linux and Mac OS.
The code might still be used in Windows, but the code is not tested in Windows.

The [Installation](Installation.md) contains the installation instructions of ResiDEM command line tool. 
Sample usage of basic command line definition can be seen in [Quick Start Guide](Quickstart.md). 
If you want to know about the [algorithm](Algorithm) or more [About](About.md) ResiDEM, you can check out our paper.

[ResiDEM: Analytical Tool for Isomorphous Difference Electron Density Maps Utilizing Dynamic Residue Identification via Density
Clustering](https://doi.org/10.1021/acs.jcim.4c00858)


[//]: # (This software was conceived in [Tama's Lab]&#40;https://sites.google.com/view/computationalbiophysicslab/&#41;)

[//]: # ( at )

[//]: # ( RIKEN Center for computational Science&#40;RCCS&#41;  by [Sriram Srinivasa Raghavan]&#40;&#41;.)

```{admonition} The Main **ResiDEM** functions in jiffy!
:class: note
```
ResiDEM is composed of a few command line utilities. A detailed description of the [command line](command.md) utilities are found [here](command.md). 
The main function of the tool to identify the residues associated with difference density in isomorphous difference density map. 
It also quantifies the difference density in single and multiple time periods collected at different experimental time delays.


- ``residem`` is a command line utility which takes three inputs for calculation: Refined reference structure,
structure factor amplitudes or intensity profile of reference and triggered states. 
This command line utility computes the isomorphous difference density map and identifies the difference density associated
to particular residues. It generates results which can be further used for network visualization and other analyses. 

- ``residem_svd_csv`` is a command line utility which takes a reference state pdb and other command line arguments and 
generates plots associated with difference density. It plots 1D,2D representation of 
normalized quantified difference density which is further used for network representation. 

- ``residem_coot`` is a command line utility which requires two arguments: one that corresponds to reference state structure(pdb), 
and the other being the generated isomorphous difference density (DED) map. 

A detailed version of command line argument can be seen in [command line](command) section. 



<script>
document.querySelectorAll('a[href^="http"]').forEach(link => {
    link.setAttribute('target', '_blank');
});
</script>