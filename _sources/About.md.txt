# About ResiDEM






Time-resolved serial femtosecond crystallography (TR-SFX) is a powerful technique that enables the observation of 
time-evolved dynamics of residual motion across crystal structures, providing insight into the complex structural
changes that occur during dynamic processes in biological systems. The information gathered from TR-SFX has significant 
implications in various fields, including structural biology, biochemistry, and drug discovery, as it sheds light 
on the intricate structural changes that occur during dynamic processes in biological systems. However, interpreting 
the results of TR-SFX can be challenging, because it can be difficult to differentiate signals
from noise while analyzing electron density maps. In this regard, new tools are required to elucidate the residual 
contributions involved in the dynamics. ResiDEM is one such tool that can be used to analyze and visualize the electron
density from TR-SFX with a focus on identifying structural changes during the reaction dynamics. This tool can generate a
three-dimensional difference electron density map using various weighting approaches. Additionally, it employs a 
clustering approach for analyzing difference electron density maps in order to identify the key residues involved in the
reaction dynamics of the target protein. The tool also provides a network representation for monitoring how the electron sum
changes with various time intermediates. In conclusion, this tool helps researchers gain a deeper understanding of the
molecular mechanisms involved in short-lived reaction intermediates captured by time-resolved TR-SFX data.


ResiDEM is a Python tool that uses various functions from the cctbx toolbox. It is designed to analyze isomorphous
difference densities from TR-SFX data through the data processing of structure factors. It generates difference electron
density maps with different weight schemes to reduce noise, enabling the analysis of changes in the electron density 
distribution within protein structures at different time points, which helps us to discover structural changes and protein
movements. This tool performs clustering and feature extraction using the DBSCAN clustering method, which
identifies and groups regions with similar changes in the electron density. This helps identify the most relevant features
and motions within the protein structure. Finally, the extracted features were visualized and interpreted to comprehend the
complex structural changes occurring within the protein over time.

For more information, kindly please see the [published paper](https://doi.org/10.1021/acs.jcim.4c00858). 

If the tool was useful please cite us: 

Raghavan, S. S. & Miyashita, O. [ResiDEM: Analytical Tool for Isomorphous Difference Electron Density Maps Utilizing Dynamic Residue Identification via Density Clustering. Journal of Chemical Information and Modeling (2024). DOI: 10.1021/acs.jcim.4c00858](https://doi.org/10.1021/acs.jcim.4c00858).


For any help or queries with respect to the tool, contact [me](mailto:hypowergravity@gmail.com).

<script>
document.querySelectorAll('a[href^="http"]').forEach(link => {
    link.setAttribute('target', '_blank');
});
</script>
