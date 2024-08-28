# Limitations

```{toctree}
    :titleonly:
    :maxdepth: 1
    :caption: "Contents:"
    :hidden:
        self

```

Isomorphous difference density (DED)maps can be quite [noisy](https://doi.org/10.1063/4.0000196). 
Different weighting schemes are applied to mitigate and supress the noise in DED maps. 
If the difference density (blobs) are continuous and stretch across several residues, then our tool only assigns such
continuous blob to one particular residue. In such cases, the sigma cutoff can be increased and clustering can be redone. 
Even manual interpretation of such densities are complex, which might be associated to large movement of the residues. 
Atom wise representation might then be useful. Generally, proteins have large amount of atoms and such large movements can 
inherently introduce complexity in DED analysis.

The main assumption in DED analysis is that the change in phase {math}`\Delta \phi_{|\Delta F|,\# }` of the triggered system is [small](https://doi.org/10.1063/4.0000196) as
compared to the reference state. In such cases the analysis might be useful. 

















<script>
document.querySelectorAll('a[href^="http"]').forEach(link => {
    link.setAttribute('target', '_blank');
});
</script>