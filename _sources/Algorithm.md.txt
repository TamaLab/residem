# Algorithm
```{toctree}
    :titleonly:
    :maxdepth: 1
    :caption: "Contents:"
    :hidden:
        self

```
The initial purpose of the tool is to analyse isomorphous difference density map (DED).  DED maps are generated from 
difference structure factor amplitudes, obtained by subtracting reference structure factor amplitudes of 
reference {math}`|F_{obs}^{ref}|` and triggered state {math}`|F_{obs}^{t}|`, measured at time t after reaction initiation.
Detailed description of DED maps can be found [elsewhere](https://doi.org/10.1063/4.0000196). 
The time-evolved difference density map describes the difference density profile as functions of time Δσ(r,t), 
where r denotes a discrete set of grid points (voxel) over the space of an asymmetric unit cell. 
In each difference map, the voxel value is subject to the dynamics of the structural intermediates. 
These difference maps are often inspected manually to identify the associated features. 
Here, we employed the Density-Based Spatial Clustering of Applications with Noise [(DBSCAN)](https://dl.acm.org/doi/10.5555/3001460.3001507) 
algorithm to separately group voxels with positive and negative values. Two sets of parameters are given for the `DBSCAN` algorithm: 
the maximum distance between two points that are considered to be part of the same cluster (eps) and the minimum set of
points to form a cluster or dense region (minPts). In our implementation, we set the eps to 0.73 Å, which is
the covalent radius of [carbon](https://link.springer.com/article/10.1007/BF00907656), 
Minimum points (minPts) are set to 5 voxels that can be grouped together.
The workflow of `DBSCAN` implemented in our tool is as given below. A clearer expanded version can be found [here](DBSCAN.md).

```{mermaid}
flowchart LR
    
    B[inputs] --> C{preprocess}
    C --> J{Are all points processed?}

    J -- No --> K{Visited?}
    K -- No --> L[Mark as visited]
    L --> M[Find neighbors within epsilon using distance matrix]
    M --> N{Neighbors >= minpts?}

    N -- Yes --> O[Start a new cluster, increment C]
    O --> P[Assign current point to cluster C]
    P --> Q[Expand Cluster]
    subgraph Workflow 2
    Q
    end
    Q --> J

    N -- No --> R[Mark as noise]
    R --> J

    K -- Yes --> J
    J -- Yes --> S[End]

```
    Workflow 1 of DBSCAN algorithm. 
    
    The main workflow iterates over each point to check if it has been visited.
    If not visited, it marks the point as visited, finds its neighbors, and 
    determines whether to start a new cluster or mark the point as noise.
    If starting a new cluster, it expands the cluster by checking neighbors.

```{mermaid}
flowchart LR
    subgraph Workflow 1
    T
    end
    T[Expand Cluster] --> T1{Are there unvisited neighbors?}
    T1 -- Yes --> T2[Mark neighbor as visited]
    T2 --> T3{Does neighbor have enough neighbors?}
    T3 -- Yes --> T4[Add neighbors to cluster C]
    T4 --> T5[Update cluster index to C]
    T5 --> T1
    T3 -- No --> T6[Continue checking neighbors]
    T6 --> T1
    T1 -- No --> T7[Finish expanding cluster C]

```

    Workflow 2 of DBSCAN algorithm for Expanding cluster.
    
    This workflow expands the cluster by recursively checking neighbors.
    If a neighbor has enough points around it, it gets added to the cluster, 
    and the process repeats until no more neighbors can be added.

The grid points within the protein are generally selected as the mask regions for cluster identification. 
After clustering, the sum of difference electron density values is calculated for each cluster by summing 
the value of the voxels within the cluster if the values are above the σ value (±3.0 σ as the default value) 
in the difference electron density map. Then, an atom that is closest to and within a 2 Å radius of the maximum peak 
within the cluster is identified. This atom, which is associated with the cluster, is further used to examine the 
density changes in individual structural intermediates. In addition, for each of the identified clustered regions, 
hereafter called electron density blobs, specific features such as their volumes and their electron sums are calculated. 
The blob volume is evaluated by providing the coordinate of the obtained cluster in “maptbx.peak_volume_estimate” function 
implemented in cctbx tool. Additionally, to visualize changes in difference electron density consistently throughout the time series data,
we normalized the values to set their range between 0 and 1, using the minimum and maximum values within the time series data.



In a previous study, [Wickstrand et al.](https://doi.org/10.1063/1.5126921) suggested an approach to associate atoms with 
different electron densities implemented in Maptool, which portrays the three-dimensional changes in electron density
as a one-dimensional representation plot associated with each residue.
Recently, [Dasgupta](https://github.com/asmit3/eden) published the same methodology as Maptool for estimating the average amplitude using the cctbx module. 
Our tool has also additionally implemented this method, where we estimate the average difference density
within a 2 Å radius of each atom for all residues in the protein. However, this approach may not provide information 
about the distribution surrounding the atoms because it is bound by a radius cutoff. The clustering approach was used 
to quantify the electron density features and delineate the residual contributions associated with difference densities 
that were not directly on top of the atoms. 






<script>
document.querySelectorAll('a[href^="http"]').forEach(link => {
    link.setAttribute('target', '_blank');
});
</script>











<script>
document.querySelectorAll('a[href^="http"]').forEach(link => {
    link.setAttribute('target', '_blank');
});
</script>