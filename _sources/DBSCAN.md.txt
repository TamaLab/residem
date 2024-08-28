# Work Flow of DBSCAN implemented in ResiDEM 



```{toctree}
    :titleonly:
    :maxdepth: 1
    :caption: "Contents:"
    :hidden:
        self

```
## Explanation

Initialization Steps involves screening coordinates Input data points which corresponds to voxel with difference peak greater than +/- 3.0 sigma value.
Epsilon, Defines the radius within which neighbors are considered and is set 0.77.
Minpts, Minimum number of points required to form a cluster is set to 5. 
After which unique data points are obtained with number of point n estimated. 
Distance matrix between all points are calculated. Then initial parameters such as 
visited (Boolean array to track visited points), noise (Boolean array to mark noise points), 
idx:( Array to store cluster indices for each point) and 
C (Counter to keep track of the current cluster ID are set to False).



```{mermaid}
flowchart TD
    
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
If not visited, it marks the point as visited, finds its neighbors, and determines whether to start a new cluster or mark the point as noise.
If starting a new cluster, it expands the cluster by checking neighbors.

```{mermaid}
flowchart TD
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
If a neighbor has enough points around it, it gets added to the cluster, and the process repeats until no more neighbors can be added.

