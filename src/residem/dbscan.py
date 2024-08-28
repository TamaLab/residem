import numpy as np
import sys
import psutil
from multiprocessing import Pool


class DBSCAN:
    def __init__(self, coord, epsilon, minpts, idx=(0, 0, 0)):
        unique_arr = np.unique(coord, axis=0) # only keeps unique points
        self.coord = unique_arr
        self.n = len(unique_arr) # get the value of n unique points
        self.dist = self.equili_dis(False) # calculate distance matrix between points.
        self.visited = np.full((self.n), False) # sets all visted points to False intially
        self.noise = np.full((self.n), False) # sets all noise points to False intially
        self.epsilon = epsilon # The maximum distance between two points to be considered neighbors 0.75 set here
        self.minpts = minpts # minimum number of points 5
        self.idx = np.full((self.n), 0) # array that stores the cluster index for each point, initialize with zero
        self.C = 0 # A counter to keep track of the current cluster ID, initialize with zero
        self.idx1 = idx

    def equili_dis(self, single=False):
        """ This method computes the Euclidean distance between points.
        It uses a meshgrid to efficiently compute the distance matrix between all pairs of points."""

        p, q = np.meshgrid(np.arange(self.n), np.arange(self.n))
        if single:
            dist = np.sqrt(np.sum(((self.coord[p] - self.idx1) ** 2), 2))[0]
        else:
            p, q = np.meshgrid(np.arange(self.n), np.arange(self.n))
            dist = np.sqrt(np.sum(((self.coord[p] - self.coord[q]) ** 2), 2))
        return dist


    def run(self):
        for i in range(len(self.coord)):
            if self.visited[i] == False:
                self.visited[i] = True
                self.neighbors = self.regionQuery(i)
                if len(self.neighbors) >= self.minpts:
                    self.C += 1
                    self.expandCluster(i)
                else:
                    self.noise[i] = True
        return self.idx, self.noise

    def regionQuery(self, i):
        g = self.dist[i, :] < self.epsilon
        Neighbors = np.where(g)[0].tolist()
        return Neighbors

    def expandCluster(self, i):
        self.idx[i] = self.C
        k = 0

        while True:
            if len(self.neighbors) <= k: return
            j = self.neighbors[k]
            if self.visited[j] != True:
                self.visited[j] = True

                self.neighbors2 = self.regionQuery(j)
                v = [self.neighbors2[i] for i in np.where(self.idx[self.neighbors2] == 0)[0]]

                if len(self.neighbors2) >= self.minpts:
                    self.neighbors = self.neighbors + v

            if self.idx[j] == 0: self.idx[j] = self.C
            k += 1

    def sort(self):

        cnum = np.max(self.idx)
        self.cluster = []
        self.noise = []
        for i in range(cnum):
            k = np.where(self.idx == (i + 1))[0].tolist()
            self.cluster.append([self.coord[k, :]])

        self.noise = self.coord[np.where(self.idx == 0)[0].tolist(), :]
        return self.cluster, self.noise


class EquiDist(DBSCAN):
    def __init__(self, array, output, idx=(0, 0, 0)):
        super().__init__(array, idx)
        self.output = output

    def save_dist(self):
        np.savetxt(self.output + '.txt', self.equili_dis(True))


def apply_dbscan(coord, epsilon, minpts):
    dbscan = DBSCAN(coord, epsilon, minpts)
    dbscan.run()
    g_cluster, _ = dbscan.sort()
    return g_cluster

def chunked_dbscan(coord, epsilon, minpts, available_memory):
    total_size = len(coord)
    available_memory = psutil.virtual_memory().available
    proportion = min(0.1, max(0.01, available_memory / (coord.itemsize * coord.shape[1] * len(coord) * 10)))
    # Calculate the estimated chunk size based on available memory and total size
    estimated_chunk_size = int(total_size * proportion)
    chunk_size = min(total_size, estimated_chunk_size)

    chunks = [coord[i:i + chunk_size] for i in range(0, total_size, chunk_size)]
    cluster_results = []

    for chunk in chunks:
        g_cluster = apply_dbscan(chunk, epsilon, minpts)
        cluster_results.extend(g_cluster)

    return cluster_results



# def main(map_grid_file):
#     coord = np.loadtxt(map_grid_file)
#     dbscan = DBSCAN(coord, 0.7, 5)
#     dbscan.run()
#     g_cluster, ncluster = dbscan.sort()
#     return g_cluster

def main(map_grid_file):
    coord = np.loadtxt(map_grid_file)

    epsilon = 0.75
    minpts = 5

    # Estimate available memory
    available_memory = psutil.virtual_memory().available

    # Calculate the proportion based on available memory
    proportion = min(0.1, max(0.01, available_memory / (coord.itemsize * coord.shape[1] * len(coord) * 10)))

    try:
        g_cluster = apply_dbscan(coord, epsilon, minpts)
    except MemoryError:
        g_cluster = chunked_dbscan(coord, epsilon, minpts, available_memory)

    return g_cluster


def main_2(g_cluster, output):
    for ran in range(len(g_cluster)):
        with open(output + ".txt", 'a+') as f:
            np.savetxt(f, g_cluster[ran][0], fmt='%.3f')
            f.write(f"\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="array file")
    parser.add_argument("output", help="output dir")
    args = parser.parse_args()
    g_cluster = main(args.input, args.output)
    main_2(g_cluster, args.output)
