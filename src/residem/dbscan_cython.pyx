# cython: boundscheck=False, wraparound=False, nonecheck=False
import numpy as np
cimport numpy as cnp
import sys

cdef class DBSCAN:
    cdef:
        cnp.ndarray coord
        int n
        cnp.ndarray dist
        cnp.ndarray visited
        cnp.ndarray noise
        float epsilon
        int minpts
        cnp.ndarray idx
        int C
        cnp.ndarray idx1
        list neighbors

    def __init__(self, coord, epsilon, minpts, idx=(0, 0, 0)):
        self.coord = coord
        self.n = len(coord)
        self.dist = self.equili_dis(False)
        self.visited = np.full((self.n), False)
        self.noise = np.full((self.n), False)
        self.epsilon = epsilon
        self.minpts = minpts
        self.idx = np.full((self.n), 0)
        self.C = 0
        self.idx1 = idx

    cpdef equili_dis(self, bint single=False):
        cdef:
            cnp.ndarray p, q
        p, q = np.meshgrid(np.arange(self.n), np.arange(self.n))
        if single:
            dist = np.sqrt(np.sum(((self.coord[p] - self.idx1) ** 2), 2))[0]
        else:
            dist = np.sqrt(np.sum(((self.coord[p] - self.coord[q]) ** 2), 2))
        return dist

    cpdef run(self):
        cdef int i
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

    cpdef regionQuery(self, int i):
        cdef cnp.ndarray g
        g = self.dist[i, :] < self.epsilon
        Neighbors = np.where(g)[0].tolist()
        return Neighbors

    cpdef expandCluster(self, int i):
        cdef int k, j
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

    cpdef sort(self):
    cdef:
        int cnum, i
        list cluster, noise
        cnp.ndarray k

    cnum = np.max(self.idx)
    cluster = []
    noise = []

    for i in range(cnum):
        k = np.where(self.idx == (i + 1))[0].tolist()
        cluster.append([self.coord[k, :]])

    noise = self.coord[np.where(self.idx == 0)[0].tolist(), :]
    return cluster, noise
