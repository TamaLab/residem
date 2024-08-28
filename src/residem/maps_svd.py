from residem.multiprocess_func import make_parallel
import os
from collections import defaultdict
import pathlib
import numpy as np
import tqdm
from residem.MapCoeffSet import MapCofSet
from cctbx.array_family import flex
import tensorflow as tf


# class SvdMap:
#     def __init__(self, array, size=4):
#         self.array = array
#         self.size = size if array.shape[0] > size else array.shape[0]
#         self.A1 = np.asarray([A[x].flatten() for x in range(A.shape[0])])
#         self.X, self.U, self.S, self.V = self.calc_svd(self.A1)
#
#     @staticmethod
#     def svd(X):
#         U, s, V = tf.linalg.svd(tf.cast(X, tf.float32))
#         # U, S, Vt = np.linalg.svd(X, full_matrices=True)
#         # S = np.zeros(X.shape)
#         # np.fill_diagonal(S, s)
#         return U, s, V
#
#     def calc_svd(self, X):
#         U, S, V = self.svd(X)
#         return X, U, S, V
#
#     def svd_arrays(self):
#         self.map_svd = defaultdict(int)
#         for j in range(self.size):
#             # self.map_svd[f"X_%s" % j] = self.S[j, j] * self.U[:, j][:, None] @ self.V[:, j][None, :]
#             self.map_svd[f"X_%s" % j] = np.array(U[j]).reshape(1, 1) @ np.array(s[j, j]).reshape(1, 1) @ np.array(
#                 V.numpy().transpose()[j]).reshape(1, V.shape[0])
#             print(j, self.map_svd[f"X_%s" % j].size())


class generate_array:
    def __init__(self, pdb, residue, map_files, radius=5.0, box_cushion=5.0):
        self.pdb = pdb
        self.residue = residue
        self.radius = radius
        self.box_cushion = box_cushion
        self.map = map_files
        self.check_space_group, self.file_name = self.make_arrays()
        self.map_coeff_map = self.diff_map[self.file_name[0]].map_cof_2

        self.array = self.check_sg_stack()
        np.save("data_svd.npy", self.array)

    def check_sg_stack(self):
        if len(set(self.check_space_group)) == 1:
            #stacked_array = np.array(np.asarray(list(self.diff_map_array.values())), dtype="float32")
            stacked_array = np.array(list(self.diff_map_array.values()), dtype="float32")

        else:
            print("The space group is not matching for all the given input maps ")
            stacked_array = None
        return stacked_array

    def make_arrays(self):
        self.diff_map = defaultdict(int)
        self.diff_map_array = defaultdict(int)
        map_name = []

        def func(x):
            name = pathlib.Path(x).stem
            map_name.append(name)
            self.diff_map["%s" % name] = MapCofSet(self.pdb, x, selection=self.residue, radius=self.radius, box_cushion=self.box_cushion)
            self.diff_map_array["%s" % name] = self.diff_map["%s" % name].map_cof_2.map_data().as_numpy_array()
            return self.diff_map["%s" % name].map_model.crystal_symmetry().space_group_number()

        space_group_number = make_parallel(func)(tqdm.tqdm(self.map))
        return space_group_number, map_name




def main():
    import argparse
    parser = argparse.ArgumentParser(prog='residem_svd_map', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f"""It performs map based SVD \n
        1. residem_svd_map -r 5b6v.pdb -f *.ccp4 \n
        2. residem_svd_map -r 5b6v.pdb -f *.ccp4 -s 2 
        2. residem_svd_map -r 5b6v.pdb -f *.ccp4 -s 2 -resi 'resseq 182:214' \n 
        The residue selection syntax is similar to that of the phenix for reference\n
        please visit the website: https://phenix-online.org/documentation/reference/atom_selections.html
        3. residem_svd_map -r 5b6v.pdb -f *.ccp4 -s 2 -resi 'resseq 182:214'  -b 3.0 -d 3.0\n
        4. residem_svd_map -r 5b6v.pdb -f *.ccp4 -s 2 -resi 'resseq 182 or resseq 214:216 or resseq 300'


    """)
    parser.add_argument("-r","--pdb", help="pdb_file")
    parser.add_argument("-resi","--residue", nargs='?', action='store', type=str, help="residue selection")
    parser.add_argument("-s","--size", action='store', type=int, help="number of the left singular vectors to be "
                                                                      "written ")
    parser.add_argument("-f", "--file", action='store', nargs='+',
                        help="map files for SVD calculations")
    parser.add_argument("-b", "--box", action='store',type=float,
                        help="radius of box cushion")
    parser.add_argument("-d", "--radius", action='store',type=float,
                        help="radius buffer around the chosen residue")


    args = parser.parse_args()
    data_files = [x for x in args.file]
    print("The given map files are ", data_files)
    residue = "all" if args.residue is None else args.residue
    box = args.box if args.box is not None else 5.0
    radius = args.radius if args.radius is not None else 5.0
    SVD_map = generate_array(args.pdb, residue, data_files,radius=radius, box_cushion=box)
    size = args.size if not None and args.size < len(data_files) else len(data_files)
    if SVD_map.array is not None:
        Array = SVD_map.array.reshape(len(data_files), -1).T
        s,U,V = tf.linalg.svd(tf.cast(Array, tf.float32))
        Array_reconstituted = np.array(U[:, 0].numpy().reshape(-1, 1) @ s[0].numpy().reshape(1, 1)) @ np.array(V[0, :]).reshape(1, -1)
        flex_reconstituted = flex.double(np.ascontiguousarray(Array_reconstituted.T[0]).reshape(SVD_map.diff_map_array[SVD_map.file_name[0]].shape))
        SVD_map.diff_map[SVD_map.file_name[0]].map_cof_2.map_manager().customized_copy(map_data = flex_reconstituted).write_map(
                f"SVD_map_reconstituted.ccp4")

        # get the lsv and the SVD of principle reconstructed vector.
        for x in range(size):
            flex_array = flex.double(np.ascontiguousarray(U[:,x]).reshape(SVD_map.diff_map_array[SVD_map.file_name[0]].shape))
            SVD_map.diff_map[SVD_map.file_name[0]].map_cof_2.map_manager().customized_copy(map_data = flex_array).write_map(
                f"SVD_map_reconstituted_%s_lsv.ccp4"%x)
            


if __name__ == '__main__':
    main()
