from scitbx.array_family import flex
import tqdm
from cctbx import sgtbx
from residem.cctbx_func import CctbxFunc
from residem.dbscan import DBSCAN
from cctbx import maptbx
import numpy as np
from scitbx import matrix
import pandas as pd
from iotbx.map_model_manager import map_model_manager
from mmtbx.find_peaks import peaks_holder
import warnings
from collections import defaultdict
import json
from residem.point_inside_polygon import PolyPoint
import os
from residem.multiprocess_func import make_parallel
import concurrent.futures
import os
from functools import wraps
warnings.simplefilter(action='ignore', category=FutureWarning)




def closest_near_peak(result, max_dist, min_dist):
    """selection of closest atom near a map object"""
    smallest_distances_sq = result.smallest_distances_sq
    in_box = smallest_distances_sq > 0
    not_too_far = smallest_distances_sq <= max_dist ** 2
    not_too_close = smallest_distances_sq >= min_dist ** 2
    selection = (not_too_far & not_too_close & in_box)
    iseqs_of_closest_atoms = result.i_seqs.select(selection)
    return iseqs_of_closest_atoms, selection


def trunc(values, decs=0):
    return np.trunc(values * 10 ** decs) / (10 ** decs)


class Peaks_coff:
    def __init__(self, map_model, resolution, sigma, path=None):
        self.path = path if not None else os.getcwd()
        self.map_model = map_model
        self.sigma = sigma
        self.high_resolution = resolution
        self.use_selection = self.map_model.model().get_xray_structure().all_selection()
        self.map_p = self.map_model.deep_copy()
        self.map_n = self.map_model.deep_copy()
        self.real_map_unpadded_p = self.map_p.map_as_fourier_coefficients(
            d_min=resolution).fft_map().apply_sigma_scaling().real_map_unpadded()
        self.real_map_unpadded_n = self.map_n.map_as_fourier_coefficients(
            d_min=resolution).fft_map().apply_sigma_scaling().real_map_unpadded() * -1
        self.uc = self.map_model.crystal_symmetry().unit_cell()
        self.grid_p = self.grid(self.real_map_unpadded_p, self.sigma)
        self.grid_n = self.grid(self.real_map_unpadded_n, self.sigma)
        self.selection = self.map_model.model().get_xray_structure().all_selection()
        self.site_frac_p = self.site_fraction(self.grid_p)
        self.site_frac_n = self.site_fraction(self.grid_n)
        self.site_cart_p = self.uc.orthogonalize(self.peak_near_model(self.real_map_unpadded_p, self.site_frac_p))
        self.site_cart_n = self.uc.orthogonalize(self.peak_near_model(self.real_map_unpadded_n, self.site_frac_n))
        # np.savetxt(f'{self.path}/array_n.txt', self.site_cart_n.as_numpy_array(), fmt='%0.4f')
        # np.savetxt(f'{self.path}/array_p.txt', self.site_cart_p.as_numpy_array(), fmt='%0.4f')
        self.cluster_p = self.DBscan(self.site_cart_p.as_numpy_array())
        self.cluster_n = self.DBscan(self.site_cart_n.as_numpy_array())
        print("Done Clustering")

        self.cluster_p_peak_volume = self.peak_volume(self.real_map_unpadded_p, self.cluster_p, sigma)
        self.cluster_n_peak_volume = self.peak_volume(self.real_map_unpadded_n, self.cluster_n, sigma)
        self.peak_height_p = self.peak_height(self.real_map_unpadded_p, self.cluster_p_peak_volume)
        self.peak_height_n = self.peak_height(self.real_map_unpadded_n, self.cluster_n_peak_volume)
        self.unit_call_volume = self.map_model.model().get_xray_structure().unit_cell().volume()
        self.cluster_p_peak_volume_df = pd.DataFrame(self.cluster_p_peak_volume)
        self.cluster_n_peak_volume_df = pd.DataFrame(self.cluster_n_peak_volume)
        self.cluster_n_peak_volume_df["peak_height"] = self.cluster_n_peak_volume_df["peak_height"] * -1
        self.cluster_p_site = self.peak_site_frac(self.cluster_p_peak_volume_df)
        self.cluster_n_site = self.peak_site_frac(self.cluster_n_peak_volume_df)
        self.cluster_positive_p = self.peak_list_cluster(self.cluster_p_site, self.cluster_p_peak_volume_df)
        self.cluster_positive_n = self.peak_list_cluster(self.cluster_n_site, self.cluster_n_peak_volume_df)
        self.cluster_p_list = self.cluster_list(self.cluster_positive_p)
        self.cluster_n_list = self.cluster_list(self.cluster_positive_n)
        self.cluster_density_p = self.cluster_density_add(self.cluster_p_list, self.cluster_p_peak_volume_df,
                                                          self.cluster_p)
        self.cluster_density_n = self.cluster_density_add(self.cluster_n_list, self.cluster_n_peak_volume_df,
                                                          self.cluster_n)
        self.cluster_positive_negative = self.cluster_combine(self.cluster_density_p, self.cluster_density_n)

    def __repr__(self):
        return "This object returns the peaks co-ordinates by dbscan clustering identified by the map file"

    def atom_with_in_cluster(self, string, cluster, num):
        atoms = self.map_model.model().get_hierarchy().select(self.map_model.model().selection(f"({string})"),
                                                              copy_atoms=True)
        point = np.around([atoms.xyz for atoms in atoms.atoms()][0], 3).tolist()
        cluster_point = np.around(cluster[num], 3).tolist()
        return "within" if PolyPoint(point, cluster_point).inside_polygon else "near"

    def cluster_density_add(self, res_list, df, cluster):
        test = pd.DataFrame(np.array(df.site_cart.to_list(), dtype=np.float32), columns=["X", "Y", "Z"])

        df_peak_list = pd.DataFrame(
            columns=["resn", "resi", "atom name", "X", "Y", "Z", "peak height", "number of voxel", 'Cluster number',
                     "Peak Volume of the Cluster", "electron_sum_around_the_site","peak_sum"])

        for i in range(len(res_list)):
            X, Y, Z = res_list["X"][i], res_list["Y"][i], res_list["Z"][i]
            number = test.loc[(test['X'] == X) & (test['Y'] == Y) & (test['Z'] == Z)].index[0]
            df_peak_list = df_peak_list.append({"resn": res_list["resn"].loc[i], "resi": res_list["resi"].loc[i],
                                                "atom name": res_list["atom name"].loc[i],
                                                "X": X, "Y": Y, "Z": Z,
                                                "peak height": float(res_list["peak height"].loc[i]),
                                                "number of voxel": df["number of voxel"].loc[number],
                                                "site": self.atom_with_in_cluster(f"resseq %s and name %s" % (
                                                    res_list["resn"].loc[i], res_list["atom name"].loc[i]), cluster,
                                                                                  number),
                                                'Cluster number': df["cluster_number"].loc[number],
                                                "Peak Volume of the Cluster": df["blob_volume"].loc[number],
                                                "electron_sum_around_the_site": df["electron_sum_around_the_site"].loc[
                                                    number],"peak_sum":df["peak_sum"].loc[number]},ignore_index=True)

            df_peak_list = df_peak_list.sort_values(by="electron_sum_around_the_site", ascending=False)
            df_peak_list['Period'] = df_peak_list.resn.str.cat(df_peak_list[["resi", "atom name"]],
                                                               sep='-').str.replace(' ', '')

        return df_peak_list

    def crystal_gridding(self):
        """crystal gridding based on the input resolution """
        mmm = self.map_model
        # mmm.box_all_maps_around_density_and_shift_origin()
        mc = self.map_model.map_as_fourier_coefficients(d_min=self.high_resolution)
        rm = mc.fft_map().apply_sigma_scaling().real_map_unpadded()

        symmetry_flags = sgtbx.search_symmetry_flags(
            use_space_group_symmetry=True,
            use_space_group_ltr=0,
            use_seminvariants=True,
            use_normalizer_k2l=False,
            use_normalizer_l2n=False)

        try:
            cg = maptbx.crystal_gridding(space_group_info=self.map_model.model().crystal_symmetry().space_group_info(),
                                         symmetry_flags=symmetry_flags,
                                         unit_cell=self.map_model.unit_cell_crystal_symmetry().unit_cell(),
                                         pre_determined_n_real=rm.accessor().all())

            cgt = maptbx.crystal_gridding_tags(gridding=cg)
        except RuntimeError:
            cgt = mc.crystal_gridding(d_min=self.high_resolution, symmetry_flags=symmetry_flags).tags()

        return mmm, rm, cgt

    def peak_list_cluster(self, site, df):
        use_selection = self.map_model.model().get_xray_structure().all_selection()
        result = self.map_model.model().get_xray_structure().closest_distances(sites_frac=site,
                                                                               distance_cutoff=2,
                                                                               use_selection=use_selection)

        iseqs_of_closest_atoms_near, selection_near = self.closest_near_peak(result, 2, 0)

        peaks_near = peaks_holder(heights=flex.double(df.peak_height).select(selection_near),
                                  sites=site.select(selection_near),
                                  iseqs_of_closest_atoms=iseqs_of_closest_atoms_near)

        return peaks_near

    def peak_site_frac(self, df):
        vol_site = self.uc.fractionalize(flex.vec3_double(df.site_cart))

        # for i in range(df.site_cart.size):
        #     vol_site.append(
        #         [m / n for m, n in zip((df.site_cart[i][0], df.site_cart[i][1], df.site_cart[i][2]),
        #                                self.real_map_unpadded_p.all())])
        return vol_site

    def peak_height(self, real_map, cluster_volume):
        cluster_volume["peak_height"] = []
        for n in range(len(cluster_volume["site_cart"])):
            peak_height = real_map.value_at_closest_grid_point(self.uc.fractionalize(cluster_volume["site_cart"][n]))
            cluster_volume["peak_height"].append(peak_height)
        return cluster_volume



    def grid(self, real_map_unpadded, cutoff):
        # grid_f = []
        # for i in range(real_map_unpadded.last()[0]):
        #     for j in range(real_map_unpadded.last()[1]):
        #         for k in range(real_map_unpadded.last()[2]):
        #             if 0 < cutoff < real_map_unpadded[i, j, k]:
        #                 grid_f.append((i, j, k))
        #             elif 0 > cutoff > real_map_unpadded[i, j, k]:
        #                 grid_f.append((i, j, k))
        last_element = real_map_unpadded.last()
        coordinates = [(i, j, k) for i in range(last_element[0]) for j in range(last_element[1]) for k in
                       range(last_element[2])]

        def process_element(coords):
            i, j, k = coords
            element = real_map_unpadded[i, j, k]
            if 0 < cutoff < real_map_unpadded[i, j, k]:
                return (i, j, k)
            elif 0 > cutoff > real_map_unpadded[i, j, k]:
                return (i, j, k)
            else:
                return

        grid_f = make_parallel(process_element)(coordinates)
        new_list = list(filter(lambda x: x is not None, grid_f))

        return new_list

    def site_fraction(self, grid):
        site_frac = flex.vec3_double()
        for i in range(len(grid)):
            site_frac.append(
                [m / n for m, n in zip((grid[i][0], grid[i][1], grid[i][2]), self.real_map_unpadded_p.all())])

        return site_frac

    @staticmethod
    def closest_near_peak(result, max_dist, min_dist=0):
        """selection of closest atom near a map object"""
        smallest_distances_sq = result.smallest_distances_sq
        in_box = smallest_distances_sq >= 0
        not_too_far = smallest_distances_sq <= max_dist ** 2
        not_too_close = smallest_distances_sq >= min_dist ** 2
        selection = (not_too_far & not_too_close & in_box)
        iseqs_of_closest_atoms = result.i_seqs.select(selection)
        return iseqs_of_closest_atoms, selection

    def peak_near_model(self, real_map, site_frac):

        result_near = self.map_model.model().get_xray_structure().closest_distances(sites_frac=site_frac,
                                                                                    distance_cutoff=3,
                                                                                    use_selection=self.selection)

        iseqs_of_closest_atoms_near, selection_near = self.closest_near_peak(result_near, 3, 0)
        peaks_h = flex.double(np.array([real_map.value_at_closest_grid_point(site_frac[n]) for n in
                                        range(len(site_frac))]))

        peaks_near = peaks_holder(
            heights=peaks_h.select(selection_near),
            sites=result_near.sites_frac.select(selection_near),
            iseqs_of_closest_atoms=iseqs_of_closest_atoms_near)
        return peaks_near.sites

    def cluster_list(self, peak_site):
        df_peak_list = pd.DataFrame(columns=["resn", "resi", "atom name", "X", "Y", "Z", "peak height"])
        for length in range(peak_site.sites.size()):
            for m in self.map_model.model().get_hierarchy().models():
                for chain in m.chains():
                    for residue_group in chain.residue_groups():
                        for conformers in residue_group.conformers():
                            for residue in conformers.residues():
                                for atoms in residue.atoms():
                                    if peak_site.iseqs_of_closest_atoms[length] == atoms.i_seq:
                                        peak = np.around(np.array(self.uc.orthogonalize(peak_site.sites[length])), 3)
                                        df_peak_list = df_peak_list.append({"resn": f"%s" % residue.resseq,
                                                                            "resi": f"%s" % str(chain.id + '-'+ residue.resname),
                                                                            "atom name": f"%s" % atoms.name,
                                                                            "X": peak[0],
                                                                            "Y": peak[1],
                                                                            "Z": peak[2],
                                                                            "peak height": f"%.3f" % peak_site.heights[
                                                                                length]},
                                                                           ignore_index=True)

        return df_peak_list

    def site_cluster(self, df):
        # calculate if the atom is inside polygon
        pass

    def DBscan(self, site_cart):
        """eps = 0.73 covalent radius of carbon,  number of cluster points"""
        db = DBSCAN(site_cart, 0.73, 3)
        db.run()
        g_cluster = db.sort()[0]
        return g_cluster

    def peak_volume(self, real_map_unpadded, cluster, cutoff):

        def result_cluster(data):
            result = abs(real_map_unpadded.value_at_closest_grid_point(data))
            return result

        def result_cluster_e_sum(data):
            result = abs(rho.value_at_closest_grid_point(data))
            return result

        voxel_volume = self.map_model.model().get_xray_structure().unit_cell().volume() / matrix.col(
            real_map_unpadded.all()).product()
        # scaling_factor = voxel_volume / self.map_model.model().get_xray_structure().f_000()


        # real_map_unpadded_all = self.map_p.map_as_fourier_coefficients(
        #     d_min=self.high_resolution).fft_map().real_map_unpadded()

        rho = (real_map_unpadded / (
                self.map_model.model().get_xray_structure().unit_cell().volume())) * self.map_model.model().get_xray_structure().f_000()

        cluster_number = []
        blob_volume = []
        site_cart = []
        cluster_size = []
        electron_sums_around_atoms = []
        peak_sum = []
        for x in range(len(cluster)):
            r = maptbx.peak_volume_estimate(
                map_data=real_map_unpadded,
                sites_cart=flex.vec3_double(cluster[x][0]),
                crystal_symmetry=self.map_model.crystal_symmetry(),
                cutoff=cutoff,
                atom_radius=1.5)
            cluster_number.append(x)
            cluster_size.append(len(cluster[x][0]))
            if r is not None:
                blob_volume.append(r)
            else:
                blob_volume.append(0)


            gias = maptbx.grid_indices_around_sites(
                unit_cell=self.map_model.model().get_xray_structure().unit_cell(),
                fft_n_real=real_map_unpadded.all(),
                fft_m_real=real_map_unpadded.all(),
                sites_cart=flex.vec3_double(cluster[x][0]),
                site_radii=flex.double(len(flex.vec3_double(cluster[x][0])), 2))

            negative_cluster = flex.vec3_double(np.concatenate(cluster[x], axis=0))
            # rho = maptbx.eight_point_interpolation(map_grid, negative_cluster) * scaling_factor
            #
            # electron_sums_around_atoms.append("%.2f" %
            #                                   abs(flex.sum(rho.as_1d().select(gias)) * voxel_volume))

            peak_height = 0
            id = 0
            for n in range(len(cluster[x][0])):
                if real_map_unpadded.value_at_closest_grid_point(self.uc.fractionalize(cluster[x][0][n])) > peak_height:
                    peak_height = real_map_unpadded.value_at_closest_grid_point(self.uc.fractionalize(cluster[x][0][n]))
                    id = n
            site_cart.append(list(cluster[x][0][id]))


            fraction = self.uc.fractionalize(negative_cluster)
            results_diff = make_parallel(result_cluster)(fraction)
            results_diff_e_sum = make_parallel(result_cluster_e_sum)(fraction)
            peak_sum.append(sum(results_diff))
            electron_sums_around_atoms.append("%.2f" %
                                              (np.sum(np.absolute(np.array(results_diff_e_sum))) * voxel_volume))

        return {"cluster_number": cluster_number,
                "blob_volume": blob_volume,
                "site_cart": np.around(site_cart, 3).tolist(),
                "electron_sum_around_the_site": electron_sums_around_atoms, "number of voxel": cluster_size,"peak_sum":peak_sum}

    def write_cluster(self):
        df_p = pd.DataFrame(self.cluster_p_peak_volume)
        df_p.to_csv("peak_cluster_p.txt", index=False)
        df_n = pd.DataFrame(self.cluster_p_peak_volume)
        df_n.to_csv("peak_cluster_n.txt", index=False)

    def df_to_json(self, data_frame, output="output",write=True):
        df_all_group = data_frame.groupby("resn")
        group_key_list = list(df_all_group.groups.keys())
        df3 = pd.DataFrame()
        for x in group_key_list:
            df3 = df3.append(df_all_group.get_group(x).groupby('Period', group_keys=False).apply(
                lambda df: df.loc[df[["number of voxel"]].astype(float).idxmax()]))

        dict_test = {}
        for x in group_key_list:
            dict_test[str(x)] = df3.groupby("resn").get_group(x).apply(list).to_dict()

        d = defaultdict(list)

        for k in list(dict_test.keys()):
            for key in dict_test[k].keys():
                d[key].append(dict_test[k][key])
        if write:
            with open(f"{output}.json", "w") as outfile:
                json.dump(d, outfile)

        return df3

    def cluster_combine(self, df_p, df_n):
        common = list(set(df_n["resn"].values.tolist()) & set(df_p["resn"].values.tolist()))
        common.sort()
        cluster = pd.concat([df_n, df_p])
        return cluster[cluster["resn"].isin(common)], cluster









class Peaks(CctbxFunc, Peaks_coff):
    def __init__(self, pdb_file, map_file, resolution, sigma, name="default"):
        super().__init__(pdb_file, map_file, resolution, sigma)
        os.makedirs("Data_folder", exist_ok=True)
        os.chdir("Data_folder")
        self.name = f"map_dump_%s" % name
        os.makedirs(self.name, exist_ok=True)
        self.map_model = map_model_manager(map_manager=self.map, model=self.model)
        self.use_selection = self.map_model.model().get_xray_structure().all_selection()
        self.map_p = self.map_model.deep_copy()
        self.map_n = self.map_model.deep_copy()
        self.real_map_unpadded_p = self.map_p.map_as_fourier_coefficients(
            d_min=resolution).fft_map().apply_sigma_scaling().real_map_unpadded()
        self.real_map_unpadded_n = self.map_n.map_as_fourier_coefficients(
            d_min=resolution).fft_map().apply_sigma_scaling().real_map_unpadded() * -1
        self.uc = self.map_model.crystal_symmetry().unit_cell()
        self.grid_p = self.grid(self.real_map_unpadded_p, self.sigma)
        self.grid_n = self.grid(self.real_map_unpadded_n, self.sigma)
        self.selection = self.map_model.model().get_xray_structure().all_selection()
        print("chunk_1 ")
        self.site_frac_p = self.site_fraction(self.grid_p)
        self.site_frac_n = self.site_fraction(self.grid_n)
        self.site_cart_p = self.uc.orthogonalize(self.peak_near_model(self.real_map_unpadded_p, self.site_frac_p))
        self.site_cart_n = self.uc.orthogonalize(self.peak_near_model(self.real_map_unpadded_n, self.site_frac_n))
        self.cluster_p = self.DBscan(self.site_cart_p.as_numpy_array())
        self.cluster_n = self.DBscan(self.site_cart_n.as_numpy_array())
        print("chunk_2 ")
        self.cluster_p_peak_volume = self.peak_volume(self.real_map_unpadded_p, self.cluster_p, sigma)
        self.cluster_n_peak_volume = self.peak_volume(self.real_map_unpadded_n, self.cluster_n, sigma)
        self.peak_height_p = self.peak_height(self.real_map_unpadded_p, self.cluster_p_peak_volume)
        self.peak_height_n = self.peak_height(self.real_map_unpadded_n, self.cluster_n_peak_volume)
        self.unit_call_volume = self.map_model.model().get_xray_structure().unit_cell().volume()
        self.cluster_p_peak_volume_df = pd.DataFrame(self.cluster_p_peak_volume)
        self.cluster_n_peak_volume_df = pd.DataFrame(self.cluster_n_peak_volume)
        print("chunk_3 ")
        self.cluster_n_peak_volume_df["peak_height"] = self.cluster_n_peak_volume_df["peak_height"] * -1
        self.cluster_p_site = self.peak_site_frac(self.cluster_p_peak_volume_df)
        self.cluster_n_site = self.peak_site_frac(self.cluster_n_peak_volume_df)
        self.cluster_positive_p = self.peak_list_cluster(self.cluster_p_site, self.cluster_p_peak_volume_df)
        self.cluster_positive_n = self.peak_list_cluster(self.cluster_n_site, self.cluster_n_peak_volume_df)
        self.cluster_p_list = self.cluster_list(self.cluster_positive_p)
        self.cluster_n_list = self.cluster_list(self.cluster_positive_n)
        print("chunk_4 ")
        self.cluster_density_p = self.cluster_density_add(self.cluster_p_list, self.cluster_p_peak_volume_df,
                                                          self.cluster_p)
        self.cluster_density_n = self.cluster_density_add(self.cluster_n_list, self.cluster_n_peak_volume_df,
                                                          self.cluster_n)
        self.cluster_positive_negative = self.cluster_combine(self.cluster_density_p, self.cluster_density_n)
        print("chunk_5 ")
        os.chdir('../')
        self.cluster_density_n.to_csv(f"%s/%s/map_dump_negative_full.csv" % ("Data_folder", self.name))
        self.cluster_density_p.to_csv(f"%s/%s/map_dump_positive_full.csv" % ("Data_folder", self.name))


def run(pdb, Map, resolution, sigma):
    sc = Peaks(pdb_file=pdb, map_file=Map, resolution=resolution, sigma=sigma, name="default")


def run_2_all():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="pdb_file")
    parser.add_argument("map", help="Fo-Fo_map_file in CCP4 format")
    parser.add_argument("resolution", type=float, help="resolution")
    parser.add_argument("sigma", nargs='?', action='store', type=float, help="sigma")
    args = parser.parse_args()
    sigma = 3.0 if args.sigma is None else args.sigma

    run(args.pdb, args.map, args.resolution, sigma)


if __name__ == "__main__":
    run_2_all()
