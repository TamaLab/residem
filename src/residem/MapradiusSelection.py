from mmtbx.find_peaks import peaks_holder
from cctbx import maptbx
from scitbx.array_family import flex
from residem.scan_all_points import Peaks_coff, Peaks
import numpy as np
from residem.point_inside_polygon import PolyPoint
import warnings
import pandas as pd
from pandas.core.common import SettingWithCopyWarning
from residem.multiprocess_func import make_parallel
import tqdm
import os
from residem.cctbx_func import CctbxFunc

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)


class MapPeakSelection_coff(Peaks_coff):
    def __init__(self, map_model, resolution, sigma, max_dist=2.0, min_dist=0, name="default", write=False, path=None,
                 positive=False):
        super().__init__(map_model, resolution, sigma, path)
        self.path = path if not None else os.getcwd()
        self.name = f"map_dump_%s" % name
        if not os.path.exists(self.path + f"/{self.name}"):
            os.mkdir(self.path + f"/{self.name}")
        self.map_model = map_model
        self.cgt = self.crystal_gridding()
        self.max_dist = max_dist if not None else 2.0
        self.min_dist = min_dist if not None else 0.0
        self.use_selection = self.map_model.model().get_xray_structure().all_selection()
        self.uc = self.map_model.crystal_symmetry().unit_cell()
        self.selection = self.map_model.model().get_xray_structure().all_selection()
        self.sigma = sigma
        self.map_coefficient = self.map_model.map_as_fourier_coefficients(d_min=resolution)
        self.real_map = self.map_coefficient.fft_map()
        self.real_map_unpadded = self.real_map.apply_sigma_scaling().real_map_unpadded()
        self.positive_list = self.peak_list(self.cgt[0], self.cgt[1], self.cgt[2], self.peak_search_parameters(),
                                            self.sigma, "positive")
        self.negative_list = self.peak_list(self.cgt[0], self.cgt[1] * -1, self.cgt[2], self.peak_search_parameters(),
                                            self.sigma, "negative")
        self.df_n = self.data_frame(self.negative_list)
        self.df_p = self.data_frame(self.positive_list)
        self.df_n["peak height"] = self.df_n["peak height"] * -1
        self.grouped_df = self.group_df()
        self.cluster_all = self.cluster_identfier()
        self.post_negative_common = self.postitive_negative()

        if positive:
            self.positive_max = self.cluster_density_p.loc[
                self.cluster_density_p["number of voxel"] == self.cluster_density_p["number of voxel"].max()]
            self.cluster_num = np.array(self.positive_max["Cluster number"], dtype=np.int)[0]
            self.residue = np.array(self.positive_max["resn"], dtype=np.int)[0]
        else:
            self.max_col = self.cluster_all.loc[
                self.cluster_all["number of voxel"] == self.cluster_all["number of voxel"].max()]

            self.cluster_num = np.array(self.max_col["Cluster number"], dtype=np.int)[0]
            self.residue = np.array(self.max_col["resn"], dtype=np.int)[0]

        if write:
            # Here the file has to be split based on chain for better visualization.

            self.save_json_csv(self.cluster_positive_negative[0], self.path, self.name,
                               "map_dump_common_both_postv_n_negtv")
            self.save_json_csv(self.cluster_positive_negative[1], self.path, self.name,
                               "map_dump_full_positive_negative")
            self.save_json_csv(self.cluster_density_n, self.path, self.name, "map_dump_full_negative")
            self.save_json_csv(self.cluster_density_p, self.path, self.name, "map_dump_full_positive")

            self.df_json_full = self.df_to_json(self.cluster_positive_negative[0], write=False)
            self.df_json_full_1 = self.df_to_json(self.cluster_positive_negative[1], write=False)
            self.df_json_full_1["resn"] = self.df_json_full_1["resn"].astype(int)
            self.df_json_full_1["peak height"] = self.df_json_full_1["peak height"].astype(float)
            self.df_json_full_1["number of voxel"] = self.df_json_full_1["number of voxel"].astype(int)
            # self.cluster_positive_negative[0].to_csv(
            #     f"%s/%s/map_dump_positive_negative_both.csv" % (self.path, self.name))
            #
            # self.cluster_density_n.to_csv(f"%s/%s/map_dump_negative_full.csv" % (self.path, self.name))
            # self.cluster_density_p.to_csv(f"%s/%s/map_dump_positive_full.csv" % (self.path, self.name))
            # self.cluster_positive_negative[1].to_csv(f"%s/%s/map_dump_full.csv" % (self.path, self.name))
            # self.df_json_negative = self.df_to_json(self.cluster_density_n,
            #                                         f"%s/%s/map_dump_negative" % (self.path, self.name))
            # self.df_json_positive = self.df_to_json(self.cluster_density_p,
            #                                         f"%s/%s/map_dump_positive" % (self.path, self.name))

            self.max_col = self.cluster_all.loc[
                self.cluster_all["number of voxel"] == self.cluster_all["number of voxel"].max()]
        #
        # # print(list(self.max_col["Cluster number"])[0],list(self.max_col["resn"])[0])
        #
        #  self.cluster_num = list(self.max_col["Cluster number"])[0]
        #  self.residue = list(self.max_col["resn"])[0]

        # pd.DataFrame(self.cluster_n).to_csv(f"%s/%s/Cluster_n.csv" % (self.path, self.name))
        # pd.DataFrame(self.cluster_n).loc[self.cluster_num].to_csv(
        #     f"%s/%s/Cluster_n_max.csv" % (self.path, self.name), index=False)

    def save_json_csv(self, df_frame, path, name, save_as):
        comman_map = df_frame
        comman_map[['chain', 'resi']] = comman_map["resi"].str.split('-', expand=True)
        unique_chains = list(comman_map['chain'].unique())
        for chain in unique_chains:
            chain_case = f"{chain}_U" if chain.isupper() else f"{chain}_L"
            data_frame = comman_map[comman_map['chain'] == chain]
            data_frame['resi'] = data_frame['chain'].str.cat(data_frame['resi'], sep='-')
            data_frame['resn'] = data_frame['resn'].astype(int)
            data_frame["peak height"] = data_frame["peak height"].astype(float)
            data_frame["number of voxel"] = data_frame["number of voxel"].astype(int)
            os.makedirs(f"%s/%s/chain_%s_csv" % (path, name, chain_case), exist_ok=True)
            os.makedirs(f"%s/%s/chain_%s_json" % (path, name, chain_case), exist_ok=True)
            data_frame.to_csv(f"%s/%s/chain_%s_csv/%s_chain_%s.csv" % (path, name, chain_case, save_as, chain_case))
            self.df_to_json(data_frame,
                            f"%s/%s/chain_%s_json/%s_chain_%s" % (path, name, chain_case, save_as, chain_case))

    def peak_search_parameters(self):
        """ Peak search parameter with input cutoff"""
        peak_search_parameters = maptbx.peak_search_parameters(peak_search_level=1, max_peaks=0,
                                                               peak_cutoff=self.sigma, interpolate=True,
                                                               min_distance_sym_equiv=None,
                                                               general_positions_only=False,
                                                               min_cross_distance=1.5, min_cubicle_edge=5)
        return peak_search_parameters

    def peak_list(self, map_object, real_map, crystal_grid_tag, peak_search_parameter, cutoff,
                  value):
        """Getting the peak list as list of cartesian coordinates"""
        psr = crystal_grid_tag.peak_search(parameters=peak_search_parameter,
                                           map=real_map)
        connect = maptbx.connectivity(map_data=real_map, threshold=cutoff, wrapping=False)
        uc = map_object.crystal_symmetry().unit_cell()
        use_selection = map_object.model().get_xray_structure().all_selection()

        site_grid = connect.maximum_coors()
        vol_site = flex.vec3_double()

        # print("grid to fractional co-ordinate")
        #
        # def grid_2_frac(site_grid):
        #     result = [m / n for m, n in zip((site_grid[0], site_grid[1], site_grid[2]), real_map.all())]
        #     return result
        #
        # result = make_parallel(grid_2_frac)(tqdm.tqdm(site_grid.size()))
        # vol_site = flex.vec3_double(result)
        # print(result)

        for i in range(site_grid.size()):
            vol_site.append(
                [m / n for m, n in zip((site_grid[i][0], site_grid[i][1], site_grid[i][2]), real_map.all())])

        result_near = map_object.model().get_xray_structure().closest_distances(sites_frac=vol_site,
                                                                                distance_cutoff=2,
                                                                                use_selection=use_selection)
        result_with_in = map_object.model().get_xray_structure().closest_distances(sites_frac=psr.sites(),
                                                                                   distance_cutoff=1,
                                                                                   use_selection=use_selection)

        iseqs_of_closest_atoms_near, selection_near = self.closest_near_peak(result_near, self.max_dist, self.min_dist)
        iseqs_of_closest_atoms_within, selection_within = self.closest_near_peak(result_with_in, self.max_dist,
                                                                                 self.min_dist)

        peaks_near = peaks_holder(heights=connect.maximum_values().select(selection_near),
                                  sites=result_near.sites_frac.select(selection_near),
                                  iseqs_of_closest_atoms=iseqs_of_closest_atoms_near)
        peaks_with_in = peaks_holder(heights=psr.heights().select(selection_within),
                                     sites=result_with_in.sites_frac.select(selection_within),
                                     iseqs_of_closest_atoms=iseqs_of_closest_atoms_within)

        get_lis = list(connect.maximum_values())
        file_dict = {"near": peaks_near, "within": peaks_with_in}
        # file_dict = {"within": peaks_with_in}

        peak_list = []
        for key, file in file_dict.items():
            for length in range(file.sites.size()):
                for m in map_object.model().get_hierarchy().models():
                    for chain in m.chains():
                        for residue_group in chain.residue_groups():
                            for conformers in residue_group.conformers():
                                for residue in conformers.residues():
                                    for atoms in residue.atoms():
                                        if file.iseqs_of_closest_atoms[length] == atoms.i_seq:
                                            peak = uc.orthogonalize(file.sites[length])
                                            if value == "positive":
                                                if key == "near":
                                                    peak_list.append("%s,%s,%s,%.4f,%.4f,%.4f,%.3f,%.3f,%s" % (
                                                        residue.resseq, str(chain.id + '-' + residue.resname),
                                                        atoms.name,
                                                        peak[0], peak[1], peak[2],
                                                        file.heights[length],
                                                        connect.regions()[get_lis.index(file.heights[length])], "near"))
                                                elif key == "within":
                                                    peak_list.append("%s,%s,%s,%.4f,%.4f,%.4f,%.3f,%.3f,%s" % (
                                                        residue.resseq, str(chain.id + '-' + residue.resname),
                                                        atoms.name,
                                                        peak[0], peak[1], peak[2],
                                                        file.heights[length],
                                                        1, "within"))

                                            elif value == 'negative':
                                                if key == "near":
                                                    peak_list.append("%s,%s,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%s" % (
                                                        residue.resseq, str(chain.id + '-' + residue.resname),
                                                        atoms.name,
                                                        peak[0], peak[1], peak[2],
                                                        file.heights[length],
                                                        connect.regions()[get_lis.index(file.heights[length])], "near"))
                                                elif key == "within":
                                                    peak_list.append("%s,%s,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%s" % (
                                                        residue.resseq, str(chain.id + '-' + residue.resname),
                                                        atoms.name,
                                                        peak[0], peak[1], peak[2],
                                                        file.heights[length],
                                                        1, "within"))

        return peak_list

    @staticmethod
    def data_frame(peak_list):
        """List to data frame """
        df = pd.DataFrame(peak_list)
        df = df[0].str.split(",", expand=True)
        df.columns = ["resn", "resi", "atom name", "X", "Y", "Z", "peak height", "number of voxel", "site"]
        df[["peak height", "number of voxel"]] = df[["peak height", "number of voxel"]].astype(float)
        df['Period'] = df.resn.str.cat(df[["resi", "atom name"]], sep='-').str.replace(' ', '')
        df = df.groupby('Period', group_keys=False).apply(lambda df: df.loc[df[["number of voxel"]].idxmax()])
        df = df.sort_values(["number of voxel"], ascending=False)
        df = df.reset_index(drop=True)
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        return df

    def group_df(self):
        frame = [self.df_n, self.df_p]
        df_all = pd.concat(frame)
        df_all_group = df_all.groupby("resn")
        group_key_list = list(df_all_group.groups.keys())
        df3 = pd.DataFrame()

        for x in group_key_list:
            df3 = df3.append(df_all_group.get_group(x).groupby('Period', group_keys=False).apply(
                lambda df: df.loc[df[["number of voxel"]].idxmax()]))

        cols = ["resn", "X", "Y", "Z", "peak height", "number of voxel"]
        df3[cols] = df3[cols].apply(lambda x: pd.to_numeric(x, errors='coerce'))
        return df3

    def cluster_identfier(self):
        data_cluster = pd.DataFrame(
            columns=["resn", "resi", "atom name", "X", "Y", "Z", "peak height", "number of voxel", "site",
                     "Cluster number",
                     "Peak Volume of the Cluster", "electron_sum_around_the_site", "Period", "peak_sum"])
        df_p = self.grouped_df[self.grouped_df["peak height"] >= 0]
        df_n = self.grouped_df[self.grouped_df["peak height"] <= 0]
        dict_cluster = {"positive": [df_p, self.cluster_p], "negative": [df_n, self.cluster_n]}
        for key, value in dict_cluster.items():
            for r in range(len(value[0])):
                row_1 = value[0].iloc[r]

                row = value[0].iloc[:, 3:6].values
                x, y, z = row[r][0], row[r][1], row[r][2]
                res_list = value[1]
                for cluster_number in range(len(res_list)):
                    if PolyPoint(np.round([x, y, z], 3).tolist(),
                                 np.round(res_list[cluster_number][0], 3).tolist()).inside_polygon:
                        if key == "positive":
                            df_cluster = self.cluster_p_peak_volume_df[
                                self.cluster_p_peak_volume_df["cluster_number"] == cluster_number]
                            data_cluster = data_cluster.append({"resn": row_1["resn"],
                                                                "resi": row_1["resi"],
                                                                "atom name": row_1["atom name"],
                                                                "X": x, "Y": y, "Z": z,
                                                                "peak height": df_cluster["peak_height"].values[0],
                                                                "number of voxel": row_1["number of voxel"],
                                                                "site": row_1["site"],
                                                                "Cluster number": cluster_number,
                                                                "Peak Volume of the Cluster":
                                                                    df_cluster["blob_volume"].values[
                                                                        0],
                                                                "electron_sum_around_the_site"
                                                                : df_cluster["electron_sum_around_the_site"].values[0],
                                                                "Period": row_1["Period"],
                                                                "peak_sum": df_cluster["peak_sum"].values[0]},
                                                               ignore_index=True)
                        elif key == "negative":
                            df_cluster = self.cluster_n_peak_volume_df[
                                self.cluster_n_peak_volume_df["cluster_number"] == cluster_number]
                            data_cluster = data_cluster.append({"resn": row_1["resn"],
                                                                "resi": row_1["resi"],
                                                                "atom name": row_1["atom name"],
                                                                "X": x, "Y": y, "Z": z,
                                                                "peak height": df_cluster["peak_height"].values[0],
                                                                # * -1,
                                                                "number of voxel": row_1["number of voxel"],
                                                                "site": row_1["site"],
                                                                "Cluster number": cluster_number,
                                                                "Peak Volume of the Cluster":
                                                                    df_cluster["blob_volume"].values[
                                                                        0],
                                                                "electron_sum_around_the_site"
                                                                : df_cluster["electron_sum_around_the_site"].values[
                                                                    0],
                                                                "Period": row_1["Period"],
                                                                "peak_sum": df_cluster["peak_sum"].values[0]},
                                                               ignore_index=True)
        return data_cluster

    def postitive_negative(self):
        # df_n = self.cluster_all[self.cluster_all["peak height"] <= 0]
        # df_n["electron_sum_around_the_site"] = df_n["electron_sum_around_the_site"].apply(pd.to_numeric)
        # df_n["electron_sum_around_the_site"] = df_n["electron_sum_around_the_site"].apply(lambda x: x * -1)
        # df_p = self.cluster_all[self.cluster_all["peak height"] >= 0]
        df_n = self.cluster_density_n
        df_p = self.cluster_density_p
        common = list(set(df_n["resn"].values.tolist()) & set(df_p["resn"].values.tolist()))
        common.sort()
        cluster = pd.concat([df_n, df_p])
        cluster = cluster.loc[cluster["Peak Volume of the Cluster"] > 1]
        cluster = cluster.loc[cluster["number of voxel"].gt(1) & cluster["number of voxel"].le(2000)]
        return cluster[cluster["resn"].isin(common)], cluster


class MapPeakSelection(Peaks, MapPeakSelection_coff):
    def __init__(self, pdb_file, map_file, resolution, sigma, max_dist=2.0, min_dist=0, name="default", path=None):
        super().__init__(pdb_file, map_file, resolution, sigma)
        self.path = path if not None else os.getcwd()
        self.name = f"map_dump_%s" % name
        if not os.path.exists(self.path + f"/{self.name}"):
            os.mkdir(self.path + f"/{self.name}")
        self.cgt = self.crystal_gridding()
        self.max_dist = max_dist if not None else 2.0
        self.min_dist = min_dist if not None else 0.0
        self.use_selection = self.map_model.model().get_xray_structure().all_selection()
        self.uc = self.map_model.crystal_symmetry().unit_cell()
        self.selection = self.map_model.model().get_xray_structure().all_selection()
        self.sigma = sigma
        self.map_coefficient = self.map_model.map_as_fourier_coefficients(d_min=resolution)
        self.real_map = self.map_coefficient.fft_map()
        self.real_map_unpadded = self.real_map.apply_sigma_scaling().real_map_unpadded()
        self.positive_list = self.peak_list(self.cgt[0], self.cgt[1], self.cgt[2], self.peak_search_parameters(),
                                            self.sigma, "positive")
        self.negative_list = self.peak_list(self.cgt[0], self.cgt[1] * -1, self.cgt[2], self.peak_search_parameters(),
                                            self.sigma, "negative")
        self.df_n = self.data_frame(self.negative_list)
        self.df_p = self.data_frame(self.positive_list)
        self.df_n["peak height"] = self.df_n["peak height"] * -1
        self.grouped_df = self.group_df()
        self.cluster_all = self.cluster_identfier()
        self.post_negative_common = self.postitive_negative()

        # self.df_json_positive_negative = self.df_to_json(self.post_negative_common[0],
        #                                                  f"%s/%s/map_dump" % (self.path, self.name))
        self.df_json_full = self.df_to_json(self.cluster_positive_negative[0],
                                            f"%s/%s/map_dump_full" % (self.path, self.name))
        self.cluster_positive_negative[0].to_csv(f"%s/%s/map_dump_positive_negative_both.csv" % (self.path, self.name))
        # self.post_negative_common[1].to_csv(f"%s/%s/map_dump_positive_negative_full.csv" % (self.path, self.name))
        # self.post_negative_common[1]['electron_sum_around_the_site']= self.post_negative_common[1]['electron_sum_around_the_site'].astype(float)
        # self.dfn_1 = self.post_negative_common[1][self.post_negative_common[1]['electron_sum_around_the_site'] <=0 ]
        # self.dfp_1 = self.post_negative_common[1][self.post_negative_common[1]['electron_sum_around_the_site'] >=0 ]
        self.cluster_density_n.to_csv(f"%s/%s/map_dump_negative_full.csv" % (self.path, self.name))
        self.cluster_density_p.to_csv(f"%s/%s/map_dump_positive_full.csv" % (self.path, self.name))

        # self.df_json_positive_negative = self.df_to_json(self.post_negative_common[0],
        #                                                  f"%s/%s/map_dump" % (self.path, self.name))
        # self.df_json_full = self.df_to_json(self.post_negative_common[1],
        #                                     f"%s/%s/map_dump_full" % (self.path, self.name))
        # self.post_negative_common[0].to_csv(f"%s/%s/map_dump_positive_negative_common.csv" % (self.path, self.name))
        # self.post_negative_common[1].to_csv(f"%s/%s/map_dump_positive_negative_full.csv" % (self.path, self.name))
        # self.post_negative_common[1]['electron_sum_around_the_site']= self.post_negative_common[1]['electron_sum_around_the_site'].astype(float)
        # self.dfn_1 = self.post_negative_common[1][self.post_negative_common[1]['electron_sum_around_the_site'] <=0 ]
        # self.dfp_1 = self.post_negative_common[1][self.post_negative_common[1]['electron_sum_around_the_site'] >=0 ]
        # self.dfn_1.to_csv(f"%s/%s/map_dump_negative_full.csv" % (self.path, self.name))
        # self.dfp_1.to_csv(f"%s/%s/map_dump_positive_full.csv" % (self.path, self.name))

        self.max_col = self.cluster_all.loc[
            self.cluster_all["number of voxel"] == self.cluster_all["number of voxel"].max()]

        self.cluster_num = list(self.max_col["Cluster number"])[0]
        self.residue = list(self.max_col["resn"])[0]
        #
        # pd.DataFrame(self.cluster_n).to_csv(f"%s/%s/Cluster_n.csv" % (self.path, self.name))
        # pd.DataFrame(self.cluster_n).loc[self.cluster_num].to_csv(f"%s/%s/Cluster_n_max.csv" % (self.path, self.name),
        #                                                           index=False)


def run(pdb_file, map_file, d_min, cutoff, max_dist=2.0, min_dist=0.0):
    test = MapPeakSelection(pdb_file, map_file, d_min, cutoff, max_dist, min_dist)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="pdb_file")
    parser.add_argument("map", nargs='?', help="Fo-Fo_map_file in CCP4 format")
    parser.add_argument("d_min", nargs='?', type=float, help="resolution cutoff")
    parser.add_argument("s", nargs='?', type=float, help="sigama cutoff")
    parser.add_argument("max", nargs='?', type=float, help="maximum distance cutoff")
    parser.add_argument("min", nargs='?', type=float, help="minimum cutoff")
    args = parser.parse_args()
    run(args.pdb, args.map, args.d_min, args.s, args.max, args.min)
