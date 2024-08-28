import pandas as pd
from scitbx.array_family import flex
from cctbx import maptbx, miller
import mmtbx.utils
# import cctbx_func
from residem.cctbx_func import CctbxFunc
from iotbx.map_model_manager import map_model_manager
from residem.structUtils import ModelSelection
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import mmtbx
import mmtbx.model
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager
import iotbx.pdb

matplotlib.use('Agg')
font = {'weight': 'bold','family':'sans-serif','sans-serif':['Helvetica']}

matplotlib.rc('font', **font)

import warnings
import os

warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)


class AtomPeaks:
    def __init__(self, pdb, map_model, resolution, sigma, path=os.getcwd(), atom_profile="mean", residue=None,
                 name="default", weight="default"):

        # os.makedirs("Data_folder", exist_ok=True)
        # os.chdir("Data_folder")
        self.name = f"map_dump_%s" % name
        # os.makedirs(self.name, exist_ok=True)
        self.pdb = pdb
        # self.pdb_in = iotbx.pdb.input(pdb_file)
        # self.map = map_manager(map_file)
        # self.model = mmtbx.model.manager(model_input=self.pdb_in)
        # self.map_model = map_model_manager(map_manager=self.map, model=self.model)
        # self.map_model.box_all_maps_around_density_and_shift_origin()
        self.map_model = map_model
        self.sigma = sigma
        self.path = path
        os.chdir(self.path)
        self.file_path = self.path + "/" + self.name
        os.makedirs(self.file_path, exist_ok=True)
        # self.map_model = map_model_manager(map_manager=self.map, model=self.model)
        self.xrs = self.map_model.model().get_xray_structure()
        self.use_selection = self.map_model.model().get_xray_structure().all_selection()
        self.map_p = self.map_model.deep_copy()
        self.map_n = self.map_model.deep_copy()
        self.map_v = self.map_model.deep_copy()
        self.out_type = atom_profile if atom_profile == "max" else "mean"
        self.real_map_unpadded_p = self.map_p.map_as_fourier_coefficients(
            d_min=resolution).fft_map(resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
        self.real_map_unpadded_n = self.map_n.map_as_fourier_coefficients(
            d_min=resolution).fft_map().apply_sigma_scaling().real_map_unpadded() * -1
        self.uc = self.map_model.crystal_symmetry().unit_cell()
        self.atom_list = ModelSelection(self.pdb).atom_df
        self.atom_peak = self.peak_height()
        self.atom_peak_df = pd.DataFrame(self.atom_peak, columns=['chainid','resnum', 'resname', 'atomname', 'X', 'Y', 'Z',
                                                                  'Positive_max', 'Negative_max', 'Value_at_atom','Positive_mean','Negative_mean'])

        self.atom_peak_df['Period'] = self.atom_peak_df['chainid'].astype(str) + '-' + self.atom_peak_df['resnum'].astype(str) + '-' + self.atom_peak_df['resname'].astype(str) + '-' + self.atom_peak_df['atomname'].astype(str)
        self.atom_peak_df_sorted = self.atom_peak_df.sort_values(by=['Period', 'Positive_max', 'Negative_max'], ascending=[True, False, False])
        self.atom_peak_df = self.atom_peak_df.drop_duplicates(subset='Period')

        self.residue_to_plot = residue.split(',') if residue is not None else None

        comman_map = self.atom_peak_df
        unique_chains = list(comman_map['chainid'].unique())
        for chain in unique_chains:
            data_frame = comman_map[comman_map['chainid'] == chain]
            # print(chain)
            self.plot_residue(data_frame,chain)
            self.atom_plot(data_frame,chain)

        #self.atom_peak_df.to_csv(f"%s/chain_%s_U/Residual_peak_height.csv" % (self.file_path, "A"))
        # self.residue_plot = self.plot_residue(self.atom_peak_df)
        # self.atom_plot = self.atom_plot(self.atom_peak_df)
        os.chdir('../')

    def crystal_grid(self):
        cr_gr = maptbx.crystal_gridding(unit_cell=self.xrs.unit_cell(),
                                        space_group_info=self.xrs.space_group_info(),
                                        pre_determined_n_real=self.real_map_unpadded_p.accessor().all())
        return cr_gr

    def peak_height(self):
        atom_peak_height = []
        for i in self.atom_list.index:
            x, y, z = self.atom_list.loc[i].X, self.atom_list.loc[i].Y, self.atom_list.loc[i].Z
            gias = maptbx.grid_indices_around_sites(
                unit_cell=self.uc,
                fft_n_real=self.real_map_unpadded_p.all(),
                fft_m_real=self.real_map_unpadded_p.all(),
                sites_cart=flex.vec3_double([np.array([x, y, z])]),
                site_radii=flex.double([2]))

            positive_max = max(self.real_map_unpadded_p.select(gias))
            negative_max = min(self.real_map_unpadded_p.select(gias))
            positive_array = np.array((self.real_map_unpadded_p.select(gias)))
            positive_filtered = positive_array[positive_array > self.sigma]
            negative_filtered = positive_array[positive_array < -self.sigma]
            if len(positive_filtered) >= 1:
                positive_mean = positive_filtered.mean()
            else:
                positive_mean = positive_array[positive_array > 0 ].mean()
            if len(negative_filtered) >= 1:
                negative_filtered = positive_array[positive_array < -self.sigma]
                negative_mean = negative_filtered.mean()
            else:
                negative_mean = positive_array[positive_array < 0 ].mean()
            value_at_point = self.real_map_unpadded_p.value_at_closest_grid_point(self.uc.fractionalize([x, y, z]))
            atom_peak_height.append([self.atom_list.loc[i].chainid,self.atom_list.loc[i].resnum, self.atom_list.loc[i].resname,
                                     self.atom_list.loc[i].atomname, self.atom_list.loc[i].X,
                                     self.atom_list.loc[i].Y, self.atom_list.loc[i].Z, positive_max,
                                     negative_max, value_at_point,positive_mean,negative_mean])
        return atom_peak_height

    def plot_residue(self,df_to_process,chain):
        df_plot_res_list = []
        plot_res = df_to_process
        df_group = plot_res.groupby(['chainid', 'resnum'])
        # df_group = plot_res.groupby('resnum')
        gb = df_group.groups
        for key, values in gb.items():
            data_group = plot_res[plot_res['resnum'] == key[1]]
            data_group = data_group.drop_duplicates(subset='atomname')

            if self.out_type == "max":
                df_plot_res_list.append([data_group['Positive_max'].max(),
                                         data_group['Negative_max'].min(),
                                         data_group["Value_at_atom"].mean(),
                                         (str(data_group["chainid"].unique()[0])+ '-'+str(data_group["resnum"].unique()[0]) + '-' + str(
                                             data_group.resname.unique()[0]))])


            else:
                df_plot_res_list.append([data_group['Positive_mean'].mean(),
                                         data_group['Negative_mean'].mean(),
                                         data_group["Value_at_atom"].mean(),
                                         (str(data_group["chainid"].unique()[0])+ '-'+ str(data_group["resnum"].unique()[0]) + '-' + str(
                                             data_group.resname.unique()[0]))])

        df_plot_res = pd.DataFrame(df_plot_res_list,
                                   columns=["Positive_max", "Negative_max", "Value_at_atom", "Residue_name"])

        div = round(len(df_plot_res.index) / 50)
        if div > 1:
            # df_split  = [df_plot_res.iloc[i:i + div] for i in range(0, len(df_plot_res), div)]
            df_split = np.array_split(df_plot_res, div)
        else:
            df_split = [df_plot_res]

        chain_case = f"{chain}_U" if chain.isupper() else f"{chain}_L"

        os.makedirs(f"%s/chain_%s" % (self.file_path, chain_case), exist_ok=True)
        out_pdf = f'%s/chain_{chain_case}/Residual_peak_height_max_chain_{chain_case}.pdf' % self.file_path if self.out_type == "max" else f'%s/chain_{chain_case}/Residual_peak_height_mean_chain_{chain_case}.pdf' % self.file_path
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        df_plot_res.to_csv(f"%s/chain_%s/Residual_peak_height.csv"% (self.file_path, chain_case))
        figs = plt.figure()

        for idx, df in enumerate(df_split):
            fig, ax = plt.subplots(figsize=(10, 10))  # inches figsize=(10, 10)
            plt.subplot()

            if self.out_type == "max":
                plt.plot(df["Residue_name"], df['Positive_max'], 'o-', color='g',
                         label='Positive max of residue')
                plt.plot(df["Residue_name"], df['Negative_max'], 'o-', color='r',
                         label='Negative max of residue')
            else:
                plt.plot(df["Residue_name"], df['Positive_max'], 'o-', color='g',
                         label='Positive mean of residue')
                plt.plot(df["Residue_name"], df['Negative_max'], 'o-', color='r',
                         label='Negative mean of residue')

            plt.plot(df["Residue_name"], df['Value_at_atom'], 'o-', color='b',
                     label='Average of the values at residue')

            plt.xticks(rotation=90)
            plt.xlabel('Residue',weight='bold', size=16)
            plt.ylabel('Difference density',weight='bold', size=16)
            # plt.yticks(np.arange(min(round(df['Negative_max']) - 1), max(round(df['Positive_max'])) + 1, 1))
            plt.title('Difference density profile of residues',weight='bold', size=16)
            plt.axhline(self.sigma, color='green', lw=2, alpha=0.6)
            plt.fill_between(df["Residue_name"], 0, self.sigma, where=df['Positive_max'] > 0,
                             color='green', alpha=0.2)

            plt.fill_between(df["Residue_name"], 0, -self.sigma, where=df['Negative_max'] < 0,
                             color='red', alpha=0.2)

            plt.plot(df[(df['Positive_max'] > self.sigma)]["Residue_name"],
                     df['Positive_max'][(df['Positive_max'] > self.sigma)],
                     linestyle='none', marker='D', markersize=10, markeredgecolor="orange", markeredgewidth=5,
                     alpha=0.6)

            plt.plot(df[(df['Negative_max'] < -self.sigma)]["Residue_name"],
                     df['Negative_max'][(df['Negative_max'] < -self.sigma)],
                     linestyle='none', marker='P', markersize=10, markeredgecolor="orange", markeredgewidth=5,
                     alpha=0.6)

            empty_list = []
            empty_list.append(list(df[(df['Negative_max'] < -self.sigma)]["Residue_name"]))
            empty_list.append(list(df[(df['Positive_max'] > self.sigma)]["Residue_name"]))
            empty_list = [item for sublist in empty_list for item in sublist]
            ax.set_xticklabels(list(df["Residue_name"]), weight='bold', size=8)
            ax.set_yticklabels(ax.get_yticks(), weight='bold', size=8)


            for tl in ax.get_xticklabels():
                txt = tl.get_text()
                # print(txt)
                if txt in empty_list:
                    txt += ' (!)'
                    tl.set_backgroundcolor('bisque')
                    # tl.set_fontweight('bold')
                    # tl.set_fontsize(18)
                tl.set_text(txt)
            plt.grid()
            plt.legend()
            # plt.rcParams.update({'font.size': 8})
            pdf.savefig(fig)
        pdf.close()
        return True

    def atom_plot(self,df_to_process,chain):

        def split_atomname(atomname):
            head = ''.join(filter(str.isalpha, atomname))
            tail = ''.join(filter(str.isdigit, atomname))

            # Handling edge cases where tail might be empty
            if tail == '':
                return head, 0
            else:
                return head, int(tail)


        residue_list = []
        if self.residue_to_plot is not None:
            for i in self.residue_to_plot:
                try:
                    residue_list.append(int(i))
                except ValueError:
                    [residue_list.append(x) for x in
                     [x for x in range(int(i.split('-')[0]), (int(i.split('-')[1])) + 1)]]

        else:
            [residue_list.append(x) for x in df_to_process.resnum.unique()]

        chain_case = f"{chain}_U" if chain.isupper() else f"{chain}_L"
        os.makedirs(f"%s/chain_%s" % (self.file_path, chain_case), exist_ok=True)
        out_pdf = f'%s/chain_{chain_case}/Atom_peak_height_chain_%s.pdf' % (self.file_path,chain_case)
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        figs = plt.figure()
        plot_res = df_to_process[df_to_process['resnum'].isin(residue_list)]
        df_group = plot_res.groupby(['chainid','resnum'])
        gb = df_group.groups
        for key, values in gb.items():
            data_group = plot_res[plot_res['resnum'] == key[1]]
            data_group = data_group.drop_duplicates(subset='atomname')
            data_group[['text', 'number']] = data_group['atomname'].apply(lambda x: pd.Series(split_atomname(x)))
            data_group.sort_values(['text', 'number'], inplace=True)
            data_group.drop(['text', 'number'], axis=1, inplace=True)
            data_group.reset_index(drop=True, inplace=True)
            # data_group = data_group.sort_values(by=['atomname'])
            fig, ax = plt.subplots(figsize=(10, 10))
            plt.subplot()
            if self.out_type == "max":
                plt.plot(data_group["atomname"], data_group['Positive_max'], 'o-', color='g',
                         label='Positive max around atom')
                plt.plot(data_group["atomname"], data_group['Value_at_atom'], 'o-', color='b',
                         label='Value at atom')
                plt.plot(data_group["atomname"], data_group['Negative_max'], 'o-', color='r',
                         label='Negative max around atom')
                plt.axhline(self.sigma, color='green', lw=2, alpha=0.6)
                plt.fill_between(data_group["atomname"], 0, self.sigma, where=data_group['Positive_max'] > 0,
                                 color='green', alpha=0.2)

                plt.fill_between(data_group["atomname"], 0, -self.sigma, where=data_group['Negative_max'] < 0,
                                 color='red', alpha=0.2)

                plt.plot(data_group[(data_group['Positive_max'] > self.sigma)]["atomname"],
                         data_group['Positive_max'][(data_group['Positive_max'] > self.sigma)],
                         linestyle='none', marker='D', markersize=10, markeredgecolor="orange", markeredgewidth=5,
                         alpha=0.6)

                plt.plot(data_group[(data_group['Negative_max'] < -self.sigma)]["atomname"],
                         data_group['Negative_max'][(data_group['Negative_max'] < -self.sigma)],
                         linestyle='none', marker='P', markersize=10, markeredgecolor="orange", markeredgewidth=5,
                         alpha=0.6)
            elif self.out_type == "mean":
                plt.plot(data_group["atomname"], data_group['Positive_mean'], 'o-', color='g',
                         label='Positive mean around atom')
                plt.plot(data_group["atomname"], data_group['Value_at_atom'], 'o-', color='b',
                         label='Value at atom')
                plt.plot(data_group["atomname"], data_group['Negative_mean'], 'o-', color='r',
                         label='Negative mean around atom')
                plt.axhline(self.sigma, color='green', lw=2, alpha=0.6)
                plt.fill_between(data_group["atomname"], 0, self.sigma, where=data_group['Positive_mean'] > 0,
                                 color='green', alpha=0.2)

                plt.fill_between(data_group["atomname"], 0, -self.sigma, where=data_group['Negative_mean'] < 0,
                                 color='red', alpha=0.2)

                plt.plot(data_group[(data_group['Positive_mean'] > self.sigma)]["atomname"],
                         data_group['Positive_mean'][(data_group['Positive_mean'] > self.sigma)],
                         linestyle='none', marker='D', markersize=10, markeredgecolor="orange", markeredgewidth=5,
                         alpha=0.6)

                plt.plot(data_group[(data_group['Negative_mean'] < -self.sigma)]["atomname"],
                         data_group['Negative_mean'][(data_group['Negative_mean'] < -self.sigma)],
                         linestyle='none', marker='P', markersize=10, markeredgecolor="orange", markeredgewidth=5,
                         alpha=0.6)

            plt.xticks(rotation=90)
            plt.xlabel(f'Residue %s - %s' % (data_group["chainid"].unique()[0],data_group["resnum"].unique()[0]),weight='bold', size=16)
            plt.ylabel('Difference density',weight='bold', size=16)
            plt.title(f'Difference density profile of atom in chain %s residue number %s' % (data_group["chainid"].unique()[0],data_group["resnum"].unique()[0]),weight='bold', size=16)


            empty_list = []
            empty_list.append(list(data_group[(data_group['Negative_max'] < -self.sigma)]["atomname"]))
            empty_list.append(list(data_group[(data_group['Positive_max'] > self.sigma)]["atomname"]))
            empty_list = [item for sublist in empty_list for item in sublist]
            ax.set_xticklabels(list(data_group["atomname"]), weight='bold', size=8)
            ax.set_yticklabels(ax.get_yticks(), weight='bold', size=8)


            for tl in ax.get_xticklabels():
                txt = tl.get_text()
                # print(txt)
                if txt in empty_list:
                    txt += ' (!)'
                    tl.set_backgroundcolor('bisque')
                    # tl.set_fontweight('bold')
                    # tl.set_fontsize(18)
                tl.set_text(txt)

            plt.grid()
            plt.legend()
            # plt.rcParams.update({'font.size': 8})
            pdf.savefig(fig)
        pdf.close()
        return True


# def main(pdb_file, map_file, resolution, sigma, out_type="mean", residue=None):
#     test = AtomPeaks(pdb_file=pdb_file, map_file=map_file, resolution=resolution, sigma=sigma, out_type=out_type,
#                      residue=residue)
#
#
# def run_2_all():
#     import argparse
#
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-r", "--pdb", help="pdb_file")
#     parser.add_argument("-m", "--map", type=str, help="Fo-Fo_map_file in CCP4 format")
#     parser.add_argument("-d", "--d_min", type=float, help="high resolution cutoff")
#     parser.add_argument("-s", "--sigma", type=float, help="sigma cutoff")
#     parser.add_argument("-t", "--method", type=str, help="method (max or mean) of residual peak value")
#     parser.add_argument("-n", "--residue", type=str, help="residue selection")
#     args = parser.parse_args()
#     sigma = args.sigma if args.sigma is not None else 3.0
#     method = args.method if args.method in ["max", "mean"] else "max"
#     residue = args.residue if args.residue is not None else None
#     main(pdb_file=args.pdb, map_file=args.map, resolution=args.d_min, sigma=sigma, out_type=method, residue=residue)


if __name__ == "__main__":
    run_2_all()
