import pandas as pd
from collections import defaultdict
from glob import glob
import itertools
import numpy as np
import matplotlib.colors as mcolors
from matplotlib import pyplot as plt
import itertools
import tensorflow as tf
import json

plt.rcParams['image.cmap'] = 'RdBu_r'
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.colors import ListedColormap

matplotlib.use('Agg')
font = {'family':'sans-serif','sans-serif':['Helvetica']}

matplotlib.rc('font', **font)
import seaborn as sns
from iotbx import pdb
from scipy.spatial import distance
import networkx as nx
import matplotlib.colors as colors
from residem.structUtils import ModelSelection
from collections import defaultdict
import tensorflow as tf
import pandas as pd
import seaborn as sns



# mypalette = ['#4E4C4C', '#5A5757', '#726F6F', '#8A8686', '#9E9A9A', '#AEACAC', '#C4C1C1', '#D0CACA', '#E6E4E4']
# my_cmap = ListedColormap(sns.color_palette(mypalette).as_hex())
my_cmap = plt.cm.Blues


class Network(ModelSelection):
    """ A class to plot network of the given pdb and electron density sum, a residue can be provided as a list
    and a data frame with the electron sum has to be given if not it will be generated on its own"""

    def __init__(self, pdb_file, df, distance_type="CM", low_threshold=0.05, medium_threshold=0.1, map_type="peak", subgraph_residue=None,map_colour=None,closeness=0.005, norm=True,edge_width = 0.0015):
        super().__init__(pdb_file)
        self.df = df
        df_index = self.df.index
        df_array = self.df.values
        self.norm = norm
        self.edge_with_size = edge_width
        # min_max_scaler = preprocessing.MinMaxScaler()
        # df_normalized = min_max_scaler.fit_transform(df_array)
        # preprocessing.minmax_scale
        # df_normalized = tf.keras.utils.normalize(df_array, axis=0)
        if self.norm==True:
            df_normalized = (df_array - np.min(df_array)) / (np.max(df_array) - np.min(df_array))
            self.df_normalized = pd.DataFrame(df_normalized, columns=self.df.columns)
            self.df_normalized.index = df_index
            self.df = self.df_normalized
        self.map_colour = sns.light_palette(my_cmap , as_cmap=True) if map_colour is None else ListedColormap(sns.color_palette(map_colour).as_hex())#sns.light_palette(map_colour, as_cmap=True)
        self.protein_hierarchy = self.protein()
        self.subgraph_residue = subgraph_residue
        self.ligand = self.not_protein()
        self.distance_type = distance_type if distance_type in ["CA", "SC", "CM"] else "CM"
        self.low_threshold = low_threshold
        self.medium_medium_threshold = medium_threshold
        self.subplot_closeness = closeness
        self.residue_com_dict = self.residue_dict("protein")
        self.residue_com_dict.update(self.residue_dict("ligand"))
        self.residue_number = list(self.residue_com_dict.keys())
        self.distance = self.distance_cal()
        self.electron_sum_all = self.electron_sum()
        self.map_type = map_type if map_type in ["peak", "density"] else "peak"
        if self.map_type == "peak":
            self.title = "peak"
        elif self.map_type == "density":
            self.title = "electron density"

        self.type_name = "C-alpha" if self.distance_type == "CA" else "side chain" if self.distance_type == "SC" else "centre of mass " if self.distance_type == "CM" else None


    def subplot(self,label=True):
        out_pdf = f'Difference_density_network_%s_subplot.pdf' % self.distance_type
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        low_size = 300
        medium_size = 900
        large_size = 1200

        # Set the layout algorithm
        layout_algorithm = nx.spring_layout

        time_list = list(self.electron_sum_all.keys())

        node_color_max = []

        for i, x in enumerate(time_list):
            G = nx.Graph()
            for residue_number in self.subgraph_residue:
                G.add_node(residue_number,
                           electron_density_sum=self.electron_sum_all[x][self.residue_number.index(str(residue_number))])
            node_attributes = nx.get_node_attributes(G, 'electron_density_sum')
            node_color_max.append(np.array(list(node_attributes.values())).max())

        for i, x in enumerate(time_list):
            G = nx.Graph()
            for residue_number in self.subgraph_residue:
                G.add_node(residue_number,
                           electron_density_sum=self.electron_sum_all[x][self.residue_number.index(str(residue_number))])

            for residue_1_number in self.subgraph_residue:
                for residue_2_number in self.subgraph_residue:
                    if residue_1_number != residue_2_number:
                        dist = self.distance[
                            self.residue_number.index(str(residue_1_number)), self.residue_number.index(
                                str(residue_2_number))]

                        if int(residue_1_number) in list(self.df.index) and int(residue_2_number) in list(self.df.index) and dist<10:
                            G.add_edge(residue_1_number, residue_2_number, distance=dist *2)
                        elif dist<10 :
                            G.add_edge(residue_1_number, residue_2_number, distance=dist)


            if i == 0:
                positions = layout_algorithm(G,k=self.subplot_closeness, scale=0.0001,iterations=50)
            fig, ax = plt.subplots(figsize=(3.063*1.8, 2.108*1.8))
            node_attributes = nx.get_node_attributes(G, 'electron_density_sum')
            node_colors = list(node_attributes.values())

            # Determine node sizes based on electron density sum
            node_sizes = []
            for density_sum in node_colors:
                if density_sum < self.low_threshold:
                    node_sizes.append(low_size)
                elif density_sum < self.medium_medium_threshold and density_sum > self.low_threshold:
                    node_sizes.append(medium_size)
                else:
                    node_sizes.append(large_size)

            node_res = dict(zip(G.nodes(), node_sizes))

            # Get edge attributes (distance)
            edge_attributes = nx.get_edge_attributes(G, 'distance')
            edge_colors = list(edge_attributes.values())
            # edge_width = [self.edge_with_size * G[u][v]['distance'] for u, v in G.edges()]
            edge_width = [self.edge_with_size * x for x in self.min_max_scaling([G[u][v]['distance'] for u, v in G.edges()])]
            # edge_width = [self.edge_with_size * x for x in [G[u][v]['distance'] for u, v in G.edges()]]


            # Draw nodes and edges with fixed positions
            nx.draw_networkx_nodes(G, positions, node_size=node_sizes, node_color=node_colors, cmap=self.map_colour)
            nx.draw_networkx_edges(G, positions, edgelist=G.edges(), edge_color=edge_colors, edge_cmap=self.map_colour, width=edge_width,
                                   alpha=0.3)



            # Create a color legend for the electron density sum

            norm = colors.Normalize(vmin=0, vmax=max(np.array(node_color_max).flatten()))
            sm = plt.cm.ScalarMappable(cmap=self.map_colour, norm=norm)
            sm.set_array([])
            # plt.colorbar(sm, label=f'Normalized %s sum' % self.title)


            # Add labels to the nodes
            nx.draw_networkx_labels(G, positions,font_weight='bold',) if label == True else nx.draw_networkx_labels(G, positions, labels={k: k for k, v in node_res.items() if v > low_size})

            plt.title(f'Network graph for time  %s' % (x))
            # plt.title(f'Network graph based on %s distance and normalized %s sum for %s' % (self.type_name, self.title, x))
            plt.axis('off')
            plt.tight_layout()
            pdf.savefig(fig)

        pdf.close()

    @staticmethod
    def min_max_scaling(data):
        if not data:
            return []

        min_val = min(data)
        max_val = max(data)

        if min_val == max_val:
            return [0.5 for _ in data]  # All values are the same

        scaled_data = [(x - min_val) / (max_val - min_val) for x in data]
        return scaled_data


    def plot_network(self,label=False):
        out_pdf = f'Difference_density_network_%s.pdf' % self.distance_type if label == True else f'Difference_density_network_%s_nolabel.pdf' % self.distance_type
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)

        # Define node size ranges
        low_size = 100
        medium_size = 700
        large_size = 1000

        # Set the layout algorithm
        layout_algorithm = nx.spring_layout

        time_list = list(self.electron_sum_all.keys())

        node_color_max = []

        for i, x in enumerate(time_list):
            G = nx.Graph()
            for residue_number in self.residue_number:
                G.add_node(residue_number,
                           electron_density_sum=self.electron_sum_all[x][self.residue_number.index(residue_number)],width=0)
            node_attributes = nx.get_node_attributes(G, 'electron_density_sum')
            node_color_max.append(np.array(list(node_attributes.values())).max())

        #
        # print(node_color_max)
        # print(node_color_max[0])
        # print(max(np.array(node_color_max).flatten()))

        for i, x in enumerate(time_list):
            G = nx.Graph()
            for residue_number in self.residue_number:
                G.add_node(residue_number,
                           electron_density_sum=self.electron_sum_all[x][self.residue_number.index(residue_number)],width=0)



            for residue_1_number in self.residue_number:
                for residue_2_number in self.residue_number:
                    if residue_1_number != residue_2_number:
                        dist = self.distance[
                            self.residue_number.index(residue_1_number), self.residue_number.index(
                                residue_2_number)]

                        if int(residue_1_number) in list(self.df.index) and int(residue_2_number) in list(self.df.index) and dist < 10:

                            G.add_edge(residue_1_number, residue_2_number, distance=dist*2)
                        elif dist < 10 :
                            G.add_edge(residue_1_number, residue_2_number, distance=dist)



            if i == 0:
                positions = layout_algorithm(G)
            fig, ax = plt.subplots(figsize=(30, 15))

            node_attributes = nx.get_node_attributes(G, 'electron_density_sum')
            node_colors = list(node_attributes.values())

            # Determine node sizes based on electron density sum
            node_sizes = []
            for density_sum in node_colors:
                if density_sum < self.low_threshold:
                    node_sizes.append(low_size)
                elif density_sum < self.medium_medium_threshold and density_sum > self.low_threshold:
                    node_sizes.append(medium_size)
                else:
                    node_sizes.append(large_size)

            node_res = dict(zip(G.nodes(), node_sizes))

            # Get edge attributes (distance)
            edge_attributes = nx.get_edge_attributes(G, 'distance')
            edge_colors = list(edge_attributes.values())
            # edge_width = self.min_max_scaling([G[u][v]['distance'] for u, v in G.edges()])
            edge_width = [self.edge_with_size *x for x in [G[u][v]['distance'] for u, v in G.edges()]]

            # edge_width = [self.edge_with_size * x for x in self.min_max_scaling([G[u][v]['distance'] for u, v in G.edges()])]


            # Draw nodes and edges with fixed positions
            nx.draw_networkx_nodes(G, positions, node_size=node_sizes, node_color=node_colors, cmap=self.map_colour,edgecolors="tab:gray",alpha=0.9)
            nx.draw_networkx_edges(G, positions, edge_color=edge_colors, edge_cmap=self.map_colour, width=edge_width,
                                   alpha=0.3)


            # Create a color legend for the electron density sum

            norm = colors.Normalize(vmin=0, vmax=max(np.array(node_color_max).flatten()))
            sm = plt.cm.ScalarMappable(cmap=self.map_colour, norm=norm)
            sm.set_array([])
            if self.norm==True:
                cbar_1 = plt.colorbar(sm, label=f'Normalized %s sum' % self.title)
            else:
                cbar_1 = plt.colorbar(sm, label=f'%s sum' % self.title)

            cbar_1.set_label(f'Normalized {self.title} sum' if self.norm else f'{self.title} sum',
                           fontsize=28)

            cbar_1.ax.tick_params(labelsize=24)

            # Add labels to the nodes
            nx.draw_networkx_labels(G, positions, font_size=24) if label == True else nx.draw_networkx_labels(G, positions, labels={k: k for k, v in node_res.items() if v > low_size})




            if self.norm==True:
                plt.title(f'Network graph based on %s distance and normalized %s sum for %s' % (self.type_name, self.title, x),fontsize = 30)
            else:
                plt.title(f'Network graph based on %s distance and %s sum for %s' % (
                self.type_name, self.title, x),fontsize = 28)
            plt.axis('off')

            pdf.savefig(fig)

        pdf.close()


    def electron_sum(self):
        electron_sum = defaultdict(int)
        for x in list(self.df.columns):
            electron_sum[f"%s" % x] = []
            for y in list(self.residue_com_dict.keys()):
                if int(y) in list(self.df.index):
                    y = int(y)
                    electron_sum[f"%s" % x].append(float(self.df[x].loc[y]))
                else:
                    electron_sum[f"%s" % x].append(0)

        return electron_sum

    def residue_dict(self, cal_type="protein"):
        residue_com_dict = {}

        def residue_cal(residue):
            atom_selection = residue.atoms()
            if atom_selection.size() > 0:
                residue_coords = atom_selection.extract_xyz()
                residue_com = np.mean(residue_coords.as_numpy_array(), axis=0)
                residue_number = residue.resid().strip()
                residue_com_dict[residue_number] = residue_com

        if cal_type == "protein":
            if self.distance_type == "CA":
                pdb_hierarchy = self.calphas()
                for residue in pdb_hierarchy.residue_groups():
                    if residue.atoms().size() > 0:
                        residue_cal(residue)
            elif self.distance_type == "CM":
                pdb_hierarchy = self.protein()
                for residue in pdb_hierarchy.residue_groups():
                    if residue.atoms().size() > 0:
                        residue_cal(residue)

            elif self.distance_type == "SC":
                pdb_hierarchy = self.sidechains()
                for residue in pdb_hierarchy.residue_groups():
                    if residue.atoms().size() > 0:
                        residue_cal(residue)

        elif cal_type == "ligand":
            pdb_hierarchy = self.not_protein()
            for residue in pdb_hierarchy.residue_groups():
                if residue.atoms().size() > 0:
                    residue_cal(residue)

        return residue_com_dict

    def distance_cal(self):
        distances = np.zeros((len(self.residue_number), len(self.residue_number)))
        for i, residue_1 in enumerate(self.residue_number):
            for j, residue_2 in enumerate(self.residue_number):
                if residue_1 != residue_2:
                    com_1 = self.residue_com_dict[str(residue_1)]
                    com_2 = self.residue_com_dict[str(residue_2)]
                    dist = np.linalg.norm(com_1 - com_2)
                    distances[i, j] = dist

        return distances


# normalizing the array
# df_array = df.values
# df_normalized = tf.keras.utils.normalize(df_array, axis=0)
# df_normalized = pd.DataFrame(df_normalized, columns=df.columns)
class SVD:
    """ userID(Time) and title(atom id), rating ( peak volume)"""

    def __init__(self, list_of_dataframes, cmap=None, residue=None, cutoff_value=0, map_type="peak",norm=True, time = None):
        self.data_frame_list = list_of_dataframes
        self.cmap = sns.light_palette("lightsalmon", as_cmap=True) if cmap is not None else sns.light_palette(cmap, as_cmap=True)
        self.norm = norm
        self.cutoff = cutoff_value
        self.data_dict = defaultdict(int)
        self.data_group_dict = defaultdict(int)
        self.map_type = map_type if map_type in ["peak", "density"] else "peak"
        self.to_map = "electron_sum_around_the_site" if self.map_type == "density" else "peak_sum"
        self.data_dict_updated = self.data_frame_to_dict()
        self.data_frame_group()
        self.new_list = [x for x in self.data_group_dict.keys() if "merged" in x]
        self.residue = residue
        self.time = time


        self.df = self.data_frame_merge()

        # self.df = self.df.set_index("resn")
        self.df = self.df.sort_index(ascending=True)
        self.df.index = self.df.index.astype(int)
        self.residue_to_plot = self.part_residue() if residue is not None else None
        df_array = self.df.values
        self.nof = self.df.copy()
        # min_max_scaler = preprocessing.MinMaxScaler()
        if self.norm == True:
            df_normalized = (df_array - np.min(df_array)) / (np.max(df_array) - np.min(df_array))
            # df_normalized = tf.keras.utils.normalize(df_array, axis=0)
            self.df_normalized = pd.DataFrame(df_normalized, columns=self.df.columns)
            self.df_before = self.df
            self.df = self.df_normalized
            self.df.index = self.df_before.index


        else:
            self.df_before = self.df

        if self.residue_to_plot is not None:
            filtered_list = [item for item in self.residue_to_plot if item in self.df.index]
            self.df = self.df.loc[filtered_list] if self.residue_to_plot is not None else self.df
        self.filtered = self.df[self.df >= self.cutoff]
        self.filtered = self.filtered.dropna(how="all", axis=0)
        self.filtered = self.filtered.fillna(0)
        # self.df.drop(self.df.columns[self.df >= self.cutoff], axis=1, inplace=True)
        if self.time is None:
            self.filtered.columns = [f"Time_%s" % i for i in range(1, len(self.df.columns) + 1)]
            with open("time_and_file_name.txt", "w+") as file:
                for i, item in enumerate(self.data_frame_list):
                    file.write(f"Time_%s : %s \n" % (i + 1, item))
                file.close()
        else:
            self.filtered.columns = self.time
        self.X, self.U, self.S, self.V = self.visualize_svd(self.df.to_numpy())
        self.residue_list = list(self.filtered.index)

        self.nof.columns = self.filtered.columns



    def data_frame_to_dict(self):
        for i in range(len(self.data_frame_list)):
            self.data_dict[f"Data_frame_%s" % i] = pd.read_csv(self.data_frame_list[i])
            self.data_dict[f"Data_frame_%s" % i] = self.data_dict[f"Data_frame_%s" % i].loc[:,
                                                   ~self.data_dict[f"Data_frame_%s" % i].columns.str.contains(
                                                       '^Unnamed')]
        return self.data_dict

    def data_frame_group(self):
        for i in range(len(self.data_dict_updated)):
            self.data_group_dict[f"merged_%s" % i] = pd.DataFrame(columns=["resn", "electron_sum"])
            self.data_group_dict[f"merged_%s" % i].set_index("resn")
            self.data_group_dict[f"Data_frame_%s" % i] = self.data_dict[f"Data_frame_%s" % i].groupby("resn")
            for j in self.data_group_dict[f"Data_frame_%s" % i].groups.keys():
                self.data_group_dict[f"merged_%s" % i] = self.data_group_dict[f"merged_%s" % i].append(
                    {"resn": int(j), "electron_sum":
                        self.data_group_dict[f"Data_frame_%s" % i].get_group(j)[
                            self.to_map].abs().sum()},
                    ignore_index=True)

    def data_frame_merge(self):
        new_dict = defaultdict(int)
        for x in self.new_list:
            self.data_group_dict[x] = self.data_group_dict[x].set_index('resn')
            new_dict[x] = self.data_group_dict[x]
        result = pd.concat([new_dict[x] for x in self.new_list], axis=1, join='outer')
        result = result.fillna(0)
        return result

    #     for i in range(1, len(self.data_dict_updated)):
    #         self.data_group_dict[f"merged_{0}"] = self.data_group_dict[f"merged_{0}"].merge(
    #             self.data_group_dict[f"merged_{i}"],
    #             on='resn', how='left', suffixes=(i - 1, i + 0))
    #     return self.data_group_dict[f"merged_{0}"].fillna(0)

    @staticmethod
    def svd(X):
        U, s, Vt = np.linalg.svd(X, full_matrices=True)
        S = np.zeros(X.shape)
        np.fill_diagonal(S, s)
        return np.round(U, 3), np.round(S, 3), np.round(Vt.T, 3)

    def visualize_svd(self, X, out=False):
        U, S, V = self.svd(X)
        all_ = np.r_[X.flatten(order='C'), U.flatten(order='C'),
                     S.flatten(order='C'), V.flatten(order='C')]

        if out:
            out_pdf = 'SVD.pdf'
            pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)

            if self.norm==True:

                plot_list_peak = {"Normalized Peak sum contribution per residue obtained by clustering method": X,
                                  "SVD Left sigular vector matrix": U,
                                  "SVD Right sigular vector matrix": V.T, "Diagonal matrix": S}

                plot_list_density = {"Normalized Electron sum contribution per residue obtained by clustering method": X,
                                     "SVD Left sigular vector matrix": U,
                                     "SVD Right sigular vector matrix": V.T, "Diagonal matrix": S}

            else:
                plot_list_peak = {" Peak sum contribution per residue obtained by clustering method": X,
                                  "SVD Left sigular vector matrix": U,
                                  "SVD Right sigular vector matrix": V.T, "Diagonal matrix": S}

                plot_list_density = {" Electron sum contribution per residue obtained by clustering method": X,
                                     "SVD Left sigular vector matrix": U,
                                     "SVD Right sigular vector matrix": V.T, "Diagonal matrix": S}

            if self.map_type == "peak":
                plot_list_dict = plot_list_peak
            elif self.map_type == "density":
                plot_list_dict = plot_list_density

            for k, v in plot_list_dict.items():
                fig, ax = plt.subplots()
                x_label = self.df.index.to_list()
                y_label = [X for X in range(1, len(self.data_dict_updated) + 1)]
                if k == list(plot_list_dict.keys())[0]:
                    fig, ax = plt.subplots(figsize=(20, 10))
                    ax = sns.heatmap(v.T, cmap=self.cmap, linewidths=.3, vmin=X.max(), vmax=X.min(),
                                     cbar_kws={"use_gridspec": "False", "shrink": 0.30, "location": "top"})
                    ax.set(xlabel='residue number', ylabel='Time')
                    # ax.set_xticklabels(x_label, rotation=90, size=10, ha='center')
                    # ax.set_yticklabels(y_label, rotation=90, size=10, ha='center')

                    if len(v.T) == len(y_label):
                        ax.set_xticks(range(len(x_label)))
                        ax.set_xticklabels(x_label, rotation=90, ha='center')
                        ax.set_yticks(range(len(y_label)))
                        ax.set_yticklabels(y_label, rotation=90, ha='center')


                elif k == list(plot_list_dict.keys())[1]:
                    fig, ax = plt.subplots(figsize=(30, 15))
                    ax = sns.heatmap(v, cmap=self.cmap, linewidths=.3,
                                     cbar_kws={"use_gridspec": "False", "shrink": 0.30, "location": "top"})
                    ax.set(xlabel='residue number', ylabel='residue number')

                    ax.set_yticks(np.arange(len(x_label)))
                    ax.set_yticklabels(x_label, rotation=360, ha='center')
                    ax.set_xticks(np.arange(len(x_label)))
                    ax.set_xticklabels(x_label, rotation=90, ha='center')
                    plt.xticks(ha='center')
                    plt.yticks(ha='center')
                elif k == list(plot_list_dict.keys())[2]:
                    fig, ax = plt.subplots(figsize=(30, 15))
                    ax = sns.heatmap(v, cmap=self.cmap, linewidths=.3,
                                     cbar_kws={"use_gridspec": "False", "shrink": 0.30, "location": "top"})
                    ax.set(xlabel='Time', ylabel='Time')

                    ax.set_yticks(np.arange(len(y_label)))
                    ax.set_yticklabels(y_label, rotation=360, ha='center')
                    # ax.set_xticks(np.arange(len(y_label)))
                    # ax.set_xticklabels(y_label, rotation=90, size=10)
                    plt.xticks(np.arange(len(y_label)), y_label, ha='center', rotation=90)
                    plt.xticks(ha='center')
                    plt.yticks(ha='center')

                elif k == list(plot_list_dict.keys())[3]:
                    fig, ax = plt.subplots(figsize=(30, 15))
                    ax = sns.heatmap(v, cmap=self.cmap, linewidths=.3, vmin=X.max(), vmax=X.min(),
                                     cbar_kws={"use_gridspec": "False", "shrink": 0.30, "location": "top"})
                    ax.set(xlabel='dimension', ylabel='dimension')
                    plt.xticks(ha='center')
                    plt.yticks(ha='center')
                ax.set_aspect("equal")
                plt.title(k)
                # plt.rcParams.update({'font.size': 12})
                pdf.savefig(fig)
            pdf.close()
        return X, U, S, V

    def analysis_of_singular_vectors(self):
        out_pdf = 'SVD_projection_1.pdf'
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        x_label = [str(int(x)) for x in self.df.index.to_list()]
        y_label = [X for X in range(1, len(self.data_dict_updated) + 1)]
        chunk_size = self.chunk_size()
        size = [x for x in range(len(self.residue_to_plot) - 1)] if self.residue_to_plot is not None and len(
            self.residue_to_plot) < 4 else [x for x in range(len(self.data_frame_list)-1)]

        for j in size:
            X_j = self.S[j, j] * self.U[:, j][:, None] @ self.V[:, j][None, :]
            X_chunks = np.array_split(X_j, chunk_size, axis=0)
            X_label = np.array_split(np.array(x_label), chunk_size, axis=0)

            for i, chunk in enumerate(X_chunks):
                fig, axs = plt.subplots(figsize=(30, 15))
                ax = sns.heatmap(chunk.T, cmap=self.cmap, linewidths=.3, vmin=self.X.max(), vmax=self.X.min(),
                                 cbar_kws={"use_gridspec": "False", "shrink": 0.80, "location": "right"})
                # ax.set_yticks([X for X in range(1, len(y_label) + 1)])
                # ax.set_xticks([X for X in range(1 + i * chunk_size, min((i + 1) * chunk_size + 1, len(x_label) + 1))])
                ax.set_xticklabels(list(X_label[i]), rotation=90, size=24)
                ax.set_yticklabels(y_label, rotation=0, size=24)
                title_txt = '$s_' + str(j + 1) + 'u_' + str(j + 1) + 'v_' + str(j + 1) + '^T$'
                # ax.set(xlabel='Residue', ylabel='Time')
                ax.set_xlabel('Residue', labelpad=20, fontsize=28)
                ax.set_ylabel('Time', labelpad=20, fontsize=28)
                ax.set_title(title_txt + f' (chunk {i + 1})', fontsize=28,pad=20)
                cbar = ax.collections[0].colorbar
                cbar.ax.set_ylabel(r' Normalized Peak sum $\rho$', fontsize=24, rotation=270, labelpad=30)
                cbar.ax.tick_params(labelsize=24)
                ax.set_aspect("equal")
                # plt.rcParams.update({'font.size': 24})
                pdf.savefig(fig)
        pdf.close()

    def chunk_size(self):
        if self.residue_to_plot is not None:
            if len(self.residue_to_plot) > 50 and len(self.residue_to_plot) < 100:
                chunk_size = 2
            elif len(self.residue_to_plot) > 100 and len(self.residue_to_plot) < 200:
                chunk_size = 4
            elif len(self.residue_to_plot) > 200 and len(self.residue_to_plot) < 300:
                chunk_size = 8
            else:
                chunk_size = 1
        elif len(self.filtered.index) > 50 and len(self.filtered.index) < 100:
            chunk_size = 2
        elif len(self.filtered.index) > 100 and len(self.filtered.index) < 200:
            chunk_size = 4
        elif len(self.filtered.index) > 200 and len(self.filtered.index) < 300:
            chunk_size = 6
        elif len(self.filtered.index) > 300:
            chunk_size = 8
        else:
            chunk_size = 1

        return chunk_size

    def analysis_of_original_vectors(self,norm=True):
        if self.map_type == "peak":
            title_name = "peak"
        elif self.map_type == "density":
            title_name = "electron density"

        out_pdf = 'SVD_original.pdf'
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        x_label = [str(int(x)) for x in list(self.df.index)]
        y_label = [X for X in range(1, len(self.data_dict_updated) + 1)]
        chunk_size = self.chunk_size()

        X_j = self.X
        X_chunks = np.array_split(X_j, chunk_size, axis=0)
        X_label = np.array_split(np.array(x_label), chunk_size, axis=0)

        for i, chunk in enumerate(X_chunks):
            fig, axs = plt.subplots(figsize=(20, 10))
            if norm ==True:
                ax = sns.heatmap(chunk.T, cmap=self.cmap, linewidths=.3,annot=False,vmax=1,vmin=0,
                             cbar_kws={"use_gridspec": "False", "location": "right", "shrink": 0.7},annot_kws={"size":30})
            else:
                ax = sns.heatmap(chunk.T, cmap=self.cmap, linewidths=.3,annot=False,
                             cbar_kws={"use_gridspec": "False", "location": "right", "shrink": 0.7},annot_kws={"size":30})
            # ax.set_yticks([X for X in range(1, len(y_label) + 1)])
            # ax.set_xticks([X for X in range(1 + i * chunk_size, min((i + 1) * chunk_size + 1, len(x_label) + 1))])
            ax.set_xticks([x for x in range(len(X_label[i]))])
            ax.set_xticklabels(X_label[i], rotation=90, fontsize=30)
            # print(X_label[i])
            # print(len(X_label[i]))
            ax.set_yticklabels(list(self.filtered.columns), rotation=0, fontsize=30)
            ax.set(xlabel='Residue', ylabel='Time')
            ax.set_xlabel('Residue', fontsize=32, labelpad=20)
            ax.set_ylabel('Time', fontsize=32, labelpad=20)
            cbar = ax.collections[0].colorbar
            cbar.ax.set_ylabel(r' Normalized %s sum $\rho$'%title_name, fontsize=32, rotation=270, labelpad=35)
            cbar.ax.tick_params(labelsize=30)
            if norm:
                ax.set_title(f"Normalized %s sum contribution per residue obtained by clustering method" % title_name,pad=20, fontsize=32)
            else:
                ax.set_title(f" %s sum contribution per residue obtained by clustering method" % title_name,pad=20, fontsize=32)


            ax.set_aspect("equal")
            # plt.rcParams.update({'font.size': 24})
            plt.tight_layout()
            pdf.savefig(fig)

        pdf.close()

    def correlattion_plot(self):
        out_pdf = 'correlation_plot.pdf'
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        dict_to_plot = {"Time": self.filtered.corr(), "Residue": self.filtered.T.corr()}
        for k, v in dict_to_plot.items():
            if k == "Time":
                fig, axs = plt.subplots(figsize=(10, 10))
                # corr_01 = (self.filtered.corr() + 1) / 2
                corr_01 = self.filtered.corr()
                ax = sns.heatmap(corr_01, cmap='Spectral',
                                 cbar_kws={"use_gridspec": "False", "shrink": 0.80, "location": "right",'label': 'Pearson correlation coefficient'})
                ax.set_title("Time correlation heatmap",fontsize=14,pad=20)
                x_label = list(self.filtered.columns) # [f"Time_%s" % i for i in range(1, self.filtered.corr().shape[0] + 1)]
                ax.set_xticklabels(x_label, rotation=90,fontsize=12)
                ax.set_yticklabels(x_label, rotation=0,fontsize=12)
                ax.set(xlabel='Time', ylabel='Time')
                ax.set_xlabel('Time', fontsize=14)
                ax.set_ylabel('Time', fontsize=14)
                cbar = ax.collections[0].colorbar
                cbar.ax.tick_params(labelsize=12)
                cbar.ax.set_ylabel('Pearson correlation coefficient', fontsize=12,rotation=270,labelpad=20)
                plt.tight_layout()


            elif k == "Residue":
                fig, axs = plt.subplots(figsize=(20, 20))
                # corr_02 = (self.filtered.T.corr() + 1) / 2
                corr_02 = self.filtered.T.corr()
                # df_t = pd.DataFrame(corr_02,index=[int(x) for x in list(self.filtered.index)],columns=[int(x) for x in list(self.filtered.index)])
                # df_t.to_csv("corr.csv")
                ax = sns.heatmap(corr_02, cmap='Spectral',vmin=-1,vmax=1,
                                 cbar_kws={"use_gridspec": "False", "shrink": 0.80, "location": "right",'label': 'Pearson correlation coefficient'})
                ax.set_title('Residue-Residue correlation heatmap observed for different time period',fontsize=32,pad=20)
                n_cols, n_rows = self.filtered.T.corr().shape
                ax.set_xticks([i + 0.5 for i in range(n_cols)])
                ax.set_yticks([i + 0.5 for i in range(n_rows)])
                x_label_res = [int(x) for x in list(self.filtered.index)]
                ax.set_xticklabels(x_label_res, rotation=90,fontsize=30)
                ax.set_yticklabels(x_label_res, rotation=0,fontsize=30)
                ax.set_xlabel('Residue', labelpad=20,fontsize=32)
                ax.set_ylabel('Residue', labelpad=20,fontsize=32)
                cbar = ax.collections[0].colorbar
                cbar.ax.tick_params(labelsize=24)
                cbar.ax.set_ylabel('Pearson correlation coefficient', fontsize=32,rotation=270,labelpad=25)
                plt.tight_layout()
                # ax.set(xlabel='Residue', ylabel='Residue')
            ax.set_aspect("equal")
            # plt.rcParams.update({'font.size': 12})
            pdf.savefig(fig)
        pdf.close()

    def variance_plot(self):

        out_pdf = 'variance_plot.pdf'
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        fig, ax = plt.subplots(figsize=(20, 5))
        S = self.S.diagonal()
        var_explained = np.round(S ** 2 / np.sum(S ** 2), decimals=3)
        ax = sns.barplot(x=list(range(1, len(var_explained) + 1)),
                         y=var_explained, color="lightgreen", alpha=0.6)

        plt.yticks(fontsize=24)
        plt.xticks(fontsize=24)

        plt.xlabel('singular vectors', labelpad=20,fontsize=24)
        plt.ylabel('Percent Variance Explained', labelpad=20,fontsize=24)
        ax.set_title("Percentage of variance explained by each singular vector",fontsize=24,pad=20)
        # plt.rcParams.update({'font.size': 24})
        plt.tight_layout()
        pdf.savefig(fig)
        pdf.close()

    def linear_relational_plot(self, font=32,norm=True):
        if self.map_type == "density":
            if norm:
                title = r"Variation of normalized electron sum $\rho$ over time"
                y_label = r'Normalized electron density sum $\rho$'
            else:
                title = r"Variation of electron sum $\rho$ over time"
                y_label = r' Electron sum $\rho$'
        elif self.map_type == "peak":
            if norm:
                title = r"Variation of normalized peak sum $\rho$ over time"
                y_label = r' Normalized peak sum $\rho$'
            else:
                title = r"Variation of peak sum $\rho$ over time"
                y_label = r'Peak sum $\rho$'

        out_pdf = "one_dimensional_linear_plot.pdf"
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        self.filtered.index = self.filtered.index.astype(int)
        sub_results = [self.filtered[i:i + 50] for i in range(0, len(self.filtered), 50)]
        for sub_result in sub_results:
            fig, ax = plt.subplots(figsize=(30, 15))
            cmap = plt.cm.get_cmap('Spectral', len(self.filtered.columns))
            colors = cmap(np.linspace(0, 1, len(self.filtered.columns)))
            for i, x in enumerate(list(sub_result.columns)):
                sub_result[x].plot(kind='bar', color=cmap(i), width=.5, alpha=0.5, legend=True,fontsize=32)
            plt.xlabel('Residue', fontsize=32,labelpad=20)
            plt.ylabel(y_label, fontsize=32,labelpad=20)
            ax.set_title(title,pad=20,fontsize=32)
            if self.norm == True:
                plt.ylim(0,1.1)
            plt.xticks(fontsize=font)
            plt.yticks(fontsize=font)
            ax.tick_params(axis='x', which='major', pad=15)
            # plt.rcParams.update({'font.size': font})
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': 32})
            plt.subplots_adjust(right=0.8)
            plt.tight_layout()
            pdf.savefig(fig)
        pdf.close()

    def part_residue(self):
        residues = []
        for part in self.residue.split(","):
            if "-" in part:
                start, end = map(int, part.split("-"))
                residues += list(range(start, end + 1))
            else:
                residues.append(int(part))

        residues_list = [x for x in residues if x in list(self.df.index)]

        return residues_list


def main():
    import argparse
    parser = argparse.ArgumentParser(prog='residem_svd_csv', formatter_class=argparse.RawTextHelpFormatter,
                                     description=f"""It analyse csv for electron density sum and plots \n
                                     pearson correlation analysis, one dimensional density bar plot and SVD heat map plot
        Few example to run the code: \n\n
        1. residem_svd_csv -r 5b6v.pdb -f *.csv -cmap Spectral \n\n
        2. residem_svd_csv -r 5b6v.pdb -f *.csv
        3. residem_svd_csv -r 5b6v.pdb -f *.csv -n 182,216,300
        4. residem_svd_csv -r 5b6v.pdb -f *.csv -i 10
        5. residem_svd_csv -r 5b6b.pdb -f *.csv -i 10 -n 182,216-218,300
        6. residem_svd_csv -r 5b6b.pdb -f *.csv -i 10 -n 182,216-218,300 -d CA/SC/CM
        7. residem_svd_csv -r 5b6b.pdb -f *.csv -i 10 -n 182,216-218,300 -d CA/SC/CM -t density/peak --cmap Spectral --ew 0.015 
        
        8. residem_svd_csv -r 5b6b.pdb --dict time_and_file_name.txt -n 172-614 -t density --ew 0.05
        
        For other color please visit the below website\n
        "https://matplotlib.org/stable/tutorials/colors/colormaps.html"
        "https://matplotlib.org/stable/gallery/color/named_colors.html"
        only css colors can be used of --ncmap
    """)
    parser.add_argument("-f", "--file", action='store', nargs='+',
                        help="data frame for cluster obtained from cluster analysis ")
    parser.add_argument("-c", "--cmap", action='store', type=str,
                        help="cmap colour")
    parser.add_argument("-s", "--ncmap", action='store', type=str,
                        help=" network cmap colour")
    parser.add_argument("-k", "--closeness", action='store', type=float,
                        help="subplot node closeness")
    parser.add_argument("-i", "--cutoff", action='store', type=int,
                        help="cutoff score for data frame")
    parser.add_argument("-n", "--residue", action='store', type=str,
                        help="list of residue list")
    parser.add_argument("-r", "--pdb", action="store", type=str,
                        help="pdb file")
    parser.add_argument("-d", "--dist", action="store", type=str,
                        help="calculation of distance type")
    parser.add_argument("-l", "--low", action="store", type=float,
                        help="electron density sum cutoff threshold low for plotting node size")
    parser.add_argument("-m", "--medium", action="store", type=float,
                        help="electron density sum cutoff threshold medium for plotting node size")
    parser.add_argument("-t", "--map_type", action="store", type=str,
                        help="It either chooses the normalized peak sum or electron density sum, input option, peak or density")
    parser.add_argument("-z","--norm", action="store",type=str,
                        help="if data has to be normalized")
    parser.add_argument( "--dict", action="store", type=str,
                        help="dictionary of input file")
    parser.add_argument( "--ew", action="store", type=float,
                        help="edge with size")


    args = parser.parse_args()

    edge_width_size = args.ew if args.cmap is not None else 0.150
    cmap = args.cmap if args.cmap is not None else "Spectral"
    ncmap = args.ncmap if args.ncmap is not None else None
    cutoff_value = args.cutoff if args.cutoff is not None else 0
    residue = args.residue if args.residue is not None else None
    distance = args.dist if args.dist is not None else "CA"
    low = args.low if args.low is not None else 0.05
    closeness = args.closeness if args.closeness is not None else 0.05
    time = args.dict if args.dict is not None else None
    high = args.medium if args.medium is not None else 0.1
    map_type = args.map_type if args.map_type in ["peak", "density"] else "peak"
    norm = eval([x for x in ["false","true"] if args.norm.lower() in x][0].capitalize()) if args.norm is not None else True

    if args.file is not None:
        data_files = [x for x in args.file]
        data_files.sort()
    elif args.dict is not None:
        text = json.load(open(args.dict))
        data_files = list(text.values())
        time = list(text.keys())


    SVD_map = SVD(data_files, cmap=cmap, residue=residue, cutoff_value=cutoff_value, map_type=map_type,norm=norm, time = time )
    if len(data_files) >= 2:
        SVD_map.analysis_of_singular_vectors()
        SVD_map.correlattion_plot()
        SVD_map.variance_plot()
        SVD_map.analysis_of_original_vectors(norm=norm)
    if len(data_files) >= 2:
        SVD_map.linear_relational_plot(norm=norm)
    else:
        SVD_map.linear_relational_plot(norm=norm)
    if args.pdb is not None:
        print("*"*10)
        print("Plotting network graph")
        Network_map = Network(args.pdb, SVD_map.nof, distance_type=distance,
                              low_threshold=low, medium_threshold=high, map_type=SVD_map.map_type,subgraph_residue=SVD_map.residue_list,map_colour=ncmap,closeness=closeness,norm=norm,edge_width=edge_width_size)
        Network_map.plot_network(label=True)
        Network_map.plot_network(label=False)
        Network_map.subplot()


if __name__ == '__main__':
    main()
