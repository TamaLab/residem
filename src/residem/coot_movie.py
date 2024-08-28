# I have a pdb model, a difference map and a folder containing a series of difference maps and 2fofc maps. The objective is to create a movie of the difference maps and 2fofc maps where each map a image is created and then the images are combined into a movie.

import coot
import os
from glob import glob


class CootMovie:
    def __init__(self, pdb, diff_map, data_folder="Data_folder_0", diff_map_folder_path="Extrapolate_diff_map", extrapolated_folder_path="Extrapolate_map"):
        self.pdb = pdb
        self.diff_map = diff_map
        self.data_folder = data_folder
        self.diff_map_folder = os.path.join(data_folder, diff_map_folder_path)
        self.two_fofc_folder = os.path.join(
            data_folder, extrapolated_folder_path)
        self.dff_map_list = glob(os.path.join(self.diff_map_folder, "*.map"))
        self.two_fofc_list = glob(os.path.join(self.two_fofc_folder, "*.map"))
        self.diff_map_dict = self.sort_ccp4_map(self.dff_map_list)
        self.two_fofc_dict = self.sort_ccp4_map(self.two_fofc_list)

    # function to load the pdb

    def load_pdb(self):
        coot.handle_read_draw_molecule(self.pdb)
        # coot.read_pdb(self.pdb)

    # function to load the difference map

    def load_diff_map(self):
        coot.handle_read_ccp4_map(self.diff_map, 1)
        coot.set_map_is_difference_map(1, 1)

    # function to save image of difference map
    def save_diff_map(self, name="diff_map"):
        try:
            os.makedirs("%s/images/%s" % (self.data_folder, name))
        except:
            pass
        # coot.set_view_quaternion(
        #     self.view_quaternion[0], self.view_quaternion[1], self.view_quaternion[2], self.view_quaternion[3])
        # coot.set_zoom(15)
        # coot.graphics_draw()
        coot.screendump_image("%s/images/%s/%s.png" %
                              (self.data_folder, name, name))

    # function to load the difference map folder
    def sort_ccp4_map(self, list_files):
        occupancy_list = [(int(file.split("/")[-1].split(".map")
                               [0].split("_")[-1]), file) for file in list_files]

        occupancy_list.sort(key=lambda x: x[0])
        # ccp4_map = {k: to_sort_pairs[k]
        #             for k in sorted(list(to_sort_pairs.keys()))}
        ccp4_map = dict(occupancy_list)
        return ccp4_map

    @staticmethod
    def map_molecule_list():

        map_list = []
        for i in range(coot.graphics_n_molecules()):
            if coot.is_valid_map_molecule(i) == 1:
                map_list.append(i)
        return map_list

    def create_images_diff(self, name="extra_diff_map"):

        try:
            os.makedirs("%s/images/%s" % (self.data_folder, name))
        except:
            pass

        to_make = self.diff_map_dict if name == "extra_diff_map" else self.two_fofc_dict
        list_of_maps = self.map_molecule_list()
        for k, v in to_make.items():
            coot.handle_read_ccp4_map(v, 1) if name == "extra_diff_map" else coot.handle_read_ccp4_map(
                v, 0)
            diff_model_number = list_of_maps[-1]
            if coot.molecule_name(diff_model_number) == v:
                diff_model_number = list_of_maps[-1]
            else:
                for i in range(0, 200):
                    if coot.molecule_name(i) == v:
                        diff_model_number = i
                        break
            print(diff_model_number, v)
            print("list_map", list_of_maps)
            if name == "extra_diff_map":
                coot.set_map_is_difference_map(diff_model_number, 1)
                coot.set_contour_level_in_sigma(diff_model_number, 3.0)
            else:
                coot.set_contour_level_in_sigma(diff_model_number, 1.0)
                coot.set_map_colour(diff_model_number, 0, 0, 0.7)
            coot.screendump_image(
                "%s/images/%s/%s_occu_%s.png" % (self.data_folder, name, name, k))
            coot.close_molecule(diff_model_number)
        # coot_movie.create_images_diff("extra")


def run(pdb, diff_map, data_folder=os.getcwd()):
    pdb = os.path.join(os.getcwd, pdb)
    diff_map = os.path.join(
        os.getcwd(), diff_map)

    coot_movie = CootMovie(pdb, diff_map, data_folder=data_folder)
    coot_movie.load_pdb()
    coot_movie.load_diff_map()
    return coot_movie


coot_movie = run()
