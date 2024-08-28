import gobject
import gtk
from glob import glob
import os
import time
import urllib
import numpy as np
import coot
import json
import pygtk
import threading
pygtk.require('2.0')

# urlpath = "https://raw.githubusercontent.com/hypowergravity/resmap/main/src/resmap/resmap.svg?token=GHSAT0AAAAAAB4OHG72DI65AV7MM76UQGEIY7KIO5A"
# image_filename = "image.svg"
# urllib.urlretrieve(urlpath, image_filename)


pathdirc = os.getcwd()


class CootMovie:
    def __init__(self, pdb=None, diff_map=None, data_folder=None, diff_map_folder_path="Extrapolate_diff_map", extrapolated_folder_path="Extrapolate_map",number =0):
        self.pdb = pdb
        self.diff_map = diff_map
        self.data_folder = data_folder
        self.diff_map_folder = os.path.join(data_folder, diff_map_folder_path)
        self.two_fofc_folder = os.path.join(
            data_folder, extrapolated_folder_path)
        self.dff_map_list = glob(os.path.join(self.diff_map_folder, "*.ccp4"))
        self.two_fofc_list = glob(os.path.join(self.two_fofc_folder, "*.ccp4"))
        self.diff_map_dict = self.sort_ccp4_map(self.dff_map_list)
        self.two_fofc_dict = self.sort_ccp4_map(self.two_fofc_list)

        self.directory = self.set_directory()
        self.I = number

    def set_directory(self):
        I = 0
        name = "images"
        exists = True
        while exists:
            path_folder =os.path.join(self.data_folder, name)
            if not os.path.exists( path_folder + "_%s" % I):
                # os.mkdir(name + "_%s" % I)
                exists = False
            else:
                I = I + 1

        name1 = name + "_%s" % I
        # name_path = os.getcwd() + f"/%s" % name1
        return name1

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
            os.makedirs("%s/%s/%s" % (self.data_folder, self.directory, name))
        except:
            pass
        # coot.set_view_quaternion(
        #     self.view_quaternion[0], self.view_quaternion[1], self.view_quaternion[2], self.view_quaternion[3])
        # coot.set_zoom(15)
        # coot.graphics_draw()
        coot.screendump_image("%s/%s/%s/%s.png" %
                              (self.data_folder,self.directory, name, name))

    # function to load the difference map folder
    def sort_ccp4_map(self, list_files):
        occupancy_list = [(int(file.split("/")[-1].split(".ccp4")
                               [0].split("_")[1]), file) for file in list_files]

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
            os.makedirs("%s/%s/%s" % (self.data_folder,self.directory, name))
        except:
            pass

        to_make = self.diff_map_dict if name == "extra_diff_map" else self.two_fofc_dict
        list_of_maps = self.map_molecule_list()
        if name == "extra_diff_map":
            for map in list_of_maps:
                coot.set_map_displayed(map, 0)
        else:
            for map in list_of_maps:
                coot.set_map_displayed(map, 1)
        for k, v in to_make.items():
            coot.handle_read_ccp4_map(v, 1) if name == "extra_diff_map" else coot.handle_read_ccp4_map(
                v, 0)
            diff_model_number = list_of_maps[-1]
            if coot.molecule_name(diff_model_number) == v:
                diff_model_number = list_of_maps[-1]
            else:

                x = 220 if self.I == 0 else (self.I+2) * 220
                for i in range(0, x):
                    if coot.molecule_name(i) == v:
                        diff_model_number = i
                        break
            # print(diff_model_number, v)
            # print("list_map", list_of_maps)
            if name == "extra_diff_map":
                coot.set_map_is_difference_map(diff_model_number, 1)
                coot.set_contour_level_in_sigma(diff_model_number, 3.0)
            else:
                coot.set_contour_level_in_sigma(diff_model_number, 1.0)
                coot.set_map_colour(diff_model_number, 0, 0, 0.7)
            coot.screendump_image(
                "%s/%s/%s/%s_occu_%s.png" % (self.data_folder,self.directory, name, name, k))
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


class AboutWindow(gtk.Window):
    def __init__(self):
        super(AboutWindow, self).__init__()

        self.set_title("About")
        self.set_default_size(200, 200)
        self.set_resizable(True)
        self.set_position(gtk.WIN_POS_CENTER)
        label = gtk.Label("Please visit the github page for documentaion")
        self.add(label)



class Coot_call:
    def __init__(self):
        coot.set_map_is_difference_map(1, True)
        coot.set_contour_level_in_sigma(1, 3.0)
        coot.set_map_colour(1, 0.3, 0.84, 0.24)


class Menu(gtk.Window):
    Coot_call()

    def on_about(self, widget):
        about_window = AboutWindow()
        about_window.show_all()

    def on_destroy(self, widget):
        gtk.hide()

    def on_file(self, widget):
        dlg = gtk.FileChooserDialog("Open..", None, gtk.FILE_CHOOSER_ACTION_OPEN,
                                    (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        response = dlg.run()
        if dlg.get_filename().endswith('.json'):
            print(dlg.get_filename())
            self.text.set_text(dlg.get_filename())
            CootFoFo(dlg.get_filename())
            self.data_file = dlg.get_filename()

            # Coot_call()
        else:
            print(dlg.get_filename() + " is not a JSON file")

        dlg.destroy()

    def Coot(self, widget):
        coot.set_map_is_difference_map(1, True)
        coot.set_contour_level_in_sigma(1, 3.0)
        coot.set_map_colour(1, 0.3, 0.84, 0.24)

    def on_button_clicked(self, widget):
        self.Coot(widget)

    def on_button_clicked_1(self, widget):
        data_file = self.data_file
        # path_parts = data_file.split(os.sep)
        # index = next(i for i, part in enumerate(
        #     path_parts) if "Data_folder_" in part)
        coot.molecule_name(1)
        new_path = "/".join(coot.molecule_name(1).split("/")[:-1])
        self.number = self.number + 1

        # new_path = os.sep.join(path_parts[:index+1])

        coot_movie = CootMovie(data_folder=new_path,number=self.number)
        coot_movie.save_diff_map()
        coot_movie.create_images_diff()
        coot_movie.create_images_diff("extra")

    def __init__(self):
        super(Menu, self).__init__()
        self.data_file = None
        self.set_title("Choose JSON File")
        self.set_default_size(200, 200)
        self.set_position(gtk.WIN_POS_CENTER)
        self.set_resizable(False)
        self.text = gtk.Entry()
        self.number = 0
        mb = gtk.MenuBar()
        menu1 = gtk.Menu()
        file = gtk.MenuItem("_File")
        file.set_submenu(menu1)
        acgroup = gtk.AccelGroup()
        self.add_accel_group(acgroup)
        open = gtk.ImageMenuItem(gtk.STOCK_OPEN)
        open.connect("activate", self.on_file)
        menu1.append(open)
        sep = gtk.SeparatorMenuItem()
        menu1.append(sep)
        mb.append(file)
        vbox = gtk.VBox(False, 2)
        vbox.pack_start(mb, False, False, 0)

        # Add the button
        # button = gtk.Button("Set map as Difference map")
        # button.set_size_request(100, 30)
        # button.connect("clicked", self.on_button_clicked)
        # vbox.pack_start(button, False, False, 0)
        #
        # self.add(vbox)
        menu2 = gtk.Menu()
        help = gtk.MenuItem("_Help")
        help.set_submenu(menu2)
        mb.append(help)
        about = gtk.MenuItem("_About")
        about.connect("activate", self.on_about)
        menu2.append(about)
        image = gtk.Image()
        # Image settings
        width = 350
        height = 200
        pixbuf = gtk.gdk.pixbuf_new_from_file(pathdirc+"/Resi-DEM_icon.png")
        # image.set_from_file(pathdirc+"/resmap.svg")
        scale_pixbuf = pixbuf.scale_simple(
            width, height, gtk.gdk.INTERP_BILINEAR)
        image.set_from_pixbuf(scale_pixbuf)
        # image.set_size_request(200, 200)
        vbox.pack_start(image, False, False, 0)
        self.add(vbox)
        button_1 = gtk.Button("Make images")
        button_1.set_size_request(100, 30)
        button_1.connect("clicked", self.on_button_clicked_1)
        # vbox.pack_start(button_1, False, False, 0)
        self.connect("destroy", gtk.main_quit)
        self.show_all()


# Coot_call()

class CootFoFo:

    def json_data(self, file):
        with open(file) as json_file:
            data = json.load(json_file)
            return data

    # close the window and quit
    def delete_event(self, widget, event, data=None):
        gtk.main_quit()
        return False

    def goto_xyz(self, widget, row, col):
        coot.set_zoom(20)
        # print(self.treeview.get_model()[row])
        try:
            coot.set_rotation_centre(float(self.treeview.get_model()[row][5]), float(self.treeview.get_model()[row][6]),
                                     float(self.treeview.get_model()[row][7]))
        except ValueError:
            coot.set_go_to_atom_chain_residue_atom_name(
                "A", int(self.treeview.get_model()[row][0]), "CA")

    def onSelectionChanged(self, path, model):
        dict_resi = {}
        for x in range(len(model)):
            dict_resi[model[x][0]] = model[x][1]
            Filtered = [k for k, v in dict_resi.items() if v == True]
            resi_list = sorted(list(map(int, Filtered)))
            np.savetxt('residue_list.txt', resi_list, newline=',', fmt='%d')

    def col1_toggled_cb(self, cell, path, model):
        model[path][1] = not model[path][1]

    def __init__(self, file):
        self.data = self.json_data(file)
        print(self.data)
        # Create a new window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)

        self.window.set_title("Isomorphous Difference peaks ")

        self.window.set_size_request(500, 400)

        self.window.connect("delete_event", self.delete_event)
        self.treestore = gtk.TreeStore(gobject.TYPE_STRING, gobject.TYPE_BOOLEAN, gobject.TYPE_BOOLEAN,
                                       gobject.TYPE_STRING, gobject.TYPE_STRING, gobject.TYPE_STRING,
                                       gobject.TYPE_STRING, gobject.TYPE_STRING, gobject.TYPE_STRING,
                                       gobject.TYPE_STRING, gobject.TYPE_STRING, gobject.TYPE_STRING,
                                       gobject.TYPE_STRING, gobject.TYPE_STRING)

        for k, v in self.data.items():
            if k == "resn":
                for x in range(len(self.data[k])):
                    parent_num = list(set(self.data[k][x].values()))[0]
                    parent = [str(parent_num), True, True, '', '',
                              '', '', '', '', '', '', '', '', '']
                    piter = self.treestore.append(None, parent)
                    for child in range(len(self.data[k][x].values())):
                        resn = list(self.data["resn"][x].values())[child]
                        resi = list(self.data["resi"][x].values())[child]
                        atom_name = list(
                            self.data["atom name"][x].values())[child]
                        X = list(self.data["X"][x].values())[child]
                        Y = list(self.data["Y"][x].values())[child]
                        Z = list(self.data["Z"][x].values())[child]
                        peak_height = list(
                            self.data["peak height"][x].values())[child]
                        peak_voxel = list(
                            self.data["number of voxel"][x].values())[child]
                        site = list(self.data["site"][x].values())[child]
                        cluster_number = list(
                            self.data["Cluster number"][x].values())[child]
                        cluster_volume = list(
                            self.data["Peak Volume of the Cluster"][x].values())[child]
                        electron_sum = list(
                            self.data["electron_sum_around_the_site"][x].values())[child]
                        child = [resn, False, False, resi, atom_name, round(X, 3), round(Y, 3), round(Z, 3),
                                 round(peak_height, 3),
                                 round(peak_voxel), site, cluster_number, cluster_volume, electron_sum]
                        self.treestore.append(piter, child)
        self.treeview = gtk.TreeView(self.treestore)

        self.cell_selection = gtk.CellRendererToggle()
        self.tvcolumn_selection = gtk.TreeViewColumn(
            "Select", self.cell_selection, active=1, visible=2)
        self.cell_selection.set_property('activatable', True)
        self.cell_selection.connect(
            'toggled', self.col1_toggled_cb, self.treestore)

        self.cell_resi = gtk.CellRendererText()
        self.tvcolumn_resi = gtk.TreeViewColumn("resn", self.cell_resi, text=0)
        self.cell_resi.set_property('editable', False)

        self.cell_resn = gtk.CellRendererText()
        self.tvcolumn_resn = gtk.TreeViewColumn(
            "chain-resi", self.cell_resn, text=3)
        self.cell_resn.set_property('editable', False)

        self.cell_atom_name = gtk.CellRendererText()
        self.tvcolumn_atom_name = gtk.TreeViewColumn(
            "atom name", self.cell_atom_name, text=4)
        self.cell_atom_name.set_property('editable', False)

        self.cell_X = gtk.CellRendererText()
        self.tvcolumn_X = gtk.TreeViewColumn(
            "X", self.cell_X, text=5, visible=1)
        self.cell_X.set_property('editable', False)

        self.cell_Y = gtk.CellRendererText()
        self.tvcolumn_Y = gtk.TreeViewColumn(
            "Y", self.cell_Y, text=6, visible=1)
        self.cell_Y.set_property('editable', False)

        self.cell_Z = gtk.CellRendererText()
        self.tvcolumn_Z = gtk.TreeViewColumn(
            "Z", self.cell_Z, text=7, visible=1)
        self.cell_Z.set_property('editable', False)

        self.cell_peak_height = gtk.CellRendererText()
        self.tvcolumn_peak_height = gtk.TreeViewColumn(
            "peak height", self.cell_peak_height, text=8)
        self.cell_peak_height.set_property('editable', False)

        self.cell_peak_volume = gtk.CellRendererText()
        self.tvcolumn_peak_volume = gtk.TreeViewColumn(
            "number of voxel", self.cell_peak_volume, text=9)
        self.cell_peak_volume.set_property('editable', False)

        self.cell_site = gtk.CellRendererText()
        self.tvcolumn_site = gtk.TreeViewColumn(
            "site", self.cell_site, text=10)
        self.cell_site.set_property('editable', False)

        self.cell_cluster_number = gtk.CellRendererText()
        self.tvcolumn_cluster_number = gtk.TreeViewColumn(
            "Cluster number", self.cell_cluster_number, text=11)
        self.cell_cluster_number.set_property('editable', False)
        self.treeview.connect("row-activated", self.goto_xyz)

        self.cell_cluster_volume = gtk.CellRendererText()
        self.tvcolumn_cluster_volume = gtk.TreeViewColumn("Peak Volume of the Cluster", self.cell_cluster_volume,
                                                          text=12)
        self.cell_cluster_volume.set_property('editable', False)

        self.cell_electron = gtk.CellRendererText()
        self.tvcolumn_electron = gtk.TreeViewColumn(
            "electron sum ", self.cell_electron, text=13)
        self.cell_electron.set_property('editable', False)

        self.treeview.append_column(self.tvcolumn_selection)
        self.treeview.append_column(self.tvcolumn_resi)
        self.treeview.append_column(self.tvcolumn_resn)
        self.treeview.append_column(self.tvcolumn_atom_name)
        # self.treeview.append_column(self.tvcolumn_X)
        # self.treeview.append_column(self.tvcolumn_Y)
        # self.treeview.append_column(self.tvcolumn_Z)
        self.treeview.append_column(self.tvcolumn_peak_height)
        self.treeview.append_column(self.tvcolumn_peak_volume)
        # self.treeview.append_column(self.tvcolumn_site)
        self.treeview.append_column(self.tvcolumn_cluster_number)
        self.treeview.append_column(self.tvcolumn_cluster_volume)
        self.treeview.append_column(self.tvcolumn_electron)

        self.treeview.set_reorderable(True)
        self.lbl = gtk.Label("Save residue for refinement")
        self.btn1 = gtk.Button(stock=gtk.STOCK_SAVE)
        self.btn1.connect("clicked", self.onSelectionChanged, self.treestore)

        gtk.Widget.set_size_request(self.btn1, 85, 15)
        box = gtk.VBox()
        vb = gtk.VBox()
        hb = gtk.HBox()
        sw = gtk.ScrolledWindow()
        sw.set_policy(gtk.POLICY_ALWAYS, gtk.POLICY_ALWAYS)
        sw.add_with_viewport(self.treeview)
        sw.add(self.treeview)
        sw.set_size_request(300, 250)
        hb.pack_start(sw, expand=True, fill=True, padding=5)
        vb.pack_start(self.lbl, True, True, 0)
        vb.pack_start(self.btn1, expand=True, fill=True, padding=5)

        box.add(hb)
        box.add(vb)
        self.window.add(box)
        self.window.show_all()


def main():
    #    CootFoFo("map_dump.json")
    # Coot_call()
    map_menu = Menu()
    gtk.main()


if __name__ == "__main__":
    # import argparse
    #
    # parser = argparse.ArgumentParser()
    # parser.add_argument("json", help="map_json_file")
    # args = parser.parse_args()
    main()
