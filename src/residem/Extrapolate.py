import mmtbx.f_model
import math
from cctbx import french_wilson
from cctbx import miller, maptbx
import numpy as np
import sys
from residem.master import DataInput
from cctbx.array_family import flex
from mmtbx import map_tools
import os
from scipy.stats import pearsonr
from collections import defaultdict
from residem.multiprocess_func import make_parallel
import tqdm
from libtbx.utils import null_out
from iotbx.map_model_manager import map_model_manager
import pandas as pd
import matplotlib.pyplot as plt
from cctbx.xray import structure
import iotbx.map_tools
from scipy.signal import find_peaks
from scipy.optimize import minimize, least_squares
from iotbx.data_manager import DataManager
import pathlib
import cctbx.miller
import seaborn as sns
from functools import reduce
import operator
from scipy.optimize import curve_fit, minimize_scalar
from sklearn.metrics import r2_score
import matplotlib.backends.backend_pdf
import warnings
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from residem.scaling import ScaleNativeDerivative
from kneed import KneeLocator
from residem.logger_set import setup_logger
import sympy as sp

matplotlib.use('Agg')
font = {'weight': 'normal', 'size': 8}

matplotlib.rc('font', **font)


class ExtrapoledMap(ScaleNativeDerivative):
    def __init__(self, reference_pdb, reference_mtz, Triggered_mtz, reference_label=None,
                 triggered_label=None,
                 high_resolution=None, low_resolution=None,fcalclabel=None,sigma=3.0, weight="default", grid_step=0.5, scale="iso",
                 write=True,
                 path=None, alpha=None,logger=None):
        super().__init__(reference_pdb, reference_mtz, Triggered_mtz, reference_label, triggered_label,fcalclabel,
                         high_resolution, low_resolution, weight=weight, grid_step=grid_step, scale=scale,
                         write=write,
                         path=path,alpha=alpha,logger=logger)

        self.path = path if path is not None else os.getcwd()
        os.makedirs(self.path, exist_ok=True)
        self.Plot_R_iso_CCiso(self.path)
        # wilson plot
        self.Plot_wilson(self.path)
        self.Plot_weight(self.path)

        self.F_cal_extra = defaultdict(int)
        self.weight = None if weight is None else weight.lower() if weight.lower() in ["ren", "ursby","hekstra"] else None
        self.sigma = sigma
        self.write = True
        if write:
            self.mtz_dataset_object = self.mtz_dataset(f"{os.path.join(self.path, 'Difference_map_weighted_all.mtz')}")

            ## map coefficient
            self.FOFODF = self.map_coefficient_cal(self.delta_F,sigma=flex.double(self.sigma_df), coefficient="FO", write=self.write,
                                               name=f"%s/F_obs_minus_F_obs.ccp4" % self.path, map_type="diff")
            self.FOFODF_K = self.map_coefficient_cal(self.wDF,sigma=flex.double(self.sigma_ren_df), coefficient="FO", write=self.write,
                                                 name=f"%s/F_obs_minus_F_obs_ren_weight.ccp4" % self.path,
                                                 map_type="diff")
            self.FOFODF_Q = self.map_coefficient_cal(self.rDF,sigma=flex.double(self.sigma_ursby_df), write=self.write, coefficient="FO",
                                                 name=f"%s/F_obs_minus_F_obs_ursby_weight.ccp4" % self.path,
                                                 map_type="diff")
            self.FOFODF_H = self.map_coefficient_cal(self.hDF,sigma=flex.double(self.sigma_hekstra_df), write=self.write, coefficient="FO",
                                                 name=f"%s/F_obs_minus_F_obs_hekstra_weight.ccp4" % self.path,
                                                 map_type="diff")



