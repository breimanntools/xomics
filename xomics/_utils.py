#! /usr/bin/python3
"""
Config with folder structure
"""
import os
import platform
from pathlib import Path
import itertools

import matplotlib as mpl
import numpy as np
import pandas as pd


# Helper Function
import seaborn as sns
from matplotlib import pyplot as plt


def _folder_path(super_folder, folder_name):
    """Modification of separator (OS depending)"""
    path = os.path.join(super_folder, folder_name + SEP)
    return path


# Folder
SEP = "\\" if platform.system() == "Windows" else "/"
FOLDER_PROJECT = str(Path(__file__).parent.parent).replace('/', SEP) + SEP
FOLDER_DATA = _folder_path(FOLDER_PROJECT, 'data')

COL_E_SCORE = "E-score"
COL_P_SCORE = "P-score"
COL_C_SCORE = "C-score"
COL_CS_STD = "CS_STD"
COL_PE_MEAN = "PE-mean"

STR_PVAL = "-log10 p-value"
STR_FC = "log2 FC"

# Volcano default colors
COLOR_UP = "firebrick"
COLOR_DOWN = "dodgerblue"
COLOR_NOT_SIG = "gray"
COLOR_TH = "black"
COLOR_GEM = "#69C2CA"


# Plotting settings
DPI = 300
FIG_FORMAT = "pdf"
ARGS_SAVE = dict(dpi=DPI, bbox_inches="tight")
LEGEND_FONTSIZE = 16


# General helper functions
def is_nested_list_with_strings(lst):
    for element in lst:
        if isinstance(element, list):
            if all(isinstance(sub_element, str) for sub_element in element):
                return True
    return False


def flatten_list(list_in, sep=","):
    """Flatten list the list and provide unique items"""
    if isinstance(list_in, pd.Series):
        list_in = list_in.to_list()
    elif is_nested_list_with_strings(list_in):
        list_in = list(itertools.chain.from_iterable(list_in))
    f_ = lambda l: [x.strip() for x in l] if type(l) is list else (l.strip() if type(l) is str else np.NaN)
    f = lambda x: x.split(sep) if sep in x and type(x) is str else [x]
    chained_list = list(itertools.chain.from_iterable([f_(f(x)) for x in list_in]))
    unique_items = []
    for item in chained_list:
        if item not in unique_items:
            unique_items.append(item)
    return unique_items


# Documentation function
def print_red(input_str):
    """Prints the given string in red text."""
    print(f"\033[91m{input_str}\033[0m")


# Plotting function
def plot_settings(fig_format="pdf", verbose=False, grid=False, grid_axis="y", font_scale=0.7,
                  change_size=True, weight_bold=True, adjust_elements=True,
                  short_ticks=False, no_ticks=False,
                  no_ticks_y=False, short_ticks_y=False, no_ticks_x=False, short_ticks_x=False):
    """General plot settings"""
    mpl.rcParams.update(mpl.rcParamsDefault)
    # Set embedded fonts in PDF
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["pdf.fonttype"] = 42
    if verbose:
        print(plt.rcParams.keys)    # Print all plot settings that can be modified in general
    if not change_size:
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = "Arial"
        font = {'family': 'Arial'}
        mpl.rc('font', **font)
        return
    sns.set_context("talk", font_scale=font_scale)  # Font settings https://matplotlib.org/3.1.1/tutorials/text/text_props.html
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Arial"
    if weight_bold:
        plt.rcParams["axes.labelweight"] = "bold"
        plt.rcParams["axes.titleweight"] = "bold"
    else:
        plt.rcParams["axes.linewidth"] = 1
        plt.rcParams["xtick.major.width"] = 0.8
        plt.rcParams["xtick.minor.width"] = 0.6
        plt.rcParams["ytick.major.width"] = 0.8
        plt.rcParams["ytick.minor.width"] = 0.6
    if short_ticks:
        plt.rcParams["xtick.major.size"] = 3.5
        plt.rcParams["xtick.minor.size"] = 2
        plt.rcParams["ytick.major.size"] = 3.5
        plt.rcParams["ytick.minor.size"] = 2
    if short_ticks_x:
        plt.rcParams["xtick.major.size"] = 3.5
        plt.rcParams["xtick.minor.size"] = 2
    if short_ticks_y:
        plt.rcParams["ytick.major.size"] = 3.5
        plt.rcParams["ytick.minor.size"] = 2
    if no_ticks:
        plt.rcParams["xtick.major.size"] = 0
        plt.rcParams["xtick.minor.size"] = 0
        plt.rcParams["ytick.major.size"] = 0
        plt.rcParams["ytick.minor.size"] = 0
    if no_ticks_x:
        plt.rcParams["xtick.major.size"] = 0
        plt.rcParams["xtick.minor.size"] = 0
    if no_ticks_y:
        plt.rcParams["ytick.major.size"] = 0
        plt.rcParams["ytick.minor.size"] = 0

    plt.rcParams["axes.labelsize"] = 17 #13.5
    plt.rcParams["axes.titlesize"] = 16.5 #15
    if fig_format == "pdf":
        mpl.rcParams['pdf.fonttype'] = 42
    elif "svg" in fig_format:
        mpl.rcParams['svg.fonttype'] = 'none'
    font = {'family': 'Arial', "weight": "bold"} if weight_bold else {"family": "Arial"}
    mpl.rc('font', **font)
    if adjust_elements:
        # Error bars
        plt.rcParams["errorbar.capsize"] = 10   # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.errorbar.html
        # Grid
        plt.rcParams["axes.grid.axis"] = grid_axis  # 'y', 'x', 'both'
        plt.rcParams["axes.grid"] = grid
        # Legend
        plt.rcParams["legend.frameon"] = False
        plt.rcParams["legend.fontsize"] = "medium" #"x-small"
        plt.rcParams["legend.loc"] = 'upper right'  # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
