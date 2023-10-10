"""
This is the main script for utility functions, folder structure, and constants.
Most imported modules contain checking functions for code validation.
"""
import os
import platform
from functools import lru_cache
import pandas as pd
import numpy as np
import itertools


from .config import options

# Import utility functions explicitly
from ._utils.check_data import (check_X, check_X_unique_samples, check_labels, check_match_X_labels,
                                check_array_like, check_superset_subset,
                                check_col_in_df)
from ._utils.check_models import check_mode_class, check_model_kwargs
from ._utils.check_type import (check_number_range, check_number_val, check_str, check_bool,
                                check_dict, check_tuple, check_list_like,
                                check_ax)

from ._utils.new_types import ArrayLike1D, ArrayLike2D

from ._utils.decorators import (catch_runtime_warnings, CatchRuntimeWarnings,
                                catch_convergence_warning, ClusteringConvergenceException,
                                catch_invalid_divide_warning,
                                doc_params)

from ._utils.utils_output import (print_out, print_start_progress, print_progress, print_finished_progress)


# Folder structure
def _folder_path(super_folder, folder_name):
    """Modification of separator (OS-depending)"""
    path = os.path.join(super_folder, folder_name + SEP)
    return path


SEP = "\\" if platform.system() == "Windows" else "/"
FOLDER_PROJECT = os.path.dirname(os.path.abspath(__file__))
FOLDER_DATA = _folder_path(FOLDER_PROJECT, '_data')
URL_DATA = "https://github.com/breimanntools/xomics/tree/master/xomics/data/"


# I Constants
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


# II Helper functions
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


# III MAIN FUNCTIONS
# Caching for data loading for better performance (data loaded ones)
@lru_cache(maxsize=None)
def read_excel_cached(name, index_col=None):
    """Load cached dataframe to save loading time"""
    df = pd.read_excel(name, index_col=index_col)
    return df.copy()


@lru_cache(maxsize=None)
def read_csv_cached(name, sep=None):
    """Load cached dataframe to save loading time"""
    df = pd.read_csv(name, sep=sep)
    return df.copy()


# Main check functions
def check_verbose(verbose):
    if verbose is None:
        # System level verbosity
        verbose = options['verbose']
    else:
        check_bool(name="verbose", val=verbose)
    return verbose

