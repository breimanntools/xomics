"""
Config with folder structure
"""

import os
import platform
from pathlib import Path


# Helper Function
def _folder_path(super_folder, folder_name):
    """Modification of separator (OS depending)"""
    path = os.path.join(super_folder, folder_name + SEP)
    return path


# Folder
SEP = "\\" if platform.system() == "Windows" else "/"
FOLDER_PROJECT = str(Path(__file__).parent).replace('/', SEP) + SEP
FOLDER_DATA = _folder_path(FOLDER_PROJECT, 'data')
FOLDER_RESULTS = _folder_path(FOLDER_PROJECT, 'results')

# Plotting settings
DPI = 300
FIG_FORMAT = "pdf"
ARGS_SAVE = dict(dpi=DPI, bbox_inches="tight")
LEGEND_FONTSIZE = 16

