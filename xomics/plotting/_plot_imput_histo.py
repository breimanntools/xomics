"""
This is a script for ...
"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd

import xomics
import xomics.utils as ut


# I Helper Functions


# II Main Functions
def plot_imput_histo(df_raw=None, df_imp=None, cols_quant=None, d_min=None, up_mnar=None,
                     alpha=0.75, binwidth=0.4, colors=None, y_max=None, x_max=None, **kwargs):
    """"""
    colors = xomics.plot_get_clist(n_colors=3) if colors is None else colors
    plt.figure(figsize=(6, 5))
    _args = dict(binwidth=binwidth, **kwargs)
    vals_imp = df_imp[cols_quant].values.flatten()
    vals_raw = df_raw[cols_quant].values.flatten()
    x_max = max(plt.xlim()) if x_max is None else x_max
    x_min = min(plt.xlim())
    y_max = max(plt.ylim()) if y_max is None else y_max
    y_min = min(plt.ylim())
    # Drop missing values (NaNs) from vals_imp and vals_raw
    vals_imp = vals_imp[~np.isnan(vals_imp)]
    vals_raw = vals_raw[~np.isnan(vals_raw)]
    sns.histplot(data=vals_imp, color=colors[1], alpha=alpha, **_args)
    sns.histplot(data=vals_raw, color=colors[0], alpha=1, **_args)
    plt.xlabel("log2 LFQ")
    plt.legend(labels=['Imputed', 'Raw'])
    plt.ylim((y_min, y_max))
    plt.xlim((x_min, x_max))
    # Add annotation text
    plt.axvline(d_min, color='black', ls="--")
    plt.text(d_min * 1.01, y_max*0.95, "Dmin")
    plt.axvline(up_mnar, color='black', ls="--")
    plt.text(up_mnar * 1.01, y_max*0.95, "upMNAR")
    plt.text(x_max*0.65, y_max*0.6, f"n Raw: {len(vals_raw)}\n"
                                   f"n Imputed: {len(vals_imp)}")

    sns.despine()

