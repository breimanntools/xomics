"""
This is a script for comparing raw vs imputed data by histogram.
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from typing import Optional, Union, Tuple, List

import xomics as xo
import xomics.utils as ut


# I Helper Functions


# II Main Functions
# TODO add check function, improve interface and make consistent, add tests, add tutorial

def plot_imput_histo(ax: Optional[plt.Axes] = None,
                     figsize: Tuple[int, int] = (6, 5),
                     df_raw: pd.DataFrame = None,
                     df_imp: pd.DataFrame = None,
                     cols_quant: List[str] = None,
                     d_min: Optional[float] = None,
                     up_mnar: Optional[float] = None,
                     alpha: float = 0.75,
                     binwidth: float = 0.4,
                     colors: Optional[List[str]] = None,
                     y_max: Optional[float] = None,
                     x_max: Optional[float] = None,
                     **kwargs):
    """
    Plot histogram of raw and imputed data

    Parameters
    ----------
    ax
        The axes upon which to draw the plot. If None, a new figure and axes are created.
    figsize
        The size of the figure to create.
    df_raw
        Dataframe containing the raw data.
    df_imp
        Dataframe containing the imputed data.
    cols_quant
        Columns to consider for the quantitative analysis.
    d_min
        The minimum value of the data for annotation.
    up_mnar
        Upper value for MNAR annotation.
    alpha
        Transparency level of the histogram bars for imputed data.
    binwidth
        Width of the histogram bins.
    colors
        List of colors for the histogram bars. If None, a default set of colors is used.
    y_max
        The maximum value for the y-axis.
    x_max
        The maximum value for the x-axis.
    **kwargs
        Additional keyword arguments passed to seaborn's histplot.

    Returns
    -------
    ax
        Axes object.
    """
    # Check input

    # Pre-process data
    colors = xo.plot_get_clist(n_colors=3) if colors is None else colors
    _args = dict(binwidth=binwidth, **kwargs)
    vals_imp = df_imp[cols_quant].values.flatten()
    vals_raw = df_raw[cols_quant].values.flatten()
    # Drop missing values (NaNs) from vals_imp and vals_raw
    vals_imp = vals_imp[~np.isnan(vals_imp)]
    vals_raw = vals_raw[~np.isnan(vals_raw)]

    # Plot
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.histplot(ax=ax, data=vals_imp, color=colors[1], alpha=alpha, **_args)
    sns.histplot(ax=ax, data=vals_raw, color=colors[0], alpha=1, **_args)

    # Get axis limits
    x_max = max(ax.get_xlim()) if x_max is None else x_max
    x_min = min(ax.get_xlim())
    y_max = max(ax.get_ylim()) * 1.05 if y_max is None else y_max
    y_min = min(ax.get_ylim())

    ax.set_xlabel("log2 LFQ")
    ax.set_ylim((y_min, y_max))
    ax.set_xlim((x_min, x_max))

    ax.legend(labels=['Imputed', 'Raw'], fontsize=ut.plot_gco(), labelspacing=0.2)

    # Add annotation text
    lw = ut.plot_gco(option="axes.linewidth")
    args = dict(size=ut.plot_gco()-1)
    ax.axvline(d_min, color='black', ls="--", lw=lw)
    ax.text(d_min * 1.01, y_max * 0.95, "Dmin", **args)
    ax.axvline(up_mnar, color='black', ls="--", lw=lw)
    ax.text(up_mnar * 1.01, y_max * 0.95, "upMNAR", **args)

    str_n = f"n Imputed: {len(vals_imp)}\nn Raw: {len(vals_raw)}"
    ax.text(x_max, y_max * 0.6, str_n, **args, ha="right")

    # Adjust plot
    sns.despine()
    plt.tight_layout()
    return ax
