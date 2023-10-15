"""
This is a script for comparing raw vs imputed data by histogram.
"""
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

import xomics as xo
import xomics.utils as ut


# I Helper Functions


# II Main Functions
def plot_imput_histo(ax=None, figsize=(6, 5), df_raw=None, df_imp=None, cols_quant=None, d_min=None, up_mnar=None,
                     alpha=0.75, binwidth=0.4, colors=None, y_max=None, x_max=None, **kwargs):
    """Plot histogram of raw and imputed data"""
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
    if ax is None:
        return fig, ax
    else:
        return ax
