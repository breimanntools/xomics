"""
This is a script for ...
"""
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Optional, Tuple, Union, List
from adjustText import adjust_text

import xomics as xo
import xomics.utils as ut

# Constats
COL_SIG_SIZE = "sig_size"

# ---------------------------------
"""
# TODO make correct (wrong FDR formula
def plot_fdr_curve(df: pd.DataFrame, col_fc: str, col_pval: str, c=1.5, s0=1.0):
    Hyperbolic FDR Curve Plotting Function

    Purpose:
    To plot hyperbolic FDR curves, primarily used in mass spectrometry and use_cases studies for multiple hypothesis
     testing correction. These curves represent a balance between significance (p-value) and effect size (fold change).

    Background:
    The concept of the hyperbolic threshold originates from the understanding that we may tolerate weaker p-values
    (less statistical significance) when the effect size (fold change) is large, and vice-versa. Thus, the hyperbolic
    curves represent a combination of significance and effect size, controlling the overall false discovery rate (FDR).

    Relevant Paper:
    - "Perseus: A Bioinformatics Platform for Integrative Analysis of Proteomics Data in Cancer Research" Tyanova et al.
    - Fudge factor: https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/pmic.201600132

    Challenges:
    Translate FDR curve code from Perseus (not cleare where described) or this R code:
    https://github.com/Scavetta/Volcanto_Plots_FDR..


    # Define the hyperbolic curve function
    def fdr_curve(x, c, s0):
        return c * (np.abs(x) - s0 + np.sqrt(x ** 2 + s0 ** 2))
    x = np.linspace(-np.abs(df[col_fc].max()), np.abs(df[col_fc].max()), 1000)
    y_curve = fdr_curve(x, c, s0)
    x = np.linspace(-np.abs(df[col_fc].max()), np.abs(df[col_fc].max()), 1000)
    # Plotting the curves
    plt.plot(x, y_curve, label='Hyperbolic FDR curve', linestyle='--', color='gray')
    # TODO does not work
"""
# ---------------------------------


# I Helper Functions
def check_match_df_names(df=None, list_names=None, col_names=None):
    """Check if gene in df"""
    if col_names is None:
        raise ValueError(f"'col_names' should not be None.")
    all_names = df[col_names].to_list()
    if list_names is None:
        return all_names
    else:
        wrong_names = [x for x in list_names if x not in all_names]
        if len(wrong_names):
            raise ValueError(f"Following names from 'list_names' are not in 'col_names': {wrong_names}")
    return list_names


# II Main Functions
def plot_volcano(ax: Optional[plt.Axes] = None,
                 figsize: Tuple[int, int] = (6, 6),
                 df: pd.DataFrame = None,
                 col_fc: str = None,
                 col_pval: str = None,
                 col_names: Optional[str] = None,
                 col_cbar: Optional[str] = None,
                 th_fc: float = 0.5,
                 th_pval: float = 0.05,
                 names_to_annotate: Optional[list] = None,
                 colors: Optional[Union[str, List[str]]] = None,
                 colors_pos_neg_non: Optional[Tuple[str, str, str]] = None,
                 size: int = 50,
                 sizes_pos_neg_non: Optional[Tuple[int, int, int]] = None,
                 alpha: float = 1.0,
                 edge_color: str = "white",
                 edge_width: float = 0.5,
                 cmap: Optional[str] = 'viridis',
                 label_fontdict: Optional[dict] = None,
                 label_adjust_text_dict: Optional[dict] = None,
                 loc_legend: int = 2,
                 legend: bool = True,
                 minor_ticks: bool = False,
                 ) -> plt.Axes:
    """
    Generate and display a volcano plot based on fold-change and p-value data.

    Parameters
    ----------
    ax
        Axes object for the plot. If not provided, a new Axes object is created.
    figsize
        Size of the figure.
    df
        DataFrame containing fold-change and p-values.
    col_fc
        Column name containing fold change values.
    col_pval
        Column name containing p-values.
    col_names
        Columns with protein/gene names.
    th_fc
        Threshold for fold-change, applied for negative and positive values.
    th_pval
        Threshold for p-value, -log10 transformed before applied.
    names_to_annotate
        List of specific protein/gene names to label on the plot.
    colors
        Colors to use for the plot. Either a single color or a list of colors.
    colors_pos_neg_non
        Tuple containing colors for significant positive, significant negative, and non-significant points.
    size
        Size of the points.
    sizes_pos_neg_non
        Tuple containing sizes for significant positive, significant negative, and non-significant points.
    alpha
        Alpha transparency level for points.
    edge_color
        Color of the edge of the points.
    edge_width
        Width of the edge of the points.
    label_fontdict
        Dictionary of font properties for labels.
    label_adjust_text_dict
        Dictionary of properties for adjust_text function to adjust overlapping labels.
    label_arrow
        If True, black arrows are used for annotations (can be adjusted using ``label_adjust_text_dict``)
    loc_legend
        Location index for the plot legend.
    legend
        If True, display the legend. If False, hide the legend.
    minor_ticks
        If True, shows minor ticks in plot.

    Returns
    -------
    ax
        The Axes object representing the plot.

    See Also
    --------
    - Adjust text package: `Adjust text <https://adjusttext.readthedocs.io/en/latest/>`_
    """
    # Initial parameter validation
    ut.check_ax(ax=ax, accept_none=True)
    ut.check_tuple(name="figsize", val=figsize, n=2, accept_none=True)
    df = ut.check_df(name="df", df=df, cols_req=[col_fc, col_pval])
    if col_names is not None or names_to_annotate is not None:
        ut.check_col_in_df(name_df="df", df=df, cols=col_names, name_cols="col_names")
        names_to_annotate = ut.check_list_like(name="names_to_annotate", val=names_to_annotate, accept_none=False)
        names_to_annotate = check_match_df_names(df=df, list_names=names_to_annotate, col_names=col_names)
    if col_cbar is not None:
        ut.check_col_in_df(name_df="df", df=df, cols=col_cbar, accept_nan=True)
    ut.check_number_range(name="th_fc", val=th_fc, min_val=0, just_int=False)
    ut.check_number_range(name="th_pval", val=th_pval, min_val=0, max_val=1, just_int=False)
    colors = ut.check_list_like(name="colors", val=colors, accept_none=True, accept_str=True)
    ut.check_tuple(name="colors_pos_neg_non", val=colors_pos_neg_non, accept_none=True, n=3, check_n=True)
    ut.check_number_range(name="size", val=size, min_val=1, just_int=True)
    ut.check_tuple(name="size_pos_neg_non", val=colors_pos_neg_non, accept_none=True, n=3, check_n=True)
    ut.check_number_range(name="alpha", val=alpha, min_val=0, max_val=1, just_int=False)
    ut.check_number_range(name="edge_width", val=edge_width, min_val=0, just_int=False)
    ut.check_dict(name="label_fontdict", val=label_fontdict, accept_none=True)
    ut.check_dict(name="label_adjust_text_dict", val=label_adjust_text_dict, accept_none=True)
    ut.check_bool(name="legend", val=legend)
    ut.check_bool(name="minor_ticks", val=minor_ticks)
    # Rescale p-value
    th_pval = -np.log10(th_pval)

    # Plot settings
    if colors_pos_neg_non is None:
        color_non_sig, color_sig_pos, color_sig_neg = xo.plot_get_clist(n_colors=3)
    else:
        color_sig_pos, color_sig_neg, color_non_sig = colors_pos_neg_non
    if sizes_pos_neg_non is None:
        size_sig_pos = size_sig_neg = size_non_sig = 1
        size_max = size
    else:
        size_sig_pos, size_sig_neg, size_non_sig = sizes_pos_neg_non
        size_max = max(sizes_pos_neg_non)/min(sizes_pos_neg_non)*size

    dict_color = {ut.STR_SIG_POS: color_sig_pos, ut.STR_SIG_NEG: color_sig_neg, ut.STR_NON_SIG: color_non_sig}
    dict_size = {ut.STR_SIG_POS: size_sig_pos, ut.STR_SIG_NEG: size_sig_neg, ut.STR_NON_SIG: size_non_sig}
    df[ut.COL_SIG_CLASS] = ut.get_sig_classes(df=df, col_fc=col_fc, col_pval=col_pval, th_pval=th_pval, th_fc=th_fc)
    df_plot = df.copy()
    df_plot[COL_SIG_SIZE] = [dict_size[c] for c in df[ut.COL_SIG_CLASS]]

    # Plotting
    if ax is None:
        plt.figure(figsize=figsize)
    kwargs = dict(data=df_plot,
                  x=col_fc,
                  y=col_pval,
                  edgecolor=edge_color,
                  linewidth=edge_width,
                  size=COL_SIG_SIZE,
                  sizes=(size, size_max),
                  alpha=alpha)
    if col_cbar is None:
        ax = sns.scatterplot(**kwargs,
                             hue=ut.COL_SIG_CLASS,
                             palette=dict_color)
    else:
        ax = sns.scatterplot(**kwargs,
                             hue=col_cbar,
                             palette=cmap)

        # Create mappable object for colorbar
        norm = plt.Normalize(df[col_cbar].min(), df[col_cbar].max())
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        # Add colorbar
        cbar = plt.colorbar(mappable, ax=ax, shrink=0.5)
        cbar.set_label(col_cbar)
        cbar.ax.spines['left'].set_visible(False)  # Hide the left spine
        cbar.outline.set_visible(False)
    plt.xlabel(col_fc)
    plt.ylabel(col_pval)
    sns.despine()

    # Adjust plot
    # Set customized colors
    if colors is not None:
        if len(colors) == 1:
            colors = colors * len(df)
        scatter = ax.collections[-1]  # Get the underlying scatter object
        scatter.set_facecolor(colors)

    # Add threshold lines and set plot limits
    lw = ut.plot_gco(option="axes.linewidth") * 0.8
    plt.axhline(y=th_pval, linestyle='--', color='black', linewidth=lw)
    plt.axvline(x=th_fc, linestyle='--', color='black', linewidth=lw)
    plt.axvline(x=-th_fc, linestyle='--', color='black', linewidth=lw)
    plt.xlim((df[col_fc].min() - 1, df[col_fc].max() + 1))
    plt.ylim((0, df[col_pval].max() * 1.1))

    # Minor ticks
    if minor_ticks:
        plt.minorticks_on()

    # Set annotation
    if names_to_annotate is not None:
        labels = [(row[col_names], row[col_fc], row[col_pval]) for i, row in df.iterrows()
                  if row[col_names] in names_to_annotate and not np.isnan(row[col_fc])]
        fontdict = dict(size=xo.plot_gcfs()-8)
        if label_fontdict is not None:
            fontdict.update(**label_fontdict)
        texts = [plt.text(x, y, label, fontdict=fontdict) for label, x, y in labels]
        label_adjust_text_dict = {} if label_adjust_text_dict is None else label_adjust_text_dict
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='black'), **label_adjust_text_dict)

    # Legend and Labels
    if not legend or col_cbar is not None:
        ax.legend().set_visible(False)
    else:
        xo.plot_legend(dict_color=dict_color,
                       list_cat=[ut.STR_SIG_NEG, ut.STR_SIG_POS, ut.STR_NON_SIG], ncol=1,
                       marker="o",
                       loc=loc_legend,
                       labelspacing=0.1, handletextpad=0.0)


    plt.tight_layout()
    return plt.gca()
