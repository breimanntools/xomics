"""
This is a script for plotting functions
"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

import xomics._utils as ut
from xomics.plotting._plotting import color_filter, label_filter, set_labels, set_legend_handles_labels, _adjust_text

# TODO optimize and unify plotting functions
# TODO highlight hidden gem in (P-score > 0.75, E-score <= 0.1 (not pathway))
# TODO provide simple highlighted functions for plotting


# I Helper Functions
def _check_col(df, col=None):
    """Check if col in df"""
    if col not in list(df):
        raise ValueError("{} should be one of following: {}".format(col, list(df)))


def _check_gene_list(df=None, gene_list=None):
    """Check if gene in df"""
    if gene_list is not None:
        all_genes = df["Gene_Name"].tolist()
        not_in_all_genes = [x for x in gene_list if x not in all_genes]
        if len(not_in_all_genes):
            raise ValueError("Genes from 'gene_list' not in 'df_ratio_pval': {}".format(not_in_all_genes))


def _check_gene_values(gene=None, x=None, y=None):
    """Check if x or y is missing value"""
    if str(x) == "nan" or str(y) == "nan":
        raise ValueError("Missing values for gene '{}'".format(gene))


# II Main Functions
# Volcano plot
def plot_volcano(df=None, col_fc=None, col_pval=None, gene_list=None,
                 th_filter=(0.05, 0.5), th_text=None,
                 precision=0.01, force=(0.5, 0.5, 0.25), avoid_conflict=0.25,
                 fig_format="png", verbose=False, loc_legnd=2,
                 filled_circle=True, box=True, label_bold=False, label_size=8, minor_ticks=True):
    """
     Generate and display a volcano plot based on fold-change and p-value _data.

     Parameters
     ----------
     df : DataFrame
         DataFrame containing fold-change and p-values.
     col_fc : str
         Column name containing fold change values.
     col_pval : str
         Column name containing p-values.
     gene_list : list, optional
         List of specific genes to label on the plot.
     th_filter : tuple, default=(0.05, 0.5)
         Thresholds for p-value and fold-change for filtering.
     th_text : tuple, optional
         Thresholds for labeling; format is (p_val, lower_ratio, upper_ratio).
     precision : float, default=0.01
         Precision level for text placement.
     force : tuple, default=(0.5, 0.5, 0.25)
         Repulsion forces to adjust text (points, text, objects).
     avoid_conflict : float, default=0.25
         Level of label overlap to avoid.
     fig_format : str, default="png"
         File format for saving the plot.
     verbose : bool, default=False
         If True, print additional information.
     loc_legnd : int, default=2
         Location index for the plot legend.
     filled_circle : bool, default=True
         If True, circles in the plot are filled.
     box : bool, default=True
         If True, labels will have a bounding box.
     label_bold : bool, default=False
         If True, labels will be bold.
     label_size : int, default=8
         Font size for labels.
     minor_ticks : bool, default=True
         If True, shows minor ticks in plot.

     Returns
     -------
     matplotlib.Axes
         The Axes object representing the plot.
     """
    # TODO Include labeling for selected genes (e.g., for GO terms), add FDR correction
    # Initial parameter validation (your _check_col function and _check_gene_list function)
    _check_col(df, col=col_fc)
    _check_col(df, col=col_pval)
    if gene_list:
        _check_gene_list(df=df, gene_list=gene_list)

    # Unpack thresholds
    th_p, th_ratio = th_filter
    if th_text is None:
        th_text = (th_p, -th_ratio, th_ratio)
    th_p_text, th_neg_ratio, th_pos_ratio = th_text

    # Rescale p-value if needed
    if th_p < 0.5:
        th_p = -np.log10(th_p)
    if th_p_text < 0.5:
        th_p_text = -np.log10(th_p_text)

    # Pre-processing for labels and colors
    kwargs = {}
    labels = []
    if gene_list:
        kwargs_filter = dict(df=df, col_ratio=col_fc, col_pval=col_pval, gene_list=gene_list)
        colors = color_filter(**kwargs_filter, th_p=th_p, th_ratio=th_ratio)
        labels = label_filter(**kwargs_filter, th_p_text=th_p_text, th_neg_ratio=th_neg_ratio,
                              th_pos_ratio=th_pos_ratio, avoid_conflict=avoid_conflict)
        kwargs['color'] = colors
        if not filled_circle:
            kwargs['edgecolor'] = colors
            kwargs['facecolor'] = 'none'

    # Plot settings
    df.plot(x=col_fc, y=col_pval, kind="scatter", figsize=(5, 5), **kwargs)

    # Minor ticks
    if minor_ticks:
        plt.minorticks_on()

    # Add threshold lines and set plot limits
    plt.axhline(y=th_p, linestyle='--', color='grey', linewidth=1.5)
    plt.axvline(x=th_ratio, linestyle='--', color='grey', linewidth=1.5)
    plt.axvline(x=-th_ratio, linestyle='--', color='grey', linewidth=1.5)
    plt.xlim([df[col_fc].min() - 1, df[col_fc].max() + 1])
    plt.ylim([0, df[col_pval].max() * 1.1])

    # If gene_list is provided, set labels
    if gene_list:
        objects = [plt.gca().lines[-3], plt.gca().lines[-2], plt.gca().lines[-1]]
        set_labels(labels,
                   objects=objects,
                   th_filter=th_filter,
                   fig_format=fig_format,
                   force_points=force[0],
                   force_text=force[1],
                   force_objects=force[2],
                   box=box,
                   verbose=verbose,
                   precision=precision,
                   label_size=label_size,
                   label_bold=label_bold)

    # Legend and Labels
    plt.legend(handles=[
        mpatches.Patch(color=ut.COLOR_UP, label='Up'),
        mpatches.Patch(color=ut.COLOR_DOWN, label='Down'),
        mpatches.Patch(color=ut.COLOR_NOT_SIG, label='Not Sig')
    ], loc=loc_legnd, frameon=False)

    plt.xlabel(col_fc, weight="bold")
    plt.ylabel(col_pval, weight="bold")
    return plt.gca()


# Enrichment plot
def plot_enrich_rank():
    """Plot enrichment by fold enrichment or significance with gene set size"""


def plot_enrich_map(df=None, row_colors=None, col_colors=None,
                      method="complete", figsize=(8, 7), wide=False,
                      font_scale=0.8, x_legend=1.2):
    """"""
    # TODO add to plotting functions
    ut.plot_settings(weight_bold=False, font_scale=font_scale,
                                no_ticks_y=True, short_ticks_x=True)
    fg_ratio = figsize[1]/figsize[0]
    ratio = 0.15
    dendrogram_ratio = [ratio*fg_ratio, ratio]
    linewidth = 2
    tree_kws = dict(linewidth=linewidth)
    cg = sns.clustermap(df, figsize=figsize, cmap="GnBu", method=method,
                        cbar_kws=None, mask=(df==0), vmin=0, vmax=2,
                        row_colors=row_colors, col_colors=col_colors,
                        #dendrogram_ratio=dendrogram_ratio,
                        tree_kws=tree_kws, linewidth=1,
                        linecolor="black", yticklabels=1)

    ax = plt.gca()
    cg.cax.set_visible(False)
    tick_labels = cg.ax_heatmap.xaxis.get_majorticklabels()
    plt.setp(tick_labels, rotation=45, ha="right", va="center", rotation_mode="anchor")
    labelsize = 10 if wide else 13
    ax.tick_params(labelsize=labelsize)
    # TODO dif color
    dict_color = {"Associated": ut.COLOR_GEM, "Not associated": "white"}
    set_legend_handles_labels(ax=cg.ax_heatmap, dict_color=dict_color, lw=1, edgecolor="black", y=-0.01,
                              x=x_legend, ncol=1, list_cat=dict_color.keys(),
                              fontsize=ut.LEGEND_FONTSIZE)


# EP plot (ExPlano plot)
def plot_explano(df=None, show_names=False, force_text=0.2, force_points=0.2, force_object=0.2):
    """"""

    plt.figure(figsize=(5, 5))
    args = dict(size=2, x=ut.COL_P_SCORE, y=ut.COL_E_SCORE,
                legend=False, edgecolor="black", clip_on=False)
    sns.scatterplot(data=df, color="gray", **args)
    df = df[df[ut.COL_PE_MEAN] >= 0.5]
    if show_names:
        _adjust_text(df_show=df, x=ut.COL_P_SCORE,
                     y=ut.COL_E_SCORE, z="Gene names", size=7, n=50, weight="normal",
                     force_points=force_points,
                     force_text=force_text,
                     force_object=force_object)
    sns.scatterplot(data=df, color="orange", **args)
    df = df[df[ut.COL_PE_MEAN] >= 0.6]
    sns.scatterplot(data=df, color="red", **args)
    sns.despine()
    plt.tight_layout()
    plt.xlim(0, 1.0)
    plt.ylim(-0.02, 1.0)
    f = lambda _str: _str.replace("_", " ")
    plt.ylabel(f(ut.COL_E_SCORE))
    plt.xlabel(f(ut.COL_P_SCORE))
    dict_color = {f"{f(ut.COL_PE_MEAN)}>=0.6": "red",
                  f"{f(ut.COL_PE_MEAN)}>=0.5": "orange"}
    set_legend_handles_labels(ax=plt.gca(), dict_color=dict_color, list_cat=dict_color.keys(),
                              ncol=1, fontsize=10, columnspacing=1, labelspacing=0.05,
                              y=1, x=0.8)


# Inferno plot (integration of multiple _data sets using information from volcano and enrichment plot)
def plot_inferno():
    """Create Inferno plot for multiple integrated omics _data"""
    # TODO create after integration function (in integrate) is finished (Clari)


# cImpute quality check plots
def plot_imput_histo(df_raw=None, df_imp=None, cols=None, d_min=None, up_mnar=None,
                     alpha=0.5, binwidth=0.4, **kwargs):
    """"""
    plt.figure(figsize=(6, 5))
    _args = dict(alpha=alpha, binwidth=binwidth, **kwargs)
    sns.histplot(df_imp[cols].to_numpy().flatten(), **_args)
    sns.histplot(df_raw[cols].to_numpy().flatten(), color="tab:orange", **_args)
    plt.xlabel("log2 LFQ")
    plt.legend(labels=['Imputed', 'Raw'])
    plt.axvline(d_min, color='black', ls="--")
    plt.text(d_min * 1.01, plt.ylim()[1] * 0.95, "Dmin")
    plt.axvline(up_mnar, color='black', ls="--")
    plt.text(up_mnar * 1.01, plt.ylim()[1] * 0.95, "upMNAR")
    sns.despine()
    plt.tight_layout()


def _plot_scatter(ax=None, df_plot=None, df_imp=None, cols=None, legend=False, title=None, hue=None, alpha=0.75):
    """"""
    df = df_plot[cols].T.describe().T
    df.index = df_imp.index
    df[hue] = df_imp[hue]
    ax = sns.scatterplot(ax=ax, data=df, x="std", y="mean", hue=hue, size="count", palette="RdBu", alpha=alpha,
                         legend=legend)
    ax.set_xlim(0, 8)
    ax.set_ylim(10, 35)
    ax.set_title(title)
    return ax


def plot_imput_scatter(df_raw=None, df_imp=None, cols=None, group=None, alpha=0.75):
    """"""
    args = dict(df_imp=df_imp.copy(), cols=cols, hue=f"CS_{group}", alpha=alpha)
    plt.figure(figsize=(6, 7))
    fig, ax = plt.subplots(1, 2)
    _plot_scatter(ax=ax[0], df_plot=df_raw.copy(), title="Raw", **args)
    sns.despine()
    _plot_scatter(ax=ax[1], df_plot=df_imp.copy(), title="Imputed", **args, legend=True)
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.suptitle(f"log2 LFQ group {group}", weight="bold", y=1)
