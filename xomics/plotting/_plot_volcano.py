"""
This is a script for ...
"""
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


import xomics.utils as ut
from .utils_plotting import (color_filter, label_filter, set_labels)


# I Helper Functions
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
def plot_volcano(df=None, col_fc=None, col_pval=None, gene_list=None,
                 th_filter=(0.05, 0.5), th_text=None,
                 precision=0.01, force=(0.5, 0.5, 0.25), avoid_conflict=0.25,
                 fig_format="png", verbose=False, loc_legnd=2,
                 filled_circle=True, box=True, label_bold=False, label_size=8, minor_ticks=True):
    """
     Generate and display a volcano plot based on fold-change and p-value data.

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
    ut.check_col_in_df(name_df="df", df=df, cols=[col_fc, col_pval])
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
