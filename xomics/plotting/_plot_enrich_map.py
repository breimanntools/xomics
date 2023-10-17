"""
This is a script for ...
"""
from matplotlib import pyplot as plt
import seaborn as sns

import xomics.utils as ut
from xomics.plotting._utils_plot import set_legend_handles_labels


# I Helper Functions


# II Main Functions
def plot_enrich_map(df=None, row_colors=None, col_colors=None,
                    method="complete", figsize=(8, 7), wide=False,
                    font_scale=0.8, x_legend=1.2):
    """"""
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
                              x=x_legend, ncol=1, list_cat=dict_color.keys())


def plot_enrich_rank():
    """Plot enrichment by fold enrichment or significance with gene set size"""

