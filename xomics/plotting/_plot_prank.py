"""
This is a script for ...
"""
from matplotlib import pyplot as plt
import seaborn as sns

import xomics.utils as ut
from xomics.plotting._utils_plot import set_legend_handles_labels, _adjust_text


# I Helper Functions
# TODO optimize and unify plotting functions
# TODO highlight hidden gem in (P-score > 0.75, E-score <= 0.1 (not pathway))
# TODO provide simple highlighted functions for plotting

# II Main Functions
def plot_prank(df=None, col_p_sore=None, col_e_score=None, col_names=None,
               show_names=False, force_text=0.2, force_points=0.2, force_object=0.2):
    """"""

    plt.figure(figsize=(5, 5))
    args = dict(size=2, x=ut.COL_P_SCORE, y=ut.COL_E_SCORE,
                legend=False, edgecolor="black", clip_on=False)
    sns.scatterplot(data=df, color="gray", **args)
    df = df[df[ut.COL_PE_MEAN] >= 0.5]
    if show_names:
        kwargs = dict(force_text=0.2, force_points=0.2, force_object=0.2)

        _adjust_text(df_show=df, x=ut.COL_P_SCORE,
                     y=ut.COL_E_SCORE, z="Gene names", size=7, n=50, weight="normal", **kwargs)
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


def plot_prank_scatter(df=None, col_p_sore=None, col_e_score=None, col_names=None,
                       show_names=False, force_text=0.2, force_points=0.2, force_object=0.2):
    """"""
    plt.figure(figsize=(5, 5))
    args = dict(size=2, x=col_p_sore, y=col_e_score,
                legend=False, edgecolor="black", clip_on=False)
    sns.scatterplot(data=df, color="gray", **args)
    df = df[df[ut.COL_PE_MEAN] >= 0.5]
    if show_names:
        _adjust_text(df_show=df, x=col_p_sore,
                     y=col_e_score, z=col_names, size=7, n=50, weight="normal",
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
    plt.ylabel(f(col_e_score))
    plt.xlabel(f(col_e_score))
    dict_color = {f"{f(ut.COL_PE_MEAN)}>=0.6": "red",
                  f"{f(ut.COL_PE_MEAN)}>=0.5": "orange"}
    set_legend_handles_labels(ax=plt.gca(), dict_color=dict_color, list_cat=dict_color.keys(),
                              ncol=1, fontsize=10, columnspacing=1, labelspacing=0.05,
                              y=1, x=0.8)