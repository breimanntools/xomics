"""
This is a script for ...
"""
import time
import pandas as pd
import numpy as np
import matplotlib as mpl
from adjustText import adjust_text
from matplotlib import pyplot as plt
import xomics.utils as ut


# TODO into utils and remove
def set_legend_handles_labels(ax=None, dict_color=None, list_cat=None, labels=None,
                              y=-0.2, x=0.5, ncol=3, fontsize=11, weight="normal",
                              lw=0, edgecolor=None, return_handles=False, loc=9,
                              labelspacing=0.2, columnspacing=2, title=None, fontsize_legend=None,
                              title_align_left=False, fontsize_weight="normal"):
    """Get legend handles from dict_color"""
    if list_cat is None:
        raise ValueError("'list_cat' should not be None")
    dict_leg = {cat: dict_color[cat] for cat in list_cat}
    if edgecolor is None:
        f = lambda l, c: mpl.patches.Patch(facecolor=l, label=c, lw=lw, edgecolor=l)
    else:
        f = lambda l, c: mpl.patches.Patch(facecolor=l, label=c, lw=lw, edgecolor=edgecolor)
    handles = [f(l, c) for c, l in dict_leg.items()]
    if return_handles:
        return handles, labels
    if labels is None:
        labels = list(dict_leg.keys())
    args = dict(prop={"weight": weight, "size": fontsize})
    if fontsize_legend is not None:
        args["title_fontproperties"] = {"weight": fontsize_weight, "size": fontsize_legend}
    legend = ax.legend(handles=handles, labels=labels, bbox_to_anchor=(x, y), ncol=ncol,
                       loc=loc, labelspacing=labelspacing, columnspacing=columnspacing,
                       borderpad=0, **args, title=title)
    if title_align_left:
        legend._legend_box.align = "left"
    return ax


def _adjust_text(df_show=None, size=12, n=15, n_objects=200, set_bbox=False,
                 x=None, y=None, z=None, weight="bold",
                 force_text=0.2, force_points=0.2, force_object=0.2, precision=0.01):
    """"""
    texts = []
    i = 0
    for _, row in df_show.iterrows():
        x_val = row[x]
        y_val = row[y]
        fontdict = dict(size=size, color="black", weight=weight)
        label = row[z] if i < n else "   "
        if i <= n and set_bbox:
            props = dict(boxstyle='round', facecolor='lightgray', alpha=0.75)
            texts.append(plt.text(x_val, y_val, label, fontdict=fontdict, bbox=props))
        else:
            texts.append(plt.text(x_val, y_val, label, fontdict=fontdict))
        i += 1
    # Adjust text
    adjust_text(texts[0:n-1],
                add_objects=texts[0:n_objects],
                arrowprops=dict(arrowstyle="-", color='black', lw=0.5),
                ha="right", va="bottom",
                force_points=force_points,
                force_text=force_text,
                force_objects=force_object,
                precision=precision,
                save_format="png")