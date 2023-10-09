"""
This is a script for ...
"""
import time
import pandas as pd
import numpy as np
import matplotlib as mpl
from adjustText import adjust_text
from matplotlib import pyplot as plt
import xomics._utils as ut


# I Helper functions
def _determine_color(x, y, th_filter):
    if abs(x) < th_filter[1] or y < th_filter[0]:
        return ut.COLOR_NOT_SIG
    return ut.COLOR_DOWN if x < 0 else ut.COLOR_UP


def _set_text(x, y, label, fontdict, props, box):
    color_key = "color"
    if box:
        fontdict[color_key] = "white"
        text = plt.text(x, y, label, fontdict=fontdict, bbox=props)
    else:
        fontdict[color_key] = "black"
        text = plt.text(x, y, label, fontdict=fontdict)
    return text


# II Main functions
# TODO check if necessary and simplify
def color_filter(df=None, th_p=2.0, th_ratio=0.5, col_ratio=None, col_pval=None, gene_list=None):
    """Classify significant values by color"""
    colors = []
    for i, row in df.iterrows():
        gene = row["Gene_Name"]
        ratio = row[col_ratio]
        p_val = row[col_pval]
        if gene_list is not None and gene in gene_list:
            if ratio > 0:
                colors.append(ut.COLOR_UP)
            else:
                colors.append(ut.COLOR_DOWN)
        elif p_val >= th_p:
            if ratio >= th_ratio:
                colors.append(ut.COLOR_UP)
            elif ratio <= -th_ratio:
                colors.append(ut.COLOR_DOWN)
            else:
                colors.append(ut.COLOR_NOT_SIG)
        else:
            colors.append(ut.COLOR_NOT_SIG)
    return colors


def label_filter(df=None, gene_list=None, th_p_text=2.0, th_neg_ratio=-0.5, th_pos_ratio=0.5,
                 col_ratio=None, col_pval=None, avoid_conflict=0.5):
    """Classify significant values by color"""
    # Adjust avoid factor
    if avoid_conflict > 1:
        avoid_conflict = avoid_conflict/100
    ac = 1 - avoid_conflict
    # Filter labels
    labels = []
    for i, row in df.iterrows():
        gene = row["Gene_Name"]
        ratio = row[col_ratio]
        p_val = row[col_pval]
        x, y = ratio, p_val
        if (gene_list is not None and gene in gene_list) or (p_val > th_p_text and (ratio < th_neg_ratio or ratio > th_pos_ratio)):
            label = (gene, x, y)
            labels.append(label)
        # Save significant but not to be shown genes to avoid overlapping conflict
        elif p_val > th_p_text*0.75 and (ratio < th_neg_ratio * ac or ratio > th_pos_ratio * ac):
            label = ("", x, y)
            labels.append(label)
    return labels


def set_labels(labels=None, objects=None, fig_format="png", label_size=8, precision=0.01, label_bold=False,
               force_points=0.75, force_text=0.75, force_objects=0.25, box=False, alpha=0.85, verbose=True,
               th_filter=None):
    """Automatically set labels for a plot.

    Parameters:
    - labels (list of tuple): Each tuple contains (label, x, y)
    - objects (list): Objects to consider for label placement
    - fig_format (str): Format of the output plot. Default is 'png'
    - label_size (int): Font size of labels. Default is 8
    - precision (float): Precision for text placement. Default is 0.01
    - label_bold (bool): If True, makes label text bold. Default is False
    - force_points (float): Repel force from points. Default is 0.75
    - force_text (float): Repel force from text. Default is 0.75
    - force_objects (float): Repel force from objects. Default is 0.25
    - box (bool): If True, puts a box around text. Default is False
    - alpha (float): Opacity of box. Default is 0.85
    - verbose (bool): If True, prints verbose output. Default is True
    - th_filter (tuple): Threshold for filtering labels based on x and y values

    Returns:
    None
    """
    fontdict = dict(size=label_size)
    if label_bold:
        fontdict.update(weight="bold")
    props = dict(boxstyle='round', alpha=alpha, edgecolor="white")
    texts = []
    for label, x, y in labels:
        color = _determine_color(x, y, th_filter)
        props.update(dict(facecolor=color))
        if label is not None:
            texts.append(_set_text(x, y, label, fontdict=fontdict, props=props, box=box))
    if verbose:
        print("{} elements have to be iteratively placed".format(len(texts)))
    adjust_text(texts,
                add_objects=objects,
                precision=precision,
                force_points=force_points,
                force_text=force_text,
                force_objects=force_objects,
                arrowprops=dict(arrowstyle="-", color='gray', lw=0.5),
                save_format=fig_format)


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