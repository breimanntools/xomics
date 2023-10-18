"""
This is a script for ...
"""
import time
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Settings


# I Helper Functions
def _plot_scatter(ax=None, df_plot=None, df_imputed=None, cols=None, legend=False, title=None, hue=None):
    """"""
    df = df_plot[cols].T.describe().T
    df.index = df_imputed.index
    df[hue] = df_imputed[hue]
    ax = sns.scatterplot(ax=ax, data=df, x="std", y="mean",
                         hue=hue, siza="count", palette="RdBu", alpha=0.75, legend=legend)
    ax.set_xlim(0, 8)
    ax.set_ylim(10, 35)
    ax.set_title(title)
    return ax


# II Main Functions
def plot_scatter(df_raw=None, df_imputed=None, cols=None, group=None):
    """"""
    args = dict(df_imputed=df_imputed.copy(), cols=cols, hue=f"CS_{group}")
    plt.figure(figsize=(6, 7))
    fig, ax = plt.subplots(1, 2)
    _plot_scatter(ax=ax[0], df_plot=df_raw.copy(), title="Raw", **args)
    sns.despine()
    _plot_scatter(ax=ax[1], df_plot=df_imputed.copy(), title="Imputed", **args, legend=True)
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.suptitle(f"log2 LFQ group {group}", weight="bold", y=1)

