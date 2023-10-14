"""
This is a script for ...
"""
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

import xomics.utils as ut


# I Helper Functions
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


# II Main Functions
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
