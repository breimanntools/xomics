"""
This is a script for ...
"""
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import xomics as xo
import dev_scripts._utils as ut

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions


# II Main Functions
def filter_usage():
    """"""
    df = pd.read_excel(ut.FOLDER_DATA + "raw_data_proteomics_lfq.xlsx") #.head(50)
    groups = ["A", "B", "C", "D", "E"]
    pp = xo.PreProcess(str_quant="log2 LFQ")
    #df = pp.filter_nan(df=df, groups=groups)
    df = pp.filter_groups(df=df, groups=groups, min_pct=0)
    print(df)


def imputation_usage():
    """"""
    df = pd.read_excel(ut.FOLDER_DATA + "raw_data_proteomics_lfq.xlsx") #.head(50)
    groups = ["A", "B", "C", "D", "E"]
    #d_min, up_mnar, d_max = xo.cImpute().get_limits(df=df, groups=groups)
    #df_imp = xo.cImpute().run(df=df, groups=groups, min_cs=0.5)

    pp = xo.PreProcess()
    cols_quant = pp.get_qcols(df=df, groups=groups)
    dict_qcol_group = pp.get_dict_qcol_group(df=df, groups=groups)
    dict_group_qcols = pp.get_dict_group_qcols(df=df, groups=groups)
    #df = pp.filter(df=df)
    df_fc = pp.run(df=df, groups=groups, groups_ctrl=["A"])
    df_fc = pp.add_ids(df=df_fc, list_ids=df["Protein IDs"])
    xo.plot_settings(weight_bold=True, short_ticks=True)
    names_to_label = [f"Protein000{i}" for i in range(1, 9)]

    ax = xo.plot_volcano(df=df_fc,
                         col_fc="log2 FC (C/A)",
                         col_pval="-log10 p-value (C/A)",
                         colors="gray",
                         alpha=0.5,
                         legend=False,
                         edge_color="black")
    print(df_fc)
    plt.xlim(-8, 8)
    plt.show()
    df_fc = df_fc[df_fc["Protein IDs"].isin(names_to_label)]
    xo.plot_volcano(ax=ax,
                    df=df_fc,
                    #names_to_annotate=names_to_label,
                    colors="tab:red",
                    #col_names="Protein IDs",
                    col_fc="log2 FC (C/A)",
                    col_pval="-log10 p-value (C/A)",
                    legend=False)
    plt.xlim(-8, 8)
    plt.show()
    """
    xo.plot_settings(font_scale=1, weight_bold=True)
    xo.plot_imput_histo(df_raw=df, df_imp=df_imp, cols_quant=cols_quant, d_min=d_min, up_mnar=up_mnar)
    plt.show()
    cs_scores = [1, 0.75, 0.5, 0]
    xo.plot_settings(font_scale=1, weight_bold=False)
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    fig.subplots_adjust(hspace=0.5)
    for i, cs_score in enumerate(cs_scores):
        # Calculate the subplot's position (row and column)
        row = i // 2
        col = i % 2
        # Create a subplot for the current CS score
        ax = axs[row, col]
        df_imp = xo.cImpute().run(df=df, groups=groups, min_cs=cs_score)
        xo.plot_imput_histo(ax=ax, df_raw=df, df_imp=df_imp, cols_quant=cols_quant,
                            d_min=d_min, up_mnar=up_mnar, y_max=6700, x_max=38)
        ax.set_title(f"CS: {cs_score}")
    plt.tight_layout()
    plt.show()
    """

# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    filter_usage()
    #imputation_usage()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
