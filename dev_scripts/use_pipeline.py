"""
This is a script for ...
"""
import os
import platform
from pathlib import Path
import time
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import dev_scripts._utils as ut

import xomics as xo

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe

# TODO (!)
# TODO highlight hidden gem in (P-score > 0.75, E-score <= 0.1 (not pathway)) and PE-mean range
# TODO finish volcano plots (provide same highlights as in EP plot)
# TODO define color concept
# TODO finish analysis for Burbulla
# Make enrich_clustermap

# TODO ()
# TODO provide simple highlighted functions for plotting

# Plotting settings
DPI = 300
FIG_FORMAT = "pdf"
ARGS_SAVE = dict(dpi=DPI, bbox_inches="tight")
LEGEND_FONTSIZE = 16

FILE_DAVID_DOWN = "David_DOWN.xlsx"
FILE_DAVID_UP = "David_UP.xlsx"
COL_E_SCORE = "E-score"
COL_P_SCORE = "P-score"
COL_C_SCORE = "C-score"
COL_CS_STD = "CS_STD"
COL_PE_MEAN = "PE-mean"


# I Helper Functions
def _get_scoring(df_p=None, df_e=None):
    """"""
    df_p[COL_P_SCORE] = xo.p_score(ids=df_p["Gene names"],
                                   x_fc=df_p["Log2FC"],
                                   x_pvals=df_p["pvalue"])
    df_p[COL_E_SCORE] = xo.e_score(ids=df_p["Gene names"],
                                   id_lists=df_e["Genes"],
                                   x_fe=df_e["Fold Enrichment"],
                                   x_pvals=df_e["pval"])
    df_p[COL_PE_MEAN] = df_p[[COL_E_SCORE, COL_P_SCORE]].mean(axis=1)
    df_p = df_p.sort_values(by=COL_PE_MEAN, ascending=False)
    return df_p


# II Main Functions
def impute():
    """"""
    # Load _data
    file = "data/raw_data_proteomics_lfq.xlsx"
    df_raw = pd.read_excel(file)
    # Settings
    #groups = ["A", "B", "C", "D", "E", "F"]
    groups = ["A", "B"]
    str_lfq = "log2 LFQ "
    str_ids = "Protein IDs"
    pp = xo.PreProcess(str_id=str_ids, str_lfq=str_lfq)
    # Creat imputation object
    cimp = xo.cImpute(str_id=str_ids, str_lfq=str_lfq)
    dict_group_cols = pp.get_dict_groups(df=df_raw, groups=groups)
    all_groups_col = pp.get_all_group_cols(dict_group_cols=dict_group_cols)
    # Imputation
    d_min, up_mnar, d_max = cimp.get_limits(df=df_raw.copy(),
                                            group_cols=all_groups_col)

    df_imp = cimp.run(df=df_raw.copy(),
                      dict_group_cols=dict_group_cols)
    df_imp = pp.filter(df=df_imp, drop_na=False)
    # Plot histogram
    xo.plot_settings()
    xo.plot_imput_histo(df_raw=df_raw,
                        df_imp=df_imp,
                        cols=all_groups_col,
                        d_min=d_min,
                        up_mnar=up_mnar)
    plt.show()
    col_fc = "log2 FC (A/B)"
    df = pp.run(df=df_imp, groups=groups, ids=df_imp.index.to_list())
    df = pp.filter(df=df, drop_na=True)
    print(len(df), len(df[df[col_fc] > 0]), len(df[df[col_fc] < 0]))
    _df_raw = pp.run(df=df_raw, groups=groups, ids=df_raw.index.to_list())
    _df_raw = pp.filter(df=_df_raw, drop_na=True)
    print(len(_df_raw), len(_df_raw[_df_raw[col_fc] > 0]), len(_df_raw[_df_raw[col_fc] < 0]))
    xo.plot_volcano(df=_df_raw, col_fc=col_fc, col_pval="-log10 p-value (A/B)")
    plt.title("Raw")
    plt.show()
    xo.plot_volcano(df=df, col_fc=col_fc, col_pval="-log10 p-value (A/B)")
    plt.title("Imputed")
    plt.show()

    # Plot scatter plot for each group#
    """
    df_raw_plot = df_raw.set_index(str_ids)
    df_raw_plot = df_raw_plot.sort_index()
    for group in dict_group_cols:
        cols = dict_group_cols[group]
        xo.plot_imput_scatter(df_raw=df_raw_plot, df_imp=df_imp, cols=cols, group=group)
        plt.tight_layout()
        plt.show()
        plt.close()
    """


def scoring_analysis(df_p=None, df_e=None):
    """"""
    """
    print(xo.p_score(ids=['protein1', 'protein2'], x_fc=[2.4, 1.5], x_pvals=[0.05, 0.2]))
    print(xo.e_score(ids=['protein1', 'protein2'], id_lists=[['protein1', 'protein2'], ['protein2']],
                     x_fe=[2, 1.5], x_pvals=[0.05, 0.1]))
    """
    df = _get_scoring(df_p=df_p, df_e=df_e)
    df = df[df["pvalue"] < 0.05]
    xo.plot_settings(weight_bold=False, short_ticks_x=True, short_ticks_y=True, font_scale=1)
    xo.plot_explano(df=df.copy(), show_names=True, force_points=0.5, force_text=0.3, force_object=0.2)
    plt.tight_layout()


def e_hit_analysis(df_p=None, df_e=None):
    """"""
    #df = get_scoring(df_p=df_p, df_e=df_e)
    #df = df[df["pvalue"] < 0.05]
    list_genes = ["NQO1", "STX4", "CLCN3", "ROBO2"]
    terms_sub_list = ['Cytoplasm', 'Cytosol', 'Dendrite', 'Fluid shear stress and atherosclerosis',
                      'Identical protein binding', 'NAD(P)H dehydrogenase (quinone) activity',
                      'NADPH dehydrogenase (quinone) activity', 'Negative regulation of catalytic activity',
                      'Protein binding', 'RNA binding', 'Synapse']
    df_hits = xo.e_hits(ids=df_p["Gene names"], id_lists=df_e["Genes"], terms=df_e["Term"],
                        terms_sub_list=None, n_ids=10, n_terms=10, sort_alpha=True)
    xo.plot_settings()
    xo.plot_enrich_map(df=df_hits)
    plt.savefig(ut.FOLDER_RESULTS + "test.pdf")


# III Test/Caller Functions
def xo_pipeline():
    """"""
    df_p = pd.read_excel(ut.FOLDER_DATA + FILE_BPAN)
    df_eu = pd.read_excel(ut.FOLDER_DATA + FILE_DAVID_UP)
    df_ed = pd.read_excel(ut.FOLDER_DATA + FILE_DAVID_DOWN)
    df_e = pd.concat([df_eu, df_ed], axis=0)
    # Pre-processing (filtering)
    pp = xo.PreProcess()
    print(df_p)
    df_p = pp.filter(df=df_p, drop_na=False, cols="Gene names", split_names="Gene names")
    print(df_p)
    df_e = pp.filter(df=df_e, drop_na=False, cols="Term")
    scoring_analysis(df_p=df_p, df_e=df_e)
    plt.show()
    e_hit_analysis(df_p=df_p, df_e=df_e)


# IV Main
def main():
    t0 = time.time()
    impute()
    xo_pipeline()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
