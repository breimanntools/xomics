"""
This is a script for ...
"""
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xomics as xo
import dev_scripts._utils as ut
# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions


# II Main Functions
def dev_data_pipeline():
    """"""
    # TODO add protein_id, gene_name to each df
    df_lip = xo.load_dataset(name="LIPID_DEMYLINATION") #drop_na=True, n=1000, random=True)
    df_lfq = xo.load_dataset(name="PROT_DEMYLINATION")
    groups = ["d00", "d03", "d07", "d14"]
    args_col = dict(col_fc="log2_fc_(d07/d00)", col_pval="-log10_p-value_(d07/d00)")

    args_lip = dict(col_id="lipid_id", col_name="lipid_name", str_quant="pmol")
    args_prot = dict(col_id="protein_id", col_name="gene_name", str_quant="log2_lfq")
    pp = xo.PreProcess(**args_lip)
    cols_quant = pp.get_qcols(df=df_lip, groups=groups)

    df_lip = pp.filter_duplicated_names(df=df_lip, col="lipid_id")
    df_lip = pp.filter_groups(df=df_lip, groups=groups)
    df_lip = pp.apply_log(df=df_lip, cols=cols_quant, log2=True)
    df_fc = pp.run(df=df_lip, groups=groups, groups_ctrl=["d00"])
    df_fc = pp.filter_nan(df=df_fc, cols=["log2_fc_(d07/d00)", "-log10_p-value_(d07/d00)"])
    df_fc = pp.add_significance(df=df_fc, **args_col)
    df_fc = df_fc.sort_values(by="sig_class")
    df_fc.to_excel(ut.FOLDER_RESULTS + "LIPID_DEMYLINATION_SIG.xlsx")

    df_fc = xo.pRank(**args_lip).p_score(df_fc=df_fc, **args_col)
    df_fc = df_fc.sort_values(by="p_score", ascending=False)
    names = df_fc["lipid_id"].to_list()[0:15]
    """
    xo.plot_settings()
    xo.plot_volcano(df=df_fc, **args_col, names_to_annotate=names, col_names="lipid_id")
    plt.xlim(-6, 6)
    plt.show()
    pp = xo.PreProcess(**args_prot)
    df_lfq = pp.filter_groups(df=df_lfq, groups=groups)
    df_fc = pp.run(df=df_lfq, groups=groups, groups_ctrl=["d00"])


    df_fc = pp.filter_nan(df=df_fc, cols=["log2_fc_(d07/d00)", "-log10_p-value_(d07/d00)"])
    df_fc = pp.filter_duplicated_names(df=df_fc, col="gene_name")
    df_fc = xo.pRank(**args_lip).p_score(df_fc=df_fc, **args_col)
    df_fc = pp.add_significance(df=df_fc, **args_col)
    df_fc = df_fc.sort_values(by="sig_class")
    df_fc.to_excel(ut.FOLDER_RESULTS + "PROT_DEMYLINATION_SIG.xlsx")
    df_fc = df_fc.sort_values(by="p_score", ascending=False)
    names = df_fc["gene_name"].to_list()[0:15]
    xo.plot_settings()
    df_fc = xo.pRank(**args_prot).p_score(df_fc=df_fc, **args_col)

    xo.plot_volcano(df=df_fc, **args_col, names_to_annotate=names, col_names="gene_name")
    plt.xlim(-6, 6)
    plt.show()
    """
    """
    cimp = xo.cImpute(col_id="lipid_id", col_name="lipid_name", str_quant="pmol")

    df_imp = cimp.run(df=df_lfq, groups=groups, min_cs=0.5)

    d_min, up_mnar, d_max = cimp.get_limits(df=df_lfq, groups=groups)

    xo.plot_settings()
    xo.plot_imput_histo(df_raw=df_lfq,
                        df_imp=df_lfq,
                        d_min=d_min,
                        up_mnar=up_mnar,
                        cols_quant=cols_quant)
    plt.ylim(0, 100)
    plt.show()
    df_raw = pp.run(df=df_lfq, groups=groups, groups_ctrl=["d00"])
    df_fc = pp.run(df=df_imp, groups=groups, groups_ctrl=["d00"])
    xo.plot_settings(weight_bold=False)
    df_fc = pp.add_significance(df=df_fc,
                                col_fc="log2_fc_(d07/d00)",
                                col_pval="-log10_p-value (d07/d00)")
    df_fc = pp.add_id(df=df_fc,
                      list_ids=df_lfq[col_gene],
                      col_name_to_add=col_gene)
    df_fc["P-Score"] = xo.pRank().p_score(x_fc=df_fc["log2_fc_(d07/d00)"],
                                          x_pval=df_fc["-log10_p-value (d07/d00)"])
    df_fc["C-Score"] = df_imp["C-Score"].values
    df_fc = df_fc.sort_values(by="P-Score", ascending=False)
    list_genes = df_fc[df_fc["sig_class"] == "Up"][col_gene].to_list()[0:10]
    df_raw = pp.add_id(df=df_raw,
                       list_ids=df_lfq[col_gene],
                       col_name_to_add=col_gene)
    xo.plot_volcano(df=df_raw, figsize=(7, 5),
                    col_fc="log2_fc_(d07/d00)",
                    col_pval="-log10_p-value (d07/d00)",
                    col_names=col_gene,
                    names_to_annotate=list_genes,
                    edge_color="black",
                    alpha=0.5)
    plt.show()

    xo.plot_volcano(df=df_fc,
                    figsize=(7, 5),
                    col_fc="log2_fc_(d07/d00)",
                    col_pval="-log10_p-value (d07/d00)",
                    col_cbar="C-Score",
                    names_to_annotate=list_genes,
                    edge_color="black",
                    col_names=col_gene,
                    alpha=0.5)
    plt.show()
    xo.plot_volcano(df=df_fc,
                    figsize=(7, 5),
                    col_fc="log2_fc_(d07/d00)",
                    col_pval="-log10_p-value (d07/d00)",
                    names_to_annotate=list_genes,
                    edge_color="black",
                    col_names=col_gene,
                    alpha=0.5)
    plt.show()
    """

# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    dev_data_pipeline()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
