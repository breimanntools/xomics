"""
This is a script for pre the interface of preprocessing omics data.
"""
import pandas as pd
import numpy as np
from typing import Optional

import xomics.utils as ut
from ._backend.preprocess_run import run_preprocess
from ._backend.preprocess_filter import filter_df, filter_names

# TODO finish testing, test on real data in dev_scripts
# TODO Filter for number of quantifications


# I Helper Functions
def check_base(base=None):
    """Ensure 'base' is a valid numerical type and has an acceptable value"""
    if not isinstance(base, (int, float)) or base not in [2, 10]:
        raise ValueError("'base' must be a numerical value and either 2 or 10")


def check_match_df_ids(df=None, list_ids=None):
    """"""

    if list_ids is None:
        raise ValueError("'list_ids' should not be None")
    if len(list_ids) != len(df):
        raise ValueError(f"'ids' (len={len(list_ids)}) should have the same length as 'df' (len={len(df)}")
    if isinstance(list_ids, pd.Series):
        return list_ids.to_list()
    elif isinstance(list_ids, list):
        return list_ids
    else:
        raise ValueError("'ids' should be a list or a pandas.DataFrame column")


# II Main Functions
# Common interface
doc_param_df_groups = \
"""\
df
    DataFrame with quantifications. ``Rows`` typically correspond to proteins and ``columns`` to conditions.
groups: list of str
    List with names grouping conditions from ``df`` columns.\
"""


class PreProcess:
    def __init__(self,
                 str_id: str = "protein_id",
                 str_quant: str = "log2_lfq"
                 ):
        """
        Initialize the PreProcessor object with specific string identifiers for ID and LFQ columns.

        Parameters
        ----------
        str_id
            Identifier for the protein ID column in the dataframe.
        str_quant
            Identifier for the LFQ columns in the dataframe.
        """
        self.str_id = str_id
        self.str_quant = str_quant

    @ut.doc_params(doc_param_df_groups=doc_param_df_groups)
    def get_qcols(self,
                  df: pd.DataFrame = None,
                  groups: list = None
                  ) -> list:
        """
        Create a list with groups from df based on str_quant and given groups

        Parameters
        ----------
        {doc_param_df_groups}

        Return
        ------
        cols_quant
            List with all quantification columns across all groups
        """
        # Check input
        df = ut.check_df(df=df, accept_none=False)
        groups = ut.check_list_like(name="groups", val=groups, accept_none=False)
        ut.check_match_df_groups(groups=groups, df=df, str_quant=self.str_quant)
        # Get columns with quantifications
        cols_quant = ut.get_qcols(df=df, groups=groups, str_quant=self.str_quant)
        return cols_quant

    @ut.doc_params(doc_param_df_groups=doc_param_df_groups)
    def get_dict_qcol_group(self,
                            df: pd.DataFrame = None,
                            groups: list = None
                            ) -> dict:
        """
        Create a dictionary with quantification columns and the group they are subordinated to

        Parameters
        ----------
        {doc_param_df_groups}

        Return
        ------
        dict_qcol_group
            Dictionary assigning names of columns with quantifications (keys) to group names (values)
        """
        # Check input
        df = ut.check_df(df=df, accept_none=False)
        groups = ut.check_list_like(name="groups", val=groups, accept_none=False)
        ut.check_match_df_groups(groups=groups, df=df, str_quant=self.str_quant)
        # Get dictionary for quantifications
        dict_qcol_group = ut.get_dict_qcol_group(df=df, groups=groups, str_quant=self.str_quant)
        return dict_qcol_group

    @ut.doc_params(doc_param_df_groups=doc_param_df_groups)
    def get_dict_group_qcols(self,
                             df: pd.DataFrame = None,
                             groups: list = None
                             ) -> dict:
        """
        Create a dictionary with for groups from df and their corresponding columns with quantifications

        Parameters
        ----------
        {doc_param_df_groups}

        Return
        ------
        dict_group_qcols
            Dictionary assigning group names (keys) to list of columns with quantifications (values)
        """
        # Check input
        df = ut.check_df(df=df, accept_none=False)
        groups = ut.check_list_like(name="groups", val=groups, accept_none=False)
        ut.check_match_df_groups(groups=groups, df=df, str_quant=self.str_quant)
        # Get dictionary for quantifications
        dict_group_qcols = ut.get_dict_group_qcols(df=df, groups=groups, str_quant=self.str_quant)
        return dict_group_qcols

    @staticmethod
    def filter(df: pd.DataFrame = None,
               cols: list = None,
               drop_na: bool = True,
               ) -> pd.DataFrame:
        """
        Filter and optionally modify the provided DataFrame based on specified parameters.

        Parameters
        ----------
        df
            The DataFrame with quantifications to filter. ``Rows`` typically correspond to proteins
            and ``columns`` to conditions.
        cols
            List of columns from ``df`` to consider for filtering.
        drop_na
            Whether to drop rows containing NaN values in the specified `cols`.

        Returns
        -------
        df
            The filtered DataFrame.
        """
        ut.check_df(df=df)
        cols = ut.check_list_like(name="cols", val=cols, accept_none=True, accept_str=True)
        ut.check_col_in_df(df=df, name_df="df", cols=cols, accept_none=True, accept_nan=True)
        ut.check_bool(name="drop_na", val=drop_na)
        # Filtering
        df = filter_df(df=df, cols=cols, drop_na=drop_na)
        return df

    @staticmethod
    def filter_names(df: pd.DataFrame = None,
                     col: list = None,
                     str_split: str = ";",
                     drop_na: bool = True,
                     ) -> pd.DataFrame:
        """
        Split names from column in dataframe and filter duplicates.

        Parameters
        ----------
        df
            The DataFrame with quantifications to filter. ``Rows`` typically correspond to proteins
            and ``columns`` to conditions.
        col
            Column from ``df`` in which to perform string splitting and filtering.
        str_split
            The string character(s) to use for splitting string values in ``col``.
        drop_na
            Whether to drop rows containing NaN values in the specified `cols`.

        Returns
        -------
        df
            The modified and filtered DataFrame.
        """
        ut.check_df(df=df)
        ut.check_col_in_df(df=df, name_df="df", cols=col, accept_nan=True, accept_none=False)
        ut.check_str(name="str_split", val=str_split)
        ut.check_bool(name="drop_na", val=drop_na)
        # Filtering
        df = filter_names(df=df, col=col, str_split=str_split, drop_na=drop_na)
        return df

    @staticmethod
    def apply_log(df: pd.DataFrame = None,
                  cols: list = None,
                  log2: bool = True,
                  neg: bool = False,
                  ) -> pd.DataFrame:
        """
        Apply a logarithmic transformation to specified columns of a DataFrame.

        Parameters
        ----------
        df
            DataFrame containing data to transform.
        cols
            Names of columns to apply the logarithmic transformation to.
        log2
            If True, apply a log2 transformation. Otherwise, apply a log10 transformation.
        neg
            If True, multiply the logarithmic result by -1.

        Returns
        -------
        df
            DataFrame with specified columns log-transformed.

        Notes
        -----
        - Make sure all values in the specified columns of df are > 0 before applying this function.
        - NaN values will remain NaN after the transformation.
        """
        # Check input
        df = ut.check_df(df=df, all_positive=True)
        if cols is None:
            cols = list(df)
        cols = ut.check_list_like(name="cols", val=cols, accept_none=False, accept_str=True)
        ut.check_col_in_df(df=df, cols=cols, accept_nan=True)
        ut.check_bool(name="log2", val=log2)
        ut.check_bool(name="neg", val=neg)
        # Log transform
        f = np.log2 if log2 else np.log10
        df[cols] = df[cols].apply(f)
        if neg:
            df[cols] *= -1
        return df

    @staticmethod
    def apply_exp(df: pd.DataFrame = None,
                  cols: list = None,
                  base: int = 2.0,
                  neg: bool = False,
                  ) -> pd.DataFrame:
        """
        Apply an exponential transformation to specified columns of a DataFrame.

        Parameters
        ----------
        df
            DataFrame containing data to transform.
        cols
            Names of columns to apply the exponential transformation to.
        base
            The base of the exponential function. If ``base=2``, apply a 2**x transformation,
            otherwise apply a 10**x transformation if ``base=10``.
        neg
            If True, multiply the exponential result by -1.

        Returns
        -------
        df
            DataFrame with specified columns exponentially transformed.

        Notes
        -----
        - NaN values will remain NaN after the transformation.
        """
        # Check input
        df = ut.check_df(df=df)
        if cols is None:
            cols = list(df)
        cols = ut.check_list_like(name="cols", val=cols, accept_none=False, accept_str=True)
        ut.check_col_in_df(df=df, cols=cols, accept_nan=True)
        ut.check_bool(name="neg", val=neg)
        check_base(base=base)
        # Exponential transform
        if neg:
            df[cols] *= -1
        df[cols] = df[cols].apply(lambda x: np.power(base, x))
        return df

    @staticmethod
    def add_id(df=None, list_ids=None, col_name_to_add="Protein IDs"):
        """
        Add column with protein ids to DataFrame.

        Parameters
        ----------
        df
            DataFrame containing fold-change and p-values.
        list_ids
            List or array of protein/gene identifiers.
        col_name_to_add
            The name of the column for the protein/gene identifiers to be added.

        Returns
        -------
        df
            DataFrame with added significance column.
        """
        # Check input
        df = ut.check_df(df=df, accept_none=False)
        ut.check_list_like(name="list_ids", val=list_ids)
        list_ids = check_match_df_ids(df=df, list_ids=list_ids)
        ut.check_str(name="col_name", val=col_name_to_add, accept_none=False)
        # Add column with ids
        df.insert(0, col_name_to_add, list_ids)
        return df

    @staticmethod
    def add_significance(df: pd.DataFrame = None,
                         col_fc: str = None,
                         col_pval: str = None,
                         th_fc: float = 0.5,
                         th_pval: float = 0.05
                         ) -> pd.DataFrame:
        """Add a column indicating significance regarding threshold for fold change and p-value

        Three types of significance classes are defined:

        - **Up**: Significant hits that are 'up-regulated'(i.e., right quadrant of volcano plot)
        - **Down**: Significant hits that are 'down-regulated' (i.e., left quadrant of volcano plot)
        - **Not Sig.**: Hits that are not significant.

        Parameters
        ----------
        df
            DataFrame containing fold-change and p-values.
        col_fc
            Column name containing fold change values.
        col_pval
            Column name containing p-values.
        th_fc
            Threshold for fold-change, applied for negative and positive values.
        th_pval
            Threshold for p-value, -log10 transformed before applied.

        Returns
        -------
        df
            DataFrame with added significance column.
        """
        # Check input
        df = ut.check_df(name="df", df=df, cols_req=[col_fc, col_pval])
        ut.check_number_range(name="th_fc", val=th_fc, min_val=0, just_int=False)
        ut.check_number_range(name="th_pval", val=th_pval, min_val=0, max_val=1, just_int=False)
        # Rescale p-value
        th_pval = -np.log10(th_pval)
        # Add significant classes (Up, Down, Not Sig.)
        df[ut.COL_SIG_CLASS] = ut.get_sig_classes(df=df, col_fc=col_fc, col_pval=col_pval, th_pval=th_pval, th_fc=th_fc)
        return df

    @ut.doc_params(doc_param_df_groups=doc_param_df_groups)
    def run(self,
            df: pd.DataFrame = None,
            groups: list = None,
            groups_ctrl: list = None,
            pvals_correction: Optional[str] = None,
            pvals_neg_log10: bool = True
            ) -> pd.DataFrame:
        """
        Perform pairwise t-tests for groups to obtain -log10 p-values and log2 fold changes,
        with optional p-value correction, nan policy, and log-scale output.

        Parameters
        ----------
        {doc_param_df_groups}
        groups_ctrl
            List with names control grouping conditions from ``df`` columns.
        pvals_correction
            Correction method for t-tests {"bonferroni", "sidak", "holm", "hommel", "fdr_bh"}.
        pvals_neg_log10
            Whether to return p-values in -log10 scale.

        Returns
        -------
        df_fc
            DataFrame with p-values and log2 fold changes for each group comparison.

        Notes
        -----
        Fold changes (FC) and P-values will be computed for each group in ``groups`` compared against
        each group in ``group_ctrl`` (group/group_ctrl), where self-comparison is omitted.
        """
        # Check input
        df = ut.check_df(df=df, accept_none=False)
        groups = ut.check_list_like(name="groups", val=groups, accept_none=False)
        if groups_ctrl is None:
            groups_ctrl = groups
        groups_ctrl = ut.check_list_like(name="groups_ctrl", val=groups_ctrl)
        ut.check_str_in_list(name="pvals_method", val=pvals_correction, accept_none=True,
                             list_options=["bonferroni", "sidak", "holm", "hommel", "fdr_bh"])
        ut.check_bool(name="pvals_neg_log10", val=pvals_neg_log10)
        ut.check_match_df_groups(df=df, groups=groups, str_quant=self.str_quant)
        ut.check_match_df_groups(df=df, groups=groups_ctrl, name_groups="groups_ctrl", str_quant=self.str_quant)
        # Get the mapping dictionaries
        df_fc = run_preprocess(df=df, groups=groups, groups_ctrl=groups_ctrl,
                               pvals_method=pvals_correction, pvals_neg_log10=pvals_neg_log10,
                               str_quant=self.str_quant)
        return df_fc
