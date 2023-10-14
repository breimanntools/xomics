"""
This is a script for pre the interface of preprocessing omics data.
"""
import pandas as pd
import numpy as np
from datetime import datetime
from typing import Optional

import xomics.utils as ut
from ._backend.preprocess_run import run_preprocess
from ._backend.preprocess_filter import filter_df

# TODO finish testing, test on real data in dev_scripts


# I Helper Functions
def check_match_df_groups(df=None, groups=None):
    """"""
    groups_df = [x.split("_") for x in list(df)]
    wrong_groups = [x for x in groups if x not in groups_df]
    if len(wrong_groups) > 0:
        raise ValueError(f"The following entries from 'groups' are not in 'df': {wrong_groups}")


def check_match_df_ids(df=None, ids=None, str_id=None):
    """"""
    if ids is None:
        if str_id not in list(df):
            raise ValueError(f"'ids' should not be None if 'str_id' ({str_id}) attribute not 'df' columns: {list(df)}")
        else:
            return df[str_id].to_list()
    if len(ids) != len(df):
        raise ValueError(f"'ids' (len={len(ids)}) should have the same length as 'df' (len={len(df)}")
    if isinstance(ids, pd.Series):
        return ids.to_list()
    elif isinstance(ids, list):
        return ids
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
                 str_id: str = "Protein IDs",
                 str_quant: str = "log2 LFQ"
                 ):
        """
        Initialize the PreProcessor object with specific string identifiers for ID and LFQ columns.

        Parameters:
        -----------
        str_id
            Identifier for the protein ID column in the dataframe.
        str_quant
            Identifier for the LFQ columns in the dataframe.
        """
        self.str_id = str_id
        self.str_quant = str_quant

    @ut.doc_params(doc_param_df_groups=doc_param_df_groups)
    def get_dict_col_quant_group(self,
                                 df: pd.DataFrame = None,
                                 groups: list = None
                                 ) -> dict:
        """
        Create a dictionary with groups from df based on lfq_str and given groups

        Parameters
        ----------
        {doc_param_df_groups}

        Return
        ------
        dict_col_quant_group
            Dictionary assigning names of columns with quantifications (keys) to group names (values)
        """
        # Check input
        df = ut.check_df(df=df, accept_none=False)
        groups = ut.check_list_like(name="groups", val=groups, accept_none=False)
        check_match_df_groups(groups=groups, df=df)
        # Get dictionary for quantifications
        dict_col_quant_group = ut.get_dict_col_quant_group(df=df, groups=groups, str_quant=self.str_quant)
        return dict_col_quant_group

    @ut.doc_params(doc_param_df_groups=doc_param_df_groups)
    def get_dict_group_cols_quant(self,
                                  df: pd.DataFrame = None,
                                  groups: list = None
                                  ) -> dict:
        """
        Create a dictionary with groups from df based on lfq_str and given groups

        Parameters
        ----------
        {doc_param_df_groups}

        Return
        ------
        dict_group_cols_quant
            Dictionary assigning group names (keys) to list of columns with quantifications (values)
        """
        # Check input
        df = ut.check_df(df=df, accept_none=False)
        groups = ut.check_list_like(name="groups", val=groups, accept_none=False)
        check_match_df_groups(groups=groups, df=df)
        # Get dictionary for quantifications
        dict_group_cols_quant = ut.get_dict_group_cols_quant(df=df, groups=groups, str_quant=self.str_quant)
        return dict_group_cols_quant

    @ut.doc_params(doc_param_df_groups=doc_param_df_groups)
    def get_cols_quant(self,
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
        check_match_df_groups(groups=groups, df=df)
        # Get columns with quantifications
        cols_quant = ut.get_cols_quant(df=df, groups=groups, str_quant=self.str_quant)
        return cols_quant

    @staticmethod
    def filter(df: pd.DataFrame = None,
               cols: list = None,
               drop_na: bool = True,
               col_split_names: Optional[str] = None,
               str_split: Optional[str] = ";"
               ) -> pd.DataFrame:
        """
        Filter and optionally modify the provided DataFrame based on specified parameters.

        Parameters:
        -----------
        df
            The DataFrame to filter and optionally modify.
        cols
            List of column names to consider for filtering.
        drop_na
            Whether to drop rows containing NaN values in the specified `cols`.
        col_split_names
            Name of the column in which to perform string splitting (if applicable).
        str_split
            The string character(s) to use for splitting string values in ``col_split_names``.

        Returns:
        --------
        pd.DataFrame
            The filtered (and optionally modified) DataFrame.

        Notes:
        ------
        - The function performs a series of checks to validate the provided parameters
          before proceeding with the filtering. It uses utility methods from an object
          `ut` for these checks, which is not defined in the provided code snippet.
        - After validation, a `filter_df` function (also not defined in the provided code)
          is called to perform the actual filtering based on the checked/processed parameters.
        - Raises ValueError: if validation checks fail for `df`, `cols`, `drop_na`,
          `col_split_names`, or `str_split`.
        """
        ut.check_df(df=df)
        cols = ut.check_list_like(name="cols", val=cols, accept_none=False, accept_str=True)
        ut.check_col_in_df(df=df, name_df="df", cols=cols, accept_nan=True)
        ut.check_bool(name="drop_na", val=drop_na)
        col_split_names = ut.check_col_in_df(df=df, name_df="df", cols=col_split_names, accept_nan=True, accept_none=True)
        ut.check_str(name="str_split", val=str_split)
        # Filtering
        df = filter_df(df=df, cols=cols, drop_na=drop_na, col_split_names=col_split_names, str_split=str_split)
        return df

    @staticmethod
    def adjust_log(df: pd.DataFrame = None,
                   cols: list = None,
                   log2: bool = True,
                   neg: bool = False
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

    @ut.doc_params(doc_param_df_groups=doc_param_df_groups)
    def run(self,
            df: pd.DataFrame = None,
            groups: list = None,
            ids: ut.ArrayLike1D = None,
            drop_na: bool = False,
            pvals_method: Optional[str] = None,
            pvals_neg_log10: bool = True
            ) -> pd.DataFrame:
        """
        Perform pairwise t-tests for groups to obtain -log10 p-values and log2 fold changes,
        with optional p-value correction, nan policy, and log-scale output.

        Parameters:
        -----------
        {doc_param_df_groups}
        ids
            List or array of protein identifiers.
        drop_na
            Whether to drop rows with missing values
        pvals_method
            Correction method for t-tests {"bonferroni", "sidak", "holm", "hommel", "fdr_bh"}.
        pvals_neg_log10
            Whether to return p-values in -log10 scale.

        Returns:
        --------
        df_fc
            DataFrame with p-values and log2 fold changes for each group comparison.
        """
        # Check functions
        df = ut.check_df(df=df, accept_none=False)
        groups = ut.check_list_like(name="groups", val=groups, accept_none=False)
        ut.check_list_like(name="ids", val=ids, accept_none=True)
        ut.check_str_in_list(name="pvals_method", val=pvals_method,
                             list_options=["bonferroni", "sidak", "holm", "hommel", "fdr_bh"])
        ut.check_bool(name="pvals_neg_log10", val=pvals_neg_log10)
        ut.check_bool(name="drop_na", val=drop_na)
        check_match_df_groups(groups=groups, df=df)
        check_match_df_ids(ids=ids, df=df, str_id=self.str_id)
        # Get the mapping dictionaries
        df_fc = run_preprocess(df=df, ids=ids, groups=groups, drop_na=drop_na,
                               pvals_method=pvals_method, pvals_neg_log10=pvals_neg_log10,
                               str_id=self.str_id)
        return df_fc
