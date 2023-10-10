"""
This is a script for ...
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from datetime import datetime

import xomics._utils as ut

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe

# TODO test on real _data in dev_scripts, make volcano, use for cImpute, finish remaining pipleline


# I Helper Functions
def _check_p_correction(method=None):
    """Check p value correction methods"""
    p_corrections = ["bonferroni", "sidak", "holm", "hommel", "fdr_bh"]
    if method is not None and method not in p_corrections:
        raise ValueError("P value correction should be one of following: " + str(p_corrections))


def _correct_p_val(p_vals=None, method=None):
    """Correct p values with given methods"""
    p_vals = [p if str(p) != "nan" else 1 for p in p_vals]
    cor_p_vals = multipletests(p_vals, method=method)[1]
    p_vals = [p if str(p_vals[i]) != "nan" else np.nan for i, p in enumerate(cor_p_vals)]
    return p_vals


def _check_ids(ids=None, df=None, str_id=None):
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


def _calculate_group_stats(df, group_cols):
    """Calculate the mean for each group's columns."""
    return df[group_cols].mean(axis=1)


# II Main Functions
class PreProcess:
    def __init__(self, str_id="Protein IDs", str_lfq="log2 LFQ"):
        """
        Initialize the PreProcessor object with specific string identifiers for ID and LFQ columns.

        Parameters:
        -----------
        str_id : str, default = "Protein IDs"
            Identifier for the protein ID column in the dataframe.
        str_lfq : str, default = "log2 LFQ"
            Identifier for the LFQ columns in the dataframe.
        """
        self.str_id = str_id
        self.str_lfq = str_lfq

    @staticmethod
    def filter(df=None, cols=None, drop_na=True, split_names=None):
        """
        Filter _data frame for non-correct or missing values
        """
        # Split names for provided column and filter duplicates by keeping first occuring
        if split_names is not None:
            if split_names not in list(df):
                raise ValueError(f"'split_names' ({split_names}) should be in columns of 'df': {list(df)}")
            else:
                f = lambda x: x.split(";")[0] if type(x) is str and ";" in x else x
                df[split_names] = [f(x) for x in df[split_names]]
                df = df.drop_duplicates(subset=split_names, keep="first")
        # TODO add filter redundant elements in list
        # Filter duplicates and datatime columns
        def _filter_df(_df, _col):
            """"""
            if _col not in list(_df):
                raise ValueError(f"{_col} from 'cols' should bein columns of 'df': {list(df)}")
            _df = _df[_df[_col].apply(lambda x: not isinstance(x, (datetime, type(np.nan))))]
            _df = _df.drop_duplicates(subset=_col, keep="first")
            return _df
        if isinstance(cols, list):
            for col in cols:
                df = _filter_df(df, col)
        elif isinstance(cols, str):
            df = _filter_df(df, cols)
        # Filter missing values
        if drop_na:
            if cols is not None:
                df = df.dropna(subset=cols)
            else:
                df = df.dropna()
        return df

    @staticmethod
    def adjust_log(df=None, cols=None, auto=True):
        """"""


    def get_dict_groups(self, df=None, groups=None, group_to_col=True):
        """
        Create a dictionary with groups from df based on lfq_str and given groups

        Parameters
        ----------
        df: DataFrame
            DataFrame containing quantified features including missing values
        groups: list of str
            List with group names
        group_to_col: bool, default = True
            Decide whether group columns (values) should be returned as list

        Return
        ------
        dict_col_group: dict
            Dictionary assigning column names (keys) to group names (values) if group_to_col=False
        dict_group_cols: dict
            Dictionary assigning groups (keys) to list of column names (values) if group_to_col=True
        """
        dict_col_group = {}
        for col in list(df):
            if self.str_lfq in col:
                col_wo_lfq_str = col.replace(self.str_lfq, "")
                for group in groups:
                    if group in col_wo_lfq_str:
                        dict_col_group[col] = group
        if group_to_col:
            dict_group_cols = {g: [k for k, v in dict_col_group.items() if v == g] for g in groups}
            return dict_group_cols
        return dict_col_group

    @staticmethod
    def get_all_group_cols(dict_group_cols=None):
        """Retrieve all columns from group dictionary

        Parameters
        ----------
        dict_group_cols: dict
            Dictionary assigning groups (keys) to list of column names (values)

        Return
        ------
        all_group_cols: list
            List with all columns of input Data Frame from cImpute.get_dict_groups containing intensity values.
        """
        all_group_cols = []
        for group in dict_group_cols:
            all_group_cols.extend(dict_group_cols[group])
        return all_group_cols

    def run(self, df=None, ids=None, groups=None, drop_na=False, pvals_method=None, pvals_neg_log10=True):
        """
        Perform pairwise t-tests for groups to obtain -log10 p-values and log2 fold changes,
        with optional p-value correction, nan policy, and log-scale output.

        Parameters:
        -----------
        df : DataFrame
            The input dataframe with sample values (e.g., LFQ intensities, transcript expression)
        ids : list or array-like
            List or array of protein identifiers.
        groups : list, array-like
            List of group names
        drop_na : bool, default=False
            Whether to drop rows with missing values
        pvals_method : str, default=None
            Correction method for t-tests {None, "bonferroni", "sidak", "holm", "hommel", "fdr_bh"}.
        pvals_neg_log10 : bool, default=True
            Whether to return p-values in -log10 scale.

        Returns:
        --------
        dict
            Dictionary containing the p-values and log2 fold changes for each group comparison.
        """
        # Check functions
        _check_p_correction(method=pvals_method)
        _check_ids(ids=ids, df=df, str_id=self.str_id)
        # Get the mapping dictionaries
        dict_group_cols = self.get_dict_groups(df=df, groups=groups)
        all_group_cols = self.get_all_group_cols(dict_group_cols=dict_group_cols)

        # Remove rows with NaNs if policy is 'omit'
        df_filtered = df[all_group_cols]
        if drop_na:
            df_filtered = df_filtered.dropna()
        # Initialize result dictionary
        # Initialize lists to collect results
        log2_FC_columns = {}
        p_value_columns = {}
        ratio_pairs = set()

        for group1 in groups:
            for group2 in groups:
                # Use frozenset to ensure the pair is unique regardless of order
                pair = frozenset([group1, group2])
                if group1 == group2 or pair in ratio_pairs:
                    continue
                # Calculate mean for each group
                mean1 = _calculate_group_stats(df_filtered, dict_group_cols[group1])
                mean2 = _calculate_group_stats(df_filtered, dict_group_cols[group2])

                # Calculate log2 fold change and p-values
                fold_change = mean2 - mean1
                _, p_values = ttest_ind(df_filtered[dict_group_cols[group1]],
                                        df_filtered[dict_group_cols[group2]],
                                        axis=1,
                                        nan_policy="omit")

                # Correct p-values if method is specified
                if pvals_method is not None:
                    p_values = _correct_p_val(p_vals=p_values, method=pvals_method)
                # Create column names
                log2_FC_col_name = f"{ut.STR_FC} ({group1}/{group2})"
                p_value_col_name = f"{ut.STR_PVAL} ({group1}/{group2})"

                # Convert pandas Series or numpy arrays to lists
                fold_change_list = fold_change.tolist()
                p_value_list = (-np.log10(p_values) if pvals_neg_log10 else p_values).tolist()

                # Append the results to the respective lists
                log2_FC_columns[log2_FC_col_name] = fold_change_list
                p_value_columns[p_value_col_name] = p_value_list

                # Add to ratio pairs to avoid duplicate comparisons
                ratio_pairs.add(pair)
        # Combine both dictionaries
        results = {**log2_FC_columns, **p_value_columns}
        # Convert dictionary to DataFrame
        df = pd.DataFrame(results)
        df.insert(0, self.str_id, ids)
        return df
