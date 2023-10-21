"""
This is a script for backend of the PreProcess.run() method.
"""
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
import warnings

import xomics.utils as ut


# I Helper Functions
def _calculate_group_stats(df, group_cols):
    """Calculate the mean for each group's columns."""
    return df[group_cols].mean(axis=1)


def _correct_p_val(p_vals=None, method=None):
    """Correct p values with given methods"""
    p_vals = [p if str(p) != "nan" else 1 for p in p_vals]
    cor_p_vals = multipletests(p_vals, method=method)[1]
    p_vals = [p if str(p_vals[i]) != "nan" else np.nan for i, p in enumerate(cor_p_vals)]
    return p_vals


# II Main Functions
def run_preprocess(df=None, groups=None, groups_ctrl=None, pvals_method=None, pvals_neg_log10=True,  str_quant=None):
    """
    Perform pairwise t-tests for groups to obtain -log10 p-values and log2 fold changes,
    with optional p-value correction, nan policy, and log-scale output.
    """
    df = df.copy()
    # Get the mapping dictionaries
    dict_group_cols_quant = ut.get_dict_group_qcols(df=df, groups=groups, str_quant=str_quant)
    cols_quant = ut.get_qcols(df=df, groups=groups, str_quant=str_quant)
    # Remove rows with NaNs
    df_quant = df[cols_quant]
    # Initialize lists to collect results
    log2_FC_columns = {}
    p_value_columns = {}
    ratio_pairs = set()

    for group in groups:
        for group_ctrl in groups_ctrl:
            # Use frozenset to ensure the pair is unique regardless of order
            pair = frozenset([group, group_ctrl])
            if group == group_ctrl or pair in ratio_pairs:
                continue
            # Calculate mean for each group
            mean1 = _calculate_group_stats(df_quant, dict_group_cols_quant[group])
            mean2 = _calculate_group_stats(df_quant, dict_group_cols_quant[group_ctrl])

            # Calculate log2 fold change and p-values
            fold_change = mean2 - mean1
            # Ignore RuntimeWarning due to missing values
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                _, p_values = ttest_ind(df_quant[dict_group_cols_quant[group]],
                                        df_quant[dict_group_cols_quant[group_ctrl]],
                                        axis=1, nan_policy="omit")
            # Correct p-values if method is specified
            if pvals_method is not None:
                p_values = _correct_p_val(p_vals=p_values, method=pvals_method)
            # Create column names
            log2_FC_col_name = f"{ut.STR_FC}_({group}/{group_ctrl})"
            p_value_col_name = f"{ut.STR_PVAL}_({group}/{group_ctrl})"

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
    df_fc = pd.DataFrame(results)
    return df_fc
