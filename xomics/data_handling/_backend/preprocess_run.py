"""
This is a script for backend of the PreProcess.run() method.
"""
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind

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
def run_preprocess(df=None, ids=None, groups=None, drop_na=False, pvals_method=None, pvals_neg_log10=True, str_id=None):
    """
    Perform pairwise t-tests for groups to obtain -log10 p-values and log2 fold changes,
    with optional p-value correction, nan policy, and log-scale output.
    """
    # Get the mapping dictionaries
    dict_group_cols_quant = ut.get_dict_group_cols_quant(df=df, groups=groups)
    cols_quant = ut.get_cols_quant(df=df, groups=groups)
    # Remove rows with NaNs if policy is 'omit'
    df_filtered = df[cols_quant]
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
            mean1 = _calculate_group_stats(df_filtered, dict_group_cols_quant[group1])
            mean2 = _calculate_group_stats(df_filtered, dict_group_cols_quant[group2])

            # Calculate log2 fold change and p-values
            fold_change = mean2 - mean1
            _, p_values = ttest_ind(df_filtered[dict_group_cols_quant[group1]],
                                    df_filtered[dict_group_cols_quant[group2]], axis=1, nan_policy="omit")

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
    df.insert(0, str_id, ids)
    return df
