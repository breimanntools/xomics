"""
This is a script for the backend of the ranking methods of the pRank (protein-centric ranking) class.
"""
import pandas as pd
import numpy as np

import xomics.utils as ut

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions
# Data transformation functions
def _normalize_values(x, z_norm=True):
    """Normalize the given array."""
    # Z normalization
    if z_norm:
        mean, std = x.flatten().mean(), x.flatten().std()
        z_scores = (x - mean) / std
        return z_scores
    # Min-max normalization
    if min(x) == max(x):
        return np.full(x.shape, 0.5)    # Return array with constant value if all inputs are identical
    min_max_scores = (x - min(x)) / (max(x) - min(x))
    return min_max_scores


def _adjust_log2_xvals(x_vals):
    """Heuristic test and conduction of log2 transformation for fold changes/enrichment values
    (performed for absolute values)"""
    # Check if transformation is needed
    if 10 < max(x_vals):
        x_vals = np.log2(x_vals)
        return x_vals, True
    else:
        return x_vals, False


def _adjust_log10_pvals(x_pvals):
    """Heuristic test and conduction of -log10 transformation for p-values"""
    # Check if transformation is needed
    p_log10 = (max(x_pvals) - min(x_pvals)) > 1
    if not p_log10:
        x_pvals = -np.log10(x_pvals)
        return x_pvals, True
    else:
        return x_pvals, False


def _normalize_pvals(x_pvals=None, z_norm=True, adjust_log=True, for_fc=True, verbose=True):
    """Normalize p-values (adjust if log10 transformed)."""
    # Scale
    if adjust_log:
        x_pvals, log_transformed = _adjust_log10_pvals(x_pvals)
        if log_transformed and verbose:
            str_fc_fe = "fold change" if for_fc else "fold enrichment"
            ut.print_out(f"-log10 transformation conducted for p-values of {str_fc_fe}.")
    # Normalize
    x_pvals = _normalize_values(x_pvals, z_norm=z_norm)
    return x_pvals


def _normalize_folds(x_vals=None, adjust_log=False, z_norm=True, for_fc=True, verbose=True):
    """Normalize fold changes or fold enrichment."""
    # Scale
    if min(x_vals) < 0:
        x_vals = abs(x_vals)
    if adjust_log:
        x_vals, log_transformed = _adjust_log2_xvals(x_vals)
        if log_transformed and verbose:
            str_fc_fe = "fold change" if for_fc else "fold enrichment"
            ut.print_out(f"Log2 transformation conducted for {str_fc_fe} values.")
    # Normalize
    x_vals = _normalize_values(x_vals, z_norm=z_norm)
    return x_vals


# Ranking functions
def _p_ranking(x_fc, x_pvals, z_norm=False):
    """
    Rank proteins based on their fold changes and significance.

    Parameters:
    x_fc: array-like, fold change scores for each protein
    x_pvals: array-like, p-values associated with each protein
    z_norm: boolean, whether to apply z-normalization

    Returns:
    ranking_score: array-like, combined and normalized scores for ranking
    """
    # Combine the fold change and p-values for each protein.
    # This generates a new metric that encapsulates both variance and significance.
    x_s = x_fc + x_pvals
    # Normalize the combined scores using either Min-Max or Z-score normalization.
    # This ensures that the final scores are in a comparable range.
    ranking_score = _normalize_values(x_s, z_norm=z_norm)
    return ranking_score


def _e_ranking(x_fe, x_pvals, x_hit, z_norm=False):
    """
    Rank proteins based on their fold enrichment and significance.

    Parameters:
    x_fe: array-like, fold enrichment scores for each protein
    x_pvals: array-like, p-values associated with each protein
    x_hit: binary matrix denoting presence (1) or absence (0) of each protein
    z_norm: boolean, whether to apply z-normalization

    Returns:
    ranking_score: array-like, combined and normalized scores for ranking
    """
    # Shift fold enrichment and p-values to positive domain
    # This ensures that each hit has a positive contribution to the ranking score
    # The addition of 0.00001  avoids identical values, which can cause problems in min-max normalization
    # Values <= 0.001 have same impact (empirically tested)
    x_fe += abs(min(x_fe)) + 0.00001
    x_pvals += abs(min(x_pvals)) + 0.00001
    # Combine the shifted fold enrichment and p-values for each protein
    x_fe_p = x_fe + x_pvals
    # Weight the combined values by presence (1) or absence (0) of each protein
    x_s = np.transpose(x_hit) * x_fe_p
    # Sum the weighted scores for each protein
    x_s = x_s.sum(axis=1)
    # Normalize the ranking scores
    ranking_score = _normalize_values(x_s, z_norm=z_norm)
    return ranking_score


# II Main Functions
def p_score(ids=None, x_fc=None, x_pvals=None, adjust_log=True, verbose=True):
    """Calculate the single protein proteomics ranking score (P score)."""
    # Normalize data
    args = dict(adjust_log=adjust_log, z_norm=True, verbose=verbose)
    norm_fc = _normalize_folds(x_vals=x_fc, **args)
    norm_pvals = _normalize_pvals(x_pvals=x_pvals, **args)
    # Scoring (min-max normalized)
    p_scores = np.array(_p_ranking(norm_fc, norm_pvals))
    return p_scores


def e_score(ids=None, id_lists=None, x_fe=None, x_pvals=None, adjust_log=True, verbose=True):
    """Calculate the single protein enrichment score (E score)."""
    # Normalize data
    args = dict(z_norm=True, adjust_log=adjust_log, for_fc=False, verbose=verbose)   # For DPBM originally z_norm=False
    x_pvals = _normalize_pvals(x_pvals=x_pvals, **args)
    x_fe = _normalize_folds(x_vals=x_fe, **args)
    # Get unique protein IDs from input sets
    unique_ids = ut.flatten_list(list_in=id_lists, sep=",")
    # Create binary hit matrix to represent the presence of unique IDs in each set
    x_hit = [[int(x in id_set) for x in unique_ids] for id_set in id_lists]
    # Scoring for unique IDs (min-max normalized)
    _ranking_scores = _e_ranking(x_fe, x_pvals, x_hit)
    # Map unique IDs to their final scores
    dict_val = dict(zip(unique_ids, _ranking_scores))
    e_scores = np.array([dict_val.get(i, 0) for i in ids])
    return e_scores


def c_score(ids=None, df_imp=None, col_id=None):
    """Obtain protein proteomics confidence score (C score) from cImpute output"""
    list_ids = df_imp.index.to_list() if col_id is None else df_imp[col_id].to_list()
    dict_c_scores = dict(zip(list_ids, df_imp[ut.COL_C_SCORE]))
    c_scores = [dict_c_scores[i] for i in ids]
    return c_scores
