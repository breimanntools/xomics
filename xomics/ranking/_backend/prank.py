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
        mean, std = np.nanmean(x.flatten()), np.nanstd(x.flatten())
        z_scores = (x - mean) / std
        return z_scores
    # Min-max normalization
    min_val, max_val = np.nanmin(x), np.nanmax(x)
    if min_val == max_val:
        return np.full(x.shape, 0.5)  # Return array with constant value if all inputs are identical
    min_max_scores = (x - min_val) / (max_val - min_val)
    return min_max_scores


def _normalize_folds(x_vals=None, z_norm=True):
    """Normalize fold changes or fold enrichment."""
    # Scale
    if min(x_vals) < 0:
        x_vals = abs(x_vals)
    # Normalize
    x_vals = _normalize_values(x_vals, z_norm=z_norm)
    return x_vals


# Ranking functions
def _p_ranking(x_fc, x_pvals, z_norm=False):
    """
    Rank proteins based on their fold changes and significance.

    Parameters
    x_fc: array-like, fold change scores for each protein
    x_pvals: array-like, p-values associated with each protein
    z_norm: boolean, whether to apply z-normalization

    Returns
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

    Parameters
    x_fe: array-like, fold enrichment scores for each protein
    x_pvals: array-like, p-values associated with each protein
    x_hit: binary matrix denoting presence (1) or absence (0) of each protein
    z_norm: boolean, whether to apply z-normalization

    Returns
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


def _e_ranking_only_pvals(x_pvals, x_hit, z_norm=False):
    """
    Rank proteins based on their  significance.

    Parameters
    x_pvals: array-like, p-values associated with each protein
    x_hit: binary matrix denoting presence (1) or absence (0) of each protein
    z_norm: boolean, whether to apply z-normalization

    Returns
    ranking_score: array-like, combined and normalized scores for ranking
    """
    # This ensures that each hit has a positive contribution to the ranking score
    # The addition of 0.00001  avoids identical values, which can cause problems in min-max normalization
    # Values <= 0.001 have same impact (empirically tested)
    x_pvals += abs(min(x_pvals)) + 0.00001
    # Weight the combined values by presence (1) or absence (0) of each protein
    x_s = np.transpose(x_hit) * x_pvals
    # Sum the weighted scores for each protein
    x_s = x_s.sum(axis=1)
    # Normalize the ranking scores
    ranking_score = _normalize_values(x_s, z_norm=z_norm)
    return ranking_score


# II Main Functions
def p_score(x_fc=None, x_pvals=None):
    """Calculate the single protein use_cases ranking score (P score)."""
    # Normalize data
    norm_fc = _normalize_folds(x_vals=x_fc, z_norm=True)
    norm_pvals = _normalize_values(x_pvals, z_norm=True)
    # Scoring (min-max normalized)
    p_scores = np.array(_p_ranking(norm_fc, norm_pvals))
    return p_scores


def e_score(names=None, name_lists=None, x_fe=None, x_pval=None):
    """Calculate the single protein enrichment score (E score)."""
    # Normalize data
    norm_pvals = _normalize_values(x_pval, z_norm=True)
    norm_fe = _normalize_folds(x_vals=x_fe, z_norm=True)
    # Get unique protein IDs from input sets
    unique_ids = ut.flatten_list(list_in=name_lists, sep=",")
    # Create binary hit matrix to represent the presence of unique IDs in each set
    x_hit = [[int(x in id_set) for x in unique_ids] for id_set in name_lists]
    # Scoring for unique IDs (min-max normalized)
    _ranking_scores = _e_ranking(norm_fe, norm_pvals, x_hit)
    # Map unique IDs to their final scores
    dict_val = dict(zip(unique_ids, _ranking_scores))
    e_scores = np.array([dict_val.get(i, 0) for i in names])
    return e_scores


def e_score_only_pvals(names=None, name_lists=None, x_pval=None):
    """Calculate the single protein enrichment score (E score)."""
    # Normalize data
    norm_pvals = _normalize_values(x_pval, z_norm=True)
    # Get unique protein IDs from input sets
    unique_ids = ut.flatten_list(list_in=name_lists, sep=",")
    # Create binary hit matrix to represent the presence of unique IDs in each set
    x_hit = [[int(x in id_set) for x in unique_ids] for id_set in name_lists]
    # Scoring for unique IDs (min-max normalized)
    _ranking_scores = _e_ranking_only_pvals(norm_pvals, x_hit)
    # Map unique IDs to their final scores
    dict_val = dict(zip(unique_ids, _ranking_scores))
    e_scores = np.array([dict_val.get(i, 0) for i in names])
    return e_scores



def c_score(ids=None, df_imp=None, col_id=None):
    """Obtain protein use_cases confidence score (C score) from cImpute output"""
    list_ids = df_imp.index.to_list() if col_id is None else df_imp[col_id].to_list()
    dict_c_scores = dict(zip(list_ids, df_imp[ut.COL_C_SCORE]))
    c_scores = [dict_c_scores[i] for i in ids]
    return c_scores
