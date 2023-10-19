"""
This is a script for backend of the PreProcess.filter() method.
"""
import numpy as np
from datetime import datetime


# Helper functions
def _split_names(df=None, col=None, str_split=";"):
    """
    Split names for provided column and filter duplicates by keeping first occurring.
    """
    df[col] = df[col].apply(lambda x: x.split(str_split)[0] if isinstance(x, str) and str_split in x else x)
    return df


def _filter_invalid_and_duplicates(df, col):
    """
    Filter out invalid (datetime and nan) and duplicate entries based on a column.
    """
    if col not in df.columns:
        raise ValueError(f"{col} from 'cols' should be in columns of 'df': {list(df)}")
    df = df[df[col].apply(lambda x: not isinstance(x, (datetime, type(np.nan))))]
    return df.drop_duplicates(subset=col, keep="first")


# Main functions
def filter_duplicated_names(df=None, cols=None, split_names=False, str_split=";"):
    """
    Filter DataFrame according to specified criteria and columns.
    """
    if split_names:
        df = _split_names(df=df, col=cols, str_split=str_split)
    df = _filter_invalid_and_duplicates(df, cols)
    return df


def _compute_pct_non_nan(row, dict_groups_qcols, list_n_group):
    list_n_non_nan_group = [row[cols].notnull().sum() for cols in dict_groups_qcols.values()]
    return np.array(list_n_non_nan_group) / list_n_group


def filter_groups(df=None, groups=None, dict_groups_qcols=None, min_pct=None):
    """Filter df such that for at least one group a minimum percentage of values is given"""
    list_n_group = [len(dict_groups_qcols[g]) for g in groups]
    pct_non_nan_groups = df.apply(lambda row: _compute_pct_non_nan(row, dict_groups_qcols, list_n_group),
                                  axis=1, result_type='expand')
    mask = np.any(pct_non_nan_groups >= min_pct, axis=1)
    return df[mask]





