"""
This is a script for backend of the PreProcess.filter() method.
"""
import numpy as np
from datetime import datetime


# Helper functions
def _filter_split_names(df, split_names, str_split=";"):
    """
    Split names for provided column and filter duplicates by keeping first occurring.
    """
    if split_names not in df.columns:
        raise ValueError(f"'split_names' ({split_names}) should be in columns of 'df': {list(df)}")
    df[split_names] = df[split_names].apply(lambda x: x.split(str_split)[0] if isinstance(x, str) and str_split in x else x)
    return df.drop_duplicates(subset=split_names, keep="first")


def _filter_invalid_and_duplicates(df, col):
    """
    Filter out invalid (datetime and nan) and duplicate entries based on a column.
    """
    if col not in df.columns:
        raise ValueError(f"{col} from 'cols' should be in columns of 'df': {list(df)}")
    df = df[df[col].apply(lambda x: not isinstance(x, (datetime, type(np.nan))))]
    return df.drop_duplicates(subset=col, keep="first")


def _filter_na(df, cols):
    """
    Filter out rows based on NA values in specified columns.
    """
    return df.dropna(subset=cols)


def filter_df(df, cols=None, drop_na=True, col_split_names=None, str_split=";"):
    """
    Filter DataFrame according to specified criteria and columns.
    """
    if col_split_names is not None:
        df = _filter_split_names(df, col_split_names, str_split=str_split)
    if cols is not None:
        if isinstance(cols, list):
            for col in cols:
                df = _filter_invalid_and_duplicates(df, col)
        elif isinstance(cols, str):
            df = _filter_invalid_and_duplicates(df, cols)
    if drop_na:
        df = _filter_na(df, cols if cols is not None else df.columns)
    return df
