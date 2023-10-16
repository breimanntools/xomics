"""
This is a script for backend of the PreProcess.filter() method.
"""
import numpy as np
from datetime import datetime


# Helper functions


def _filter_invalid_and_duplicates(df, col):
    """
    Filter out invalid (datetime and nan) and duplicate entries based on a column.
    """
    if col not in df.columns:
        raise ValueError(f"{col} from 'cols' should be in columns of 'df': {list(df)}")
    df = df[df[col].apply(lambda x: not isinstance(x, (datetime, type(np.nan))))]
    return df.drop_duplicates(subset=col, keep="first")


# Main functions
def filter_df(df=None, cols=None, drop_na=True):
    """
    Filter DataFrame according to specified criteria and columns.
    """
    if cols is not None:
        if isinstance(cols, list):
            for col in cols:
                df = _filter_invalid_and_duplicates(df, col)
        elif isinstance(cols, str):
            df = _filter_invalid_and_duplicates(df, cols)
    if drop_na:
        columns = cols if cols is not None else df.columns
        df = df.dropna(subset=columns)
    return df


def filter_names(df=None, col=None, str_split=";", drop_na=True):
    """
    Split names for provided column and filter duplicates by keeping first occurring.
    """
    df[col] = df[col].apply(lambda x: x.split(str_split)[0] if isinstance(x, str) and str_split in x else x)
    if drop_na:
        df = df.drop_duplicates(subset=col, keep="first")
    return df

