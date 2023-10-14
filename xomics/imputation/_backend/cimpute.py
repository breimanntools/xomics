"""
This is a script for backend of the cImpute class
"""
import pandas as pd
import numpy as np
from scipy.stats import truncnorm
from sklearn.impute import KNNImputer

import xomics.utils as ut


# I Helper Functions
def _compute_cs(vals=None, mv_class=None, n=None):
    """Compute confidence score (CS) depending on missing value category and
    proportion of missing values"""
    count_nan = lambda val: len([x for x in val if str(x) == "nan"])
    if mv_class == ut.STR_NM:
        return 1
    elif mv_class == ut.STR_MAR:
        return 0
    elif mv_class == ut.STR_MCAR:
        return round((n - count_nan(vals)) / n, 2)
    elif mv_class == ut.STR_MNAR:
        return round(count_nan(vals) / n, 2)


def _impute_mnr(df=None, std_factor=0.5, d_min=None, up_mnar=None):
    """MinProb imputation as suggested by Lazar et al., 2016

    Arguments
    --------
    df: DataFrame
        DataFrame with missing values just classified as MNAR
    std_factor: int, default = 0.5
        Factor to control size of standard deviation of distribution relative to distance of upMNAR and Dmin.

    Notes
    -----
    std = (up_mnar - d_min) * std_factor

    See also
    --------
    https://www.rdocumentation.org/packages/imputeLCMD/versions/2.0/topics/impute.MinProb
    https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/MissingValues.html
    """
    df = df.copy()
    d1, d2 = df.shape
    # TODO optimize
    scale = (up_mnar - d_min/2) # * std_factor   # Standard deviation (spread, scale, or "width")
    # Generate random numbers using truncated (left-censored) normal distribution
    vals = truncnorm.rvs(a=0, b=1, size=d1*d2, loc=d_min, scale=scale).reshape((d1, d2))
    mask = df.isnull()
    df[mask] = vals
    return df


def _impute_mcar(df=None, n_neighbors=6):
    """KNN imputation via sklearn implementation

    Arguments
    ---------
    df: DataFrame
        DataFrame with missing values just classified as MNAR
    n neighbors: int, default=6 (Liu and Dongre, 2020)
        Number of neighboring samples to use for imputation
    """
    imputer = KNNImputer(n_neighbors=n_neighbors)
    X = np.array(df)
    index, cols = df.index, df.columns
    X = imputer.fit_transform(X)
    df = pd.DataFrame(X, columns=cols, index=index)
    return df


def _impute(df=None, mv_class=None, d_min=None, up_mnar=None, std_factor=0.5, n_neighbors=6):
    """Wrapper for imputation methods applied on an experimental group"""
    if mv_class == ut.STR_NM:
        return df
    elif mv_class == ut.STR_MAR:
        return df
    elif mv_class == ut.STR_MCAR:
        return _impute_mcar(df=df, n_neighbors=n_neighbors)
    elif mv_class == ut.STR_MNAR:
        return _impute_mnr(df=df, d_min=d_min, std_factor=std_factor, up_mnar=up_mnar)


# II Main Functions
def get_up_mnar(df=None, loc_pct_up_mnar=0.25):
    """Get upper bound for MNAR MVs for whole data set"""
    d_min = df.min().min()  # Detection limit
    d_max = df.max().max()  # Largest detected value
    dr = d_max - d_min      # Detection range
    up_mnar = d_min + loc_pct_up_mnar * dr   # Upper MNAR border
    return d_min, up_mnar


def classify_of_mvs(df_group=None, up_mnar=None):
    """Classification of missing values for given protein intensities of an experimental group"""
    n_groups = len(list(df_group))
    mv_classes = []     # NaN Classes
    for i, row in df_group.iterrows():
        n_nan = row.isnull().sum()
        n_higher_up_mnar = np.array((row > up_mnar)).sum()
        n_lower_or_equal_up_mnar = np.array((row <= up_mnar)).sum()
        # MNAR (Missing Not At Random)
        if n_lower_or_equal_up_mnar + n_nan == n_groups:
            mv_classes.append(ut.STR_MNAR)
        # MCAR (Missing Completely At Random)
        elif n_higher_up_mnar + n_nan == n_groups:
            mv_classes.append(ut.STR_MCAR)
        # NM (No Missing values)
        elif n_higher_up_mnar + n_lower_or_equal_up_mnar == n_groups:
            mv_classes.append(ut.STR_NM)
        # MAR (Missing At Random)
        else:
            mv_classes.append(ut.STR_MAR)
    return mv_classes


def compute_cs(df_group=None, mv_classes=None):
    """Computation of confidence scores depending on missing value classification and proportion of missing values"""
    n = len(list(df_group))
    dict_prot_cs = {}
    for mv_class in ut.LIST_MV_CLASSES:
        mask = [True if l == mv_class else False for l in mv_classes]
        df = df_group[mask]
        for entry, row in df.iterrows():
            cs = _compute_cs(vals=list(row), mv_class=mv_class, n=n)
            dict_prot_cs[entry] = cs
    list_cs = [dict_prot_cs[entry] for entry in df_group.index]
    return list_cs


def impute(df_group=None, mv_classes=None, list_cs=None, min_cs=0.5, d_min=None, up_mnar=None,
           n_neighbors=5, std_factor=0.5):
    """Group-wise imputation over whole data set"""
    df_group = df_group.copy()
    list_df = []
    for mv_class in ut.LIST_MV_CLASSES:
        mask = np.array([True if (l == mv_class and cs >= min_cs) else False for l, cs in zip(mv_classes, list_cs)])
        df = df_group[mask]
        df_imped = _impute(df=df, mv_class=mv_class,
                             n_neighbors=n_neighbors, std_factor=std_factor,
                             d_min=d_min, up_mnar=up_mnar)
        list_df.append(df_imped)
    df_group_imputed = pd.concat(list_df, axis=0).sort_index()
    mask = np.array([True if i in df_group_imputed.index else False for i in df_group.index])
    df_group[mask] = df_group_imputed
    return df_group


