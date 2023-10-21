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
    # TODO adjust for MNAR to go up if number of values increase again (minimum if n nan == n/2)
    count_nan = lambda val: len([x for x in val if str(x) == "nan"])
    if mv_class == ut.STR_NM:
        return 1
    elif mv_class == ut.STR_MAR:
        return 0
    elif mv_class == ut.STR_MCAR:
        return round((n - count_nan(vals)) / n, 2)
    elif mv_class == ut.STR_MNAR:
        return round(count_nan(vals) / n, 2)


def _impute_mnr(df=None, d_min=None, up_mnar=None):
    """MinProb imputation as suggested by Lazar et al., 2016

    See also
    --------
    https://www.rdocumentation.org/packages/imputeLCMD/versions/2.0/topics/impute.MinProb
    https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/MissingValues.html
    """
    df = df.copy()
    d1, d2 = df.shape
    # TODO optimize or justify scale factor (middle of MNAR so far. Could be adjusted based on std of distribution)
    # Factor to control size of standard deviation of left-censored distribution
    scale = (up_mnar - d_min/2)
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


def _impute(df=None, mv_class=None, d_min=None, up_mnar=None, n_neighbors=6, min_cs=0.5):
    """Wrapper for imputation methods applied on an experimental group"""
    if mv_class == ut.STR_NM:
        return df
    elif mv_class == ut.STR_MAR:
        # Missing at random only imputed if cs == 0
        if min_cs == 0:
            return _impute_mcar(df=df, n_neighbors=n_neighbors)
        return df
    elif mv_class == ut.STR_MCAR:
        return _impute_mcar(df=df, n_neighbors=n_neighbors)
    elif mv_class == ut.STR_MNAR:
        return _impute_mnr(df=df, d_min=d_min, up_mnar=up_mnar)


def _create_groupwise_dfs(cs_vals=None, mv_classes=None, group_dict=None, index=None, prefixes=None):
    """Creates and concatenates DataFrames for CS and NaN values per group."""
    df_cs = pd.DataFrame(cs_vals).T
    df_cs.columns = [f"{prefixes[0]}_{group}" for group in group_dict]

    df_nan = pd.DataFrame(mv_classes).T
    df_nan.columns = [f"{prefixes[1]}_{group}" for group in group_dict]

    both_dfs = pd.concat([df_cs, df_nan], axis=1)
    both_dfs.index = index
    return both_dfs


# II Main Functions
def get_up_mnar(df=None, loc_pct_upmnar=0.25):
    """Get upper bound for MNAR MVs for whole data set"""
    d_min = df.min().min()  # Detection limit
    d_max = df.max().max()  # Largest detected value
    dr = d_max - d_min      # Detection range
    up_mnar = d_min + loc_pct_upmnar * dr   # Upper MNAR border
    return d_min, up_mnar


def classify_of_mvs(df_group=None, up_mnar=None):
    """Classification of missing values for given protein intensities of an experimental group"""
    n_groups = len(list(df_group))
    mv_classes = []     # NaN Classes
    for i, row in df_group.iterrows():
        n_nan = row.isnull().sum()
        n_higher_up_mnar = np.array((row > up_mnar)).sum()
        n_lower_or_equal_up_mnar = np.array((row <= up_mnar)).sum()
        # NM (No Missing values)
        if n_nan == 0:
            mv_classes.append(ut.STR_NM)
        # MNAR (Missing Not At Random: all values are nan or lower than up_mnar)
        elif n_lower_or_equal_up_mnar + n_nan == n_groups:
            mv_classes.append(ut.STR_MNAR)
        # MCAR (Missing Completely At Random: all values are nan or higher than up_mnar)
        elif n_higher_up_mnar + n_nan == n_groups:
            mv_classes.append(ut.STR_MCAR)
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


def impute(df_group=None, mv_classes=None, list_cs=None, min_cs=0.5, d_min=None, up_mnar=None, n_neighbors=5):
    """Group-wise imputation over whole data set"""
    df_group = df_group.copy()
    list_df = []
    for mv_class in ut.LIST_MV_CLASSES:
        mask = np.array([True if (l == mv_class and cs >= min_cs) else False for l, cs in zip(mv_classes, list_cs)])
        df = df_group[mask]
        if len(df) > 0:
            df_imput = _impute(df=df, mv_class=mv_class,
                               n_neighbors=n_neighbors,
                               d_min=d_min,
                               up_mnar=up_mnar,
                               min_cs=min_cs)
            list_df.append(df_imput)
    df_group_imputed = pd.concat(list_df, axis=0).sort_index()
    mask = np.array([True if i in df_group_imputed.index else False for i in df_group.index])
    df_group[mask] = df_group_imputed
    return df_group


# TODO optimize n_neighbors, optimize for performance
# Main function
def run_cimpute(df=None, groups=None, min_cs=0.5, loc_pcat_upmnar=0.25, n_neighbors=5, str_id=None, str_quant=None):
    """Run complete cImpute pipeline"""
    df = df.copy()
    df.index = df[str_id]
    dict_group_cols_quant = ut.get_dict_group_qcols(df=df, groups=groups, str_quant=str_quant)
    cols_quant = ut.get_qcols(df=df, groups=groups, str_quant=str_quant)
    d_min, up_mnar = get_up_mnar(df=df[cols_quant], loc_pct_upmnar=loc_pcat_upmnar)
    list_df_groups = []
    list_mv_classes = []
    cs_vals = []
    for group in dict_group_cols_quant:
        cols_quant = dict_group_cols_quant[group]
        df_group = df[cols_quant]
        mv_classes = classify_of_mvs(df_group=df_group, up_mnar=up_mnar)
        list_cs = compute_cs(df_group=df_group, mv_classes=mv_classes)
        df_group = impute(df_group=df_group, mv_classes=mv_classes, list_cs=list_cs, min_cs=min_cs, d_min=d_min,
                          up_mnar=up_mnar, n_neighbors=n_neighbors)
        list_df_groups.append(df_group)
        list_mv_classes.append(mv_classes)
        cs_vals.append(list_cs)

    # Merge imputation for all groups
    df_imp = pd.concat(list_df_groups, axis=1)

    # Add aggregated CS values (mean and std)
    cs_means = np.array(cs_vals).mean(axis=0).round(2)
    cs_stds = np.array(cs_vals).std(axis=0).round(2)
    df_imp[ut.COL_C_SCORE], df_imp[ut.COL_C_STD] = cs_means, cs_stds

    # Add  CS values per group
    # Concatenate CS and NaN values per group and merge with df_imp
    df_cs_nan = _create_groupwise_dfs(cs_vals=cs_vals,
                                      mv_classes=list_mv_classes,
                                      group_dict=dict_group_cols_quant,
                                      index=df.index, prefixes=["CS", "MV"])
    df_imp = pd.concat([df_imp, df_cs_nan], axis=1)
    return df_imp
