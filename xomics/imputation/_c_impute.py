"""
cImpute (conditional Imputation) is a hybrid imputation algorithm for missing values (MVs) in (prote)omics _data.
Missing values can be distinguished into three categories as described by Lazar et al., 2016 and Wei et al., 2018
for proteomic _data sets as follows:

    a) Missing Completely At Random (MCAR): MVs due to random errors and stochastic fluctuations during process of
        _data acquisition. Since MCAR MVs can not be explained by measured intensities, they are uniformly distributed.
    b) Missing At Random (MAR): MVs due to suboptimal _data processing and conditional dependencies. MAR is a more
        general class than MCAR, where all MCAR MVs are MAR MVs. The distribution of MAR MVs can just be speculated
        and likely differs highly between experiments.
    c) Missing Not At Random (MNAR): MVs due to experimental bias (i.e., the detection limit in mass-spectrometry
        experiments). Commonly, MNAR MVs are described by a left-censored Gaussian distribution (i.e., the Gaussian
        distribution is truncated on the region of lower abundances, which is the left side of the distribution).

cImpute aims to impute only MVs matching well-defined confidence criteria and consists of following four steps:

    1. Definition of upper bound for MNAR MVs to distinguish between MNAR and MCAR
    2. Classification of MVs for detected proteins in an experimental group
    3. Computation of confidence score (CS)
    4. Group-wise imputation for proteins with CS higher than given threshold

For more details look into the README

"""
import pandas as pd
import numpy as np
import xomics._utils as ut

from scipy.stats import truncnorm
from sklearn.impute import KNNImputer

# TODO a) add check functions for interface
# TODO b) generalize (e.g., lfq -> intensities, test with other input)
# TODO c) optimize (e.g., optimize of clusters for KNN via Silhouette score)
# TODO d) testing
# TODO e) extend documentation & create google colab
# TODO f) benchmarking (ring trail, standards (Tenzer Lab), existing artificial benchmark sets (c.f. publication)
# TODO g) Extend to other omics _data

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe

# Constants
STR_MCAR = "MCAR"
STR_MNAR = "MNAR"
STR_MAR = "MAR"
STR_NM = "NM"
LIST_MV_CLASSES = [STR_MCAR, STR_MNAR, STR_MAR, STR_NM]
STR_CS = "CS"
STR_MV_LABELS = "labels"


# I Helper Functions
def _get_all_group_cols(dict_group_cols=None):
    """Retrieve all columns from group dictionary"""
    all_group_cols = []
    for group in dict_group_cols:
        all_group_cols.extend(dict_group_cols[group])
    return all_group_cols


def _compute_cs(vals=None, mv_class=None, n=None):
    """Compute confidence score (CS) depending on missing value category and
    proportion of missing values"""
    count_nan = lambda val: len([x for x in val if str(x) == "nan"])
    if mv_class == STR_NM:
        return 1
    elif mv_class == STR_MAR:
        return 0
    elif mv_class == STR_MCAR:
        return round((n - count_nan(vals)) / n, 2)
    elif mv_class == STR_MNAR:
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
    if mv_class == STR_NM:
        return df
    elif mv_class == STR_MAR:
        return df
    elif mv_class == STR_MCAR:
        return _impute_mcar(df=df, n_neighbors=n_neighbors)
    elif mv_class == STR_MNAR:
        return _impute_mnr(df=df, d_min=d_min, std_factor=std_factor, up_mnar=up_mnar)


# II Main Functions
def get_up_mnar(df=None, loc_up_mnar=0.25):
    """Get upper bound for MNAR MVs for whole _data set"""
    d_min = df.min().min()  # Detection limit
    d_max = df.max().max()  # Largest detected value
    dr = d_max - d_min      # Detection range
    up_mnar = d_min + loc_up_mnar * dr   # Upper MNAR border
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
            mv_classes.append(STR_MNAR)
        # MCAR (Missing Completely At Random)
        elif n_higher_up_mnar + n_nan == n_groups:
            mv_classes.append(STR_MCAR)
        # NM (No Missing values)
        elif n_higher_up_mnar + n_lower_or_equal_up_mnar == n_groups:
            mv_classes.append(STR_NM)
        # MAR (Missing At Random)
        else:
            mv_classes.append(STR_MAR)
    return mv_classes


def compute_cs(df_group=None, mv_classes=None):
    """Computation of confidence scores depending on missing value classification and proportion of missing values"""
    n = len(list(df_group))
    dict_prot_cs = {}
    for mv_class in LIST_MV_CLASSES:
        mask = [True if l == mv_class else False for l in mv_classes]
        df = df_group[mask]
        for entry, row in df.iterrows():
            cs = _compute_cs(vals=list(row), mv_class=mv_class, n=n)
            dict_prot_cs[entry] = cs
    list_cs = [dict_prot_cs[entry] for entry in df_group.index]
    return list_cs


def impute(df_group=None, mv_classes=None, list_cs=None, min_cs=0.5, d_min=None, up_mnar=None,
           n_neighbors=5, std_factor=0.5):
    """Group-wise imputation over whole _data set"""
    df_group = df_group.copy()
    list_df = []
    for mv_class in LIST_MV_CLASSES:
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


# Wrapper
class cImpute:
    """Hybrid imputation algorithm for missing values (MVs) in (prote)omics _data.

    Parameters
    ----------
    str_id: str, default = "Protein IDs"
        Column name of entry ids of input DataFrame for associated methods
    str_lfq: str, default = "log2 LFQ"
        Common substring of intensity columns of input DataFrame for associated methods

    """
    def __init__(self, str_id="Protein IDs", str_lfq="log2 LFQ"):
        self.list_mv_classes = LIST_MV_CLASSES
        self.str_id = str_id
        self.str_lfq = str_lfq

    @staticmethod
    def get_limits(df=None, loc_up_mnar=0.25, group_cols=None):
        """Get minimum of detected values (d_min, i.e., detection limit), upper bound of MNAR MVs (up_mnar),
        and maximum of detected values (d_max).

        Parameters
        ----------
        df: DataFrame
            DataFrame containing quantified features including missing values
        loc_up_mnar: int, default=0.1, [0-1]
            Location factor for the upMNAR given as relative proportion of the detection range
        group_cols: list or array-like
            List of all group columns

        Return
        ------
        d_min: int
            Minimum of detected values
        up_mnar: int
            upper bound of MNAR MVs
        d_max: int
            Maximum of detected values
        """
        df = df.copy()
        d_min, up_mnar = get_up_mnar(df=df[group_cols], loc_up_mnar=loc_up_mnar)
        d_max = df[group_cols].max().max()
        return d_min, up_mnar, d_max

    def run(self, df=None, dict_group_cols=None, loc_up_mnar=0.25, min_cs=0.5, std_factor=0.5, n_neighbors=5):
        """Hybrid method for imputation of omics _data called conditional imputation (cImpute)
        using MinProb for MNAR (Missing Not at Random) missing values and KNN imputation for
        MCAR (Missing completely at Random) missing values.

        Parameters
        ----------
        df: DataFrame
            DataFrame containing quantified features including missing values
        dict_group_cols: dict
            Dictionary for groups to list of columns containing intensity values for group.
        min_cs: int, default 0.5 [0-1]
            Minimum of confidence score used for selecting values for protein in groups to apply imputation on.
        loc_up_mnar: int, default 0.25 [0-1]
            Factor to determine the location of the upper detection limit bound. In percent of total value range.
        std_factor: int, default = 0.5 (MinProb parameter)
            Factor to control size of standard deviation of left-censored distribution
            relative to distance of upMNAR and Dmin.
        n_neighbors : int, default=5 (KNN imputation parameter)
            Number of neighboring samples to use for imputation.

        Return
        ------
        df_imp: DataFrame
            DataFrame with (a) imputed intensities values and (b) group-wise confidence score and NaN classification.
        """
        # TODO refactor, simplify, optimize n clusters, optimize for performance, use arrays
        df = df.copy()
        df.index = df[self.str_id]
        df = df.sort_index()
        all_group_cols = _get_all_group_cols(dict_group_cols=dict_group_cols)
        d_min, up_mnar = get_up_mnar(df=df[all_group_cols], loc_up_mnar=loc_up_mnar)
        # TODO change to numpy Arrays & compute summary statistic (n MVs per class and group)
        list_df_groups = []
        list_mv_classes = []
        cs_vals = []

        for group in dict_group_cols:
            cols = dict_group_cols[group]
            df_group = df[cols]
            mv_classes = classify_of_mvs(df_group=df_group, up_mnar=up_mnar)
            list_cs = compute_cs(df_group=df_group, mv_classes=mv_classes)
            df_group = impute(df_group=df_group, mv_classes=mv_classes, list_cs=list_cs,
                              min_cs=min_cs,
                              d_min=d_min, up_mnar=up_mnar,
                              std_factor=std_factor,
                              n_neighbors=n_neighbors)
            list_df_groups.append(df_group)
            list_mv_classes.append(mv_classes)
            cs_vals.append(list_cs)

        # Merge imputation for all groups
        df_imp = pd.concat(list_df_groups, axis=1)
        # Add aggregated CS values (mean and std)
        cs_means = np.array(cs_vals).mean(axis=0).round(2)
        cs_stds = np.array(cs_vals).std(axis=0).round(2)
        df_imp.insert(len(list(df_imp)), ut.COL_C_SCORE, cs_means)
        df_imp.insert(len(list(df_imp)), ut.COL_CS_STD, cs_stds)
        # Add  CS values per group
        df_cs = pd.DataFrame(cs_vals).T
        df_cs.columns = [f"CS_{group}" for group in dict_group_cols]
        df_cs.index = df.index
        df_nan = pd.DataFrame(list_mv_classes).T
        df_nan.columns = [f"NaN_{group}" for group in dict_group_cols]
        df_nan.index = df.index
        df_imp = pd.concat([df_imp, df_cs, df_nan], axis=1)
        return df_imp
