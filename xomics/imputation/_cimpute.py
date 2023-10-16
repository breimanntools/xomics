"""
This is a script for interface of the cImpute (conditional Imputation) class.
"""
import pandas as pd
import numpy as np
import xomics.utils as ut

from ._backend.cimpute import run_cimpute, get_up_mnar


# TODO a) generalize (e.g., lfq -> intensities, test with other input)
# TODO b) benchmarking (ring trail, standards (Tenzer Lab), existing artificial benchmark sets (c.f. publication)
# TODO c) optimize (e.g., optimize of clusters for KNN via Silhouette score)
# TODO d) Extend to other omics data

# I Helper Functions


# II Main Functions
class cImpute:
    """Transparent hybrid data imputation method.

    `cImpute` (conditional Imputation) is a transparent hybrid imputation algorithm designed to address missing values
    (MVs) in (prote)omics data. The types of missing values can be broadly categorized into three groups based
    on their nature and the reasons behind their occurrence, as detailed in [Lazar16] and [Wei18]:

    - **Missing Completely At Random (MCAR)**: MVs resulting from random errors in data acquisition. Due to the inherent
      randomness of MCAR MVs, they cannot be explained purely by measured intensities and are uniformly distributed.
    - **Missing At Random (MAR)**: MVs due to data processing flaws and variable dependencies. While MAR includes all
      MCAR MVs, its distribution is speculative and varies significantly across different experiments.
    - **Missing Not At Random (MNAR)**: MVs caused by experimental biases like detection limits in mass-spectrometry.
      They often follow a left-censored Gaussian distribution, indicating truncation at lower abundances.

    Parameters
    ----------
    str_id: str, default = "Protein IDs"
        Column name of entry ids of input DataFrame for associated methods
    str_quant: str, default = "log2 LFQ"
        Common substring of intensity columns of input DataFrame for associated methods

    Notes
    -----
    The primary goal of `cImpute` is to focus on the imputation of MVs that align with well-defined confidence criteria.
    This approach comprises four main steps:

    1. Establishing the upper bound for MNAR MVs to distinguish between MNAR and MCAR.
    2. Categorizing MVs for detected proteins within a specific experimental group.
    3. Calculating a confidence score (CS) for each protein.
    4. Performing group-wise imputation for proteins whose CS exceeds a predefined threshold.

    """
    def __init__(self,
                 str_id: str = "protein_id",
                 str_quant: str = "log2_lfq"
                 ):
        self.list_mv_classes = ut.LIST_MV_CLASSES
        self.str_id = str_id
        self.str_quant = str_quant

    def get_limits(self, df=None, loc_pct_up_mnar=0.25, groups=None, cols_quant=None):
        """Get minimum of detected values (d_min, i.e., detection limit), upper bound of MNAR MVs (up_mnar),
        and maximum of detected values (d_max).

        Parameters
        ----------
        df: DataFrame
            DataFrame containing quantified features including missing values
        loc_pct_up_mnar: int, default=0.1, [0-1]
            Location factor for the upMNAR given as relative proportion (percentage) of the detection range
        groups: list or array-like
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
        if cols_quant is None:
            cols_quant = ut.get_qcols(df=df, groups=groups, str_quant=self.str_quant)
        d_min, up_mnar = get_up_mnar(df=df[cols_quant], loc_pct_up_mnar=loc_pct_up_mnar)
        d_max = df[cols_quant].max().max()
        return d_min, up_mnar, d_max

    def run(self, df=None, groups=None, min_cs=0.5, loc_up_mnar=0.25, n_neighbors=5):
        """Hybrid method for imputation of omics data called conditional imputation (cImpute)
        using MinProb for MNAR (Missing Not at Random) missing values and KNN imputation for
        MCAR (Missing completely at Random) missing values.

        Parameters
        ----------
        df: DataFrame
            DataFrame containing quantified features including missing values
        groups
            List of quantification group (substrings of columns in ``df``)
        min_cs: int, default 0.5 [0-1]
            Minimum of confidence score used for selecting values for protein in groups to apply imputation on.
        loc_up_mnar: int, default 0.25 [0-1]
            Factor to determine the location of the upper detection limit bound. In percent of total value range.

        n_neighbors : int, default=5 (KNN imputation parameter)
            Number of neighboring samples to use for imputation.

        Return
        ------
        df_imp: DataFrame
            DataFrame with (a) imputed intensities values and (b) group-wise confidence score and NaN classification.

        Notes
        -----
        MAR is only imputed if ``min_cs=0`` using the imputation for MCAR.
        """
        # Check input

        # Run cImpute algorithm
        df_imp = run_cimpute(df=df, groups=groups,
                             min_cs=min_cs, loc_up_mnar=loc_up_mnar, n_neighbors=n_neighbors,
                             str_quant=self.str_quant, str_id=self.str_id)
        return df_imp



