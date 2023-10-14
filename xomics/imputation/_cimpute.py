"""
cImpute (conditional Imputation) is a hybrid imputation algorithm for missing values (MVs) in (prote)omics data.
Missing values can be distinguished into three categories as described by Lazar et al., 2016 and Wei et al., 2018
for proteomic data sets as follows:

    a) Missing Completely At Random (MCAR): MVs due to random errors and stochastic fluctuations during process of
        data acquisition. Since MCAR MVs can not be explained by measured intensities, they are uniformly distributed.
    b) Missing At Random (MAR): MVs due to suboptimal data processing and conditional dependencies. MAR is a more
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
import xomics.utils as ut

from ._backend.cimpute import run_cimpute, get_up_mnar


# TODO a) generalize (e.g., lfq -> intensities, test with other input)
# TODO b) benchmarking (ring trail, standards (Tenzer Lab), existing artificial benchmark sets (c.f. publication)
# TODO c) optimize (e.g., optimize of clusters for KNN via Silhouette score)
# TODO d) Extend to other omics data

# I Helper Functions


# II Main Functions
class cImpute:
    """Hybrid imputation algorithm for missing values (MVs) in (prote)omics data.

    Parameters
    ----------
    str_id: str, default = "Protein IDs"
        Column name of entry ids of input DataFrame for associated methods
    str_quant: str, default = "log2 LFQ"
        Common substring of intensity columns of input DataFrame for associated methods

    """
    def __init__(self, str_id="Protein IDs", str_quant="log2 LFQ"):
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
            cols_quant = ut.get_cols_quant(df=df, groups=groups, str_quant=self.str_quant)
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



