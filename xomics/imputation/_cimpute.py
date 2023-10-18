"""
This is a script for interface of the cImpute (conditional Imputation) class.
"""
import pandas as pd
import numpy as np
import xomics.utils as ut
from typing import Tuple

from ._backend.cimpute import run_cimpute, get_up_mnar


# TODO a) generalize (e.g., lfq -> intensities, test with other input)
# TODO b) benchmarking (ring trail, standards (Tenzer Lab), existing artificial benchmark sets (c.f. publication)
# TODO c) optimize (e.g., optimize of clusters for KNN via Silhouette score)
# TODO d) Extend to other omics data

# I Helper Functions
doc_param_df_groups_upmnar = \
"""\
df
    DataFrame containing quantified values with MVs. ``Rows`` typically correspond to proteins  and ``columns``  to conditions.
groups
    List of quantification group (substrings of columns in ``df``).
loc_pct_upmnar
    Location factor [0-1] for the upper MNAR limit (upMNAR) given as relative proportion (percentage) of the detection range.\
"""


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
        """
        Parameters
        ----------
        str_id
            Column name of entry ids of input DataFrame for associated methods
        str_quant
            Common substring of intensity columns of input DataFrame for associated methods
        """
        self.list_mv_classes = ut.LIST_MV_CLASSES
        self.str_id = str_id
        self.str_quant = str_quant

    @ut.doc_params(doc_param_df_groups_upmnar=doc_param_df_groups_upmnar)
    def get_limits(self,
                   df: pd.DataFrame = None,
                   groups: ut.ArrayLike1D = None,
                   loc_pct_upmnar: float = 0.25,
                   cols_quant: ut.ArrayLike1D = None
                   ) -> Tuple[float, float, float]:
        """
        Get minimum of detected values (d_min, i.e., detection limit), upper bound of MNAR MVs (up_mnar),
        and maximum of detected values (d_max).

        Parameters
        ----------
        {doc_param_df_groups_upmnar}
        cols_quant
            Column names with quantification data in ``df``.

        Return
        ------
        d_min
            Minimum of detected values
        up_mnar
            upper bound of MNAR MVs
        d_max
            Maximum of detected values
        """
        # Check input
        cols_quant = ut.check_list_like(name="cols_quant", val=cols_quant, accept_none=True)
        df = ut.check_df(df=df, accept_none=False, cols_req=cols_quant)
        groups = ut.check_list_like(name="groups", val=groups, accept_none=False)
        ut.check_match_df_groups(groups=groups, df=df, str_quant=self.str_quant)
        ut.check_number_range(name="loc_pct_upmnar", val=loc_pct_upmnar, min_val=0, max_val=1,
                              just_int=False, accept_none=False)
        # Compute limits
        if cols_quant is None:
            cols_quant = ut.get_qcols(df=df, groups=groups, str_quant=self.str_quant)
        d_min, up_mnar = get_up_mnar(df=df[cols_quant], loc_pct_upmnar=loc_pct_upmnar)
        d_max = df[cols_quant].max().max()
        return d_min, up_mnar, d_max

    @ut.doc_params(doc_param_df_groups_upmnar=doc_param_df_groups_upmnar)
    def run(self,
            df: pd.DataFrame = None,
            groups: ut.ArrayLike1D = None,
            loc_pct_upmnar: float = 0.25,
            min_cs: float = 0.5,
            n_neighbors: int = 5
            ) -> pd.DataFrame:
        """
        Run cImpute algorithm.

        Hybrid method for imputation of omics data called conditional imputation (cImpute)
        using MinProb for MNAR (Missing Not at Random) missing values and KNN imputation for
        MCAR (Missing completely at Random) missing values.

        Parameters
        ----------
        {doc_param_df_groups_upmnar}
        min_cs
            Minimum of confidence score [0-1] used for selecting values for protein in groups to apply imputation on.
        n_neighbors
            Number of neighboring samples to use for MCAR imputation by KNN.

        Return
        ------
        df_imp
            DataFrame with (a) imputed intensities values and (b) group-wise confidence score and NaN classification.

        Notes
        -----
        - MAR is only imputed if ``min_cs=0`` using the imputation for MCAR.
        """
        # Check input
        df = ut.check_df(df=df, accept_none=False)
        groups = ut.check_list_like(name="groups", val=groups, accept_none=False)
        ut.check_match_df_groups(groups=groups, df=df, str_quant=self.str_quant)
        ut.check_number_range(name="loc_pct_upmnar", val=loc_pct_upmnar, min_val=0, max_val=1,
                              just_int=False, accept_none=False)
        ut.check_number_range(name="min_cs", val=min_cs, min_val=0, max_val=1,
                              just_int=False, accept_none=False)
        ut.check_number_range(name="n_neighbors", val=n_neighbors, min_val=1,
                              just_int=True, accept_none=False)
        # Run imputation
        df_imp = run_cimpute(df=df, groups=groups,
                             min_cs=min_cs, loc_pcat_upmnar=loc_pct_upmnar, n_neighbors=n_neighbors,
                             str_quant=self.str_quant, str_id=self.str_id)
        return df_imp



