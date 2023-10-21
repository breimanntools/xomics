"""
This is a script for interface of the pRank (protein-centric ranking) class.
"""
import pandas as pd
import numpy as np
from typing import Optional, List

import xomics.utils as ut
from ._backend.prank import p_score, e_score, c_score, e_score_only_pvals
from ._backend.ehits import e_hits


# I Helper Functions
# TODO assess whether check functions can be replaced by ut.check functions.
# TODO refactor check functions
# Check functions
def check_all_non_negative(val_list, name=None):
    """Check if all values in the list are non-negative numbers."""
    if not all(isinstance(x, (int, float)) for x in val_list):
        raise ValueError(f"All elements in '{name}' should be integers or floats.")

    if any(x < 0 for x in val_list):
        raise ValueError(f"All elements in '{name}' should be non-negative.")


def check_input_scoring_match(*args):
    """Check if all input lists/arrays have the same length."""
    # Check whether input lists of same length
    list_n = [len(x) for x in args]
    str_error = f"All input lists must have the same number of items."
    if max(list_n) != min(list_n):
        raise ValueError(str_error)
    f = lambda x: x if isinstance(x, np.ndarray) or ut.is_nested_list_with_strings(x) else np.array(x)
    return [f(x) for x in args]


def check_numeric_elements(x, name=None):
    """"""
    if isinstance(x, list):
        if not not all(isinstance(x, (int, float)) for x in x):
            raise ValueError(f"'{name}' should only contain numerical values")
    elif isinstance(x, np.ndarray):
        if not np.all(np.isreal(x)):
            raise ValueError(f"'{name}' should only contain numerical values")
    else:
        raise ValueError(f"'{name}' should be list or numpy.array")


def _check_all_strings(name=None, lst=None):
    """"""
    non_string_elements = [x for x in lst if not isinstance(x, str)]
    if len(non_string_elements) > 0:
        raise ValueError(f"All elements in '{name}' should not be string. Following are not: {non_string_elements}")


def _check_duplicates(name=None, lst=None):
    """"""
    seen = set()
    list_duplicates = [item for item in lst if item in seen or seen.add(item)]
    if len(lst) != len(set(lst)):
        raise ValueError(f"{name} contains the following duplicated elements: {list_duplicates}")


def _check_emtpy_input(ids=None, id_lists=None, list_terms=None):
    """"""
    if ids is None or len(ids) == 0:
        raise ValueError(f"'ids' should not be empty.")
    if id_lists is None or len(id_lists) == 0:
        raise ValueError(f"'id_lists' should not be empty.")
    if list_terms is None or len(list_terms) == 0:
        raise ValueError(f"'term_list' should not be empty.")


def _check_id_lists_term_list_lengths(id_lists, term_list):
    """Check if the length of 'id_lists' is equal to the length of 'term_list'."""
    if len(id_lists) != len(term_list):
        raise ValueError(f"The length of 'id_lists' ({len(id_lists)}) must be equal to the length of 'term_list' ({len(term_list)}).")


def _check_id_lists_is_nested_list(id_lists):
    """Check if 'id_lists' is a nested list of lists."""
    if isinstance(id_lists, pd.Series):
        id_lists = [[ids] for ids in id_lists]
    if not ut.is_nested_list_with_strings(id_lists):
        raise ValueError("'id_lists' should be a nested list of lists or pd.Series.")
    else:
        return id_lists


def _check_list_input(name=None, lst=None):
    """"""
    if ut.is_nested_list_with_strings(lst):
        raise ValueError(f"'{name}' should not be nested")
    if isinstance(lst, pd.Series):
        return lst.to_list()
    elif isinstance(lst, list):
        return lst
    else:
        raise ValueError(f"'{name}' should be list or pd.Series")


def _check_terms_sub_list(name=None, terms_sub_list=None, terms=None):
    """"""
    if terms_sub_list is None:
        return terms_sub_list
    if ut.is_nested_list_with_strings(terms_sub_list):
        raise ValueError(f"'{name}' should not be nested")
    if isinstance(terms_sub_list, pd.Series):
        terms_sub_list = terms_sub_list.to_list()
    elif isinstance(terms_sub_list, list):
        pass
    else:
        raise ValueError(f"'{name}' should be list or pd.Series")
    not_in_terms = [x for x in terms_sub_list if x not in terms]
    if len(not_in_terms) > 0:
        raise ValueError(f"Following elements of '{name}' are not in 'terms': {not_in_terms}")
    _check_all_strings(name="terms_sub_list", lst=terms_sub_list)
    _check_duplicates(name="terms_sub_list", lst=terms_sub_list)
    return terms_sub_list


# TODO finsih docuemtnation, typing, refactor check functions (include log check), testing, tutorial
# II Main Functions
class pRank:
    """
    Hybrid imputation algorithm for missing values (MVs) in (prote)omics data.
    """
    def __init__(self,
                 col_id: str = ut.COL_PROT_ID,
                 col_name: str = ut.COL_GENE_NAME,
                 str_quant: str = ut.STR_QUANT,
                 ):
        """
        Parameters
        ----------
        col_id
            Name of column with identifiers in DataFrame.
        col_name
            Name of column with sample names in DataFrame.
        str_quant
            Identifier for the LFQ columns in the DataFrame.
        """
        ut.check_str(name="col_id", val=col_id, accept_none=False)
        ut.check_str(name="col_name", val=col_name, accept_none=False)
        ut.check_str(name="str_quant", val=str_quant, accept_none=False)
        self.col_id = col_id
        self.col_name = col_name
        self.str_quant = str_quant

    @staticmethod
    def p_score(df_fc: pd.DataFrame = None,
                col_fc: str = None,
                col_pval: str = None,
                ) -> pd.DataFrame:
        """
        Calculate the single protein use_cases ranking score (P score) by first z-normalizing fold change scores
        and p-values, and then integrating them protein-wise to obtain min-max normalized ranking scores.

        Parameters
        ----------
        df_fc
            DataFrame with fold-change and p-values.
        col_fc
            Name of column from ``df`` with fold change values for each protein (log2 fold recommended).
        col_pval
            Name of column from ``df`` with p-values for each protein (-log10 fold recommended).

        Returns
        -------
        df_fc
            Input DataFrame with p-score for each protein given in 'P-Score' column.
        """
        # Checking functions
        df_fc = ut.check_df(name="df_fc", df=df_fc)
        ut.check_col_in_df(df=df_fc, name_df="df_fc", cols=[col_fc, col_pval], name_cols=["col_fc", "col_pval"])
        # Get arrays with values
        x_fc = df_fc[col_fc].values
        x_pval = df_fc[col_pval].values
        # Check columns values
        x_fc, x_pval = check_input_scoring_match(x_fc, x_pval)
        check_numeric_elements(x_fc, name="x_fc")
        check_numeric_elements(x_pval, name="x_pvals")
        # Get P-score
        p_scores = p_score(x_fc=x_fc, x_pvals=x_pval)
        df_fc[ut.COL_P_SCORE] = p_scores
        return df_fc

    @staticmethod
    def e_score(df_fc: pd.DataFrame = None,
                col_name: str = None,
                df_enrich: pd.DataFrame = None,
                col_fe: str = None,
                col_pval: str = None,
                col_name_lists: str = None,
                ) -> pd.DataFrame:
        """
        Calculate the single protein enrichment score (E score) by first z-normalizing fold enrichment scores and
        p-values, and then integrating them protein-wise to obtain a min-max normalized ranking score.

        Parameters
        ----------
        df_fc
            DataFrame with fold-change and p-values.
        col_name
            Name of column from ``df_fc`` with protein names.
        df_enrich
            DataFrame with fold enrichment and p-values for each enrichment term.
        col_fe
            Name of column from ``df_enrich`` with fold enrichment values for each enrichment term (log2 fold recommended).
        col_pval
            Name of column from ``df_enrich`` with p-values for each term (-log10 fold recommended).
        col_name_lists
            Name of column from ``df_enrich`` with protein name lists. Lists should contain names from ``col_names``.

        Returns
        -------
        df_fc
            Input DataFrame with E-score for each protein given in 'P-Score' column.
        """
        # TODO check for duplicated Term
        # Checking functions
        df_fc = ut.check_df(name="df_fc", df=df_fc)
        ut.check_str(name="col_name", val=col_name)
        ut.check_col_in_df(df=df_fc, name_df="df_fc", cols=col_name, name_cols="col_id")
        df_enrich = ut.check_df(name="df_enrich", df=df_enrich)
        ut.check_str(name="col_pval", val=col_pval)
        ut.check_str(name="col_fe", val=col_fe, accept_none=True)
        ut.check_str(name="col_name_lists", val=col_name_lists)
        ut.check_col_in_df(df=df_enrich, name_df="df_enrich", cols=[col_pval, col_name_lists], name_cols=["col_pval", "col_id_lists"])
        # Get arrays with values
        names = df_fc[col_name].values
        x_pval = df_enrich[col_pval].values
        name_lists = df_enrich[col_name_lists].values
        # Check columns values
        if col_fe is not None:
            ut.check_col_in_df(df=df_enrich, name_df="df_enrich", cols=[col_fe], name_cols=["col_fe"])
            x_fe = df_enrich[col_fe].values
            _, x_fe, x_pval = check_input_scoring_match(name_lists, x_fe, x_pval)
            check_numeric_elements(x_fe, name="col_fe")
        else:
            _, x_pval = check_input_scoring_match(name_lists, x_pval)
            check_numeric_elements(x_pval, name="col_pvals")
        # Get E-score
        if col_fe is not None:
            e_scores = e_score(names=names, name_lists=name_lists, x_fe=x_fe, x_pval=x_pval)
        else:
            e_scores = e_score_only_pvals(names=names, name_lists=name_lists, x_pval=x_pval)
        df_fc[ut.COL_E_SCORE] = e_scores
        return df_fc

    @staticmethod
    def c_score(df_imp: pd.DataFrame = None,
                ids: ut.ArrayLike1D = None,
                col_id: str = None
                ) -> pd.DataFrame:
        """Obtain protein use_cases confidence score (C score) from cImpute output

        Parameters
        ----------
        df_imp
            Data Frame with imputed values from ``cImpute``.
        ids
            List or array of protein identifiers.
        col_id
            Name of id column from 'df_imp'. If None, index will be considered for ids.

        Returns
        -------
        c_scores
            Array of confidence scores (C scores) from imputation for each protein.
        """
        c_scores = c_score(ids=ids, df_imp=df_imp, col_id=col_id)
        df_imp[ut.COL_C_SCORE] = c_scores
        return df_imp

    @staticmethod
    def e_hits(ids=None,
               id_lists=None,
               terms=None,
               terms_sub_list=None,
               n_ids=None,
               n_terms=None,
               sort_alpha=False
               ) -> pd.DataFrame:
        """
        Get association matrix for protein ids and enrichment terms.

        Get matrix with associations between protein/gene ids and id sets representing protein/gene lists
        associated with specific biological terms obtained from an enrichment analysis (referred to as 'enrichment terms')
        such as GO or KEGG pathway terms.

        Parameters
        ----------
        ids : array-like
            Array of protein identifiers.
        id_lists : list of lists
            List of protein identifier sets from enrichment analysis (e.g., set of proteins linked to specific GO term).
        terms : list or array-like
            List of enrichment terms matching to id_lists
        terms_sub_list : list or array-like, default = None
            Sublist of enrichment terms (must be subset of 'terms'). If not None, terms will be used to filter output
        n_ids : integer, default = None
            Filter results for 'n_ids' genes/proteins from 'ids' with the highest number of associations if not None
        n_terms : integer, default = None
            Filter results for 'n_terms' terms from 'term_list' with the highest number of associations if not None
        sort_alpha : bool, default = False
            Sort falues in alphabetically (if True) or in descending order of hit counts (if False)

        Returns
        -------
        df_e_hit : pandas.DataFrame
            Data frame with links between gene/protein ids and 'enrichment' terms

        Examples
        --------
        >>> e_hits(ids=['gene1', 'gene2', 'gene3'], id_lists=[['gene1', 'gene2'], ['gene2', 'gene3']],
        ... terms=['term1', 'term2'], n_ids=2, n_terms=1)
               gene1  gene2
        term1      1      1
        """
        # Check input
        _check_emtpy_input(ids=ids, id_lists=id_lists, list_terms=terms)
        id_lists = _check_id_lists_is_nested_list(id_lists)
        ids = _check_list_input(name="ids", lst=ids)
        terms = _check_list_input(name="list_terms", lst=terms)
        terms_sub_list = _check_terms_sub_list(name="terms_sub_list", terms_sub_list=terms_sub_list, terms=terms)
        _check_all_strings(name="ids", lst=ids)
        _check_duplicates(name="ids", lst=ids)
        _check_all_strings(name="terms", lst=terms)
        _check_duplicates(name="terms", lst=terms)
        ut.check_number_range(name="n_ids", val=n_ids, min_val=1, accept_none=True, just_int=True)
        ut.check_number_range(name="n_terms", val=n_terms, min_val=1, accept_none=True, just_int=True)
        # Obtain gene/protein associations with enrichment terms
        df_e_hits = e_hits(ids=ids, id_lists=id_lists, terms=terms, terms_sub_list=terms_sub_list, n_ids=n_ids,
                           n_terms=n_terms, sort_alpha=sort_alpha)
        return df_e_hits

