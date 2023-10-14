"""
This is a script for interface of the pRank (protein-centric ranking) class.
"""
import pandas as pd
import numpy as np

import xomics.utils as ut
from ._backend.prank import p_score, e_score, c_score
from ._backend.ehits import e_hits


# I Helper Functions
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


def _check_non_negative_number(name=None, val=None, min_val=0, max_val=None, accept_none=False, just_int=True):
    """Check if value of given name variable is non-negative integer"""
    check_types = [int] if just_int else [float, int]
    str_check = "non-negative int" if just_int else "non-negative float or int"
    add_str = f"n>{min_val}" if max_val is None else f"{min_val}<=n<={max_val}"
    if accept_none:
        add_str += " or None"
    error = f"'{name}' ({val}) should be {str_check} n, where " + add_str
    if accept_none and val is None:
        return None
    if type(val) not in check_types:
        raise ValueError(error)
    if val < min_val:
        raise ValueError(error)
    if max_val is not None and val > max_val:
        raise ValueError(error)


# II Main Functions
class pRank:
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

    @staticmethod
    def p_score(ids=None, x_fc=None, x_pvals=None, adjust_log=True, verbose=True):
        """
        Calculate the single protein proteomics ranking score (P score) by first z-normalizing fold change scores
        and p-values, and then integrating them protein-wise to obtain min-max normalized ranking scores.

        Parameters
        ----------
        ids : list or array-like
            List or array of protein identifiers.
        x_fc : array-like
            Array of fold changes for each protein (log2 fold recommanded).
        x_pvals : array-like
            Array of p-values for each protein (log10 fold recommanded).
        adjust_log : bool, default=True
            Whether to automatically test and conduct
            - log2 transformormation for fold enrichment values
            - -log10 transformation for p-values.
        verbose : bool, default = True
            Verbosity mode on or off

        Returns
        -------
        p_scores : numpy.ndarray
            Array of proteomics ranking scores (P scores) for each protein.

        Examples
        --------
        >>> p_score(ids=['protein1', 'protein2'], x_fc=[2.4, 1.5], x_pvals=[0.05, 0.2])
        array([1.0, 0.])
        """
        # Checking functions
        ids, x_fc, x_pvals = check_input_scoring_match(ids, x_fc, x_pvals)
        check_numeric_elements(x_fc, name="x_fc")
        check_numeric_elements(x_pvals, name="x_pvals")
        ut.check_bool(name="verbose", val=verbose)
        # Get P-score
        p_scores = p_score(ids=ids, x_fc=x_fc, x_pvals=x_pvals, adjust_log=adjust_log, verbose=verbose)
        return p_scores

    @staticmethod
    def e_score(ids=None, id_lists=None, x_fe=None, x_pvals=None, adjust_log=True, verbose=True):
        """
        Calculate the single protein enrichment score (E score) by first z-normalizing fold enrichment scores and
        p-values, and then integrating them protein-wise to obtain a min-max normalized ranking score.

        Parameters
        ----------
        ids : list or array-like
            List or array of protein identifiers.
        id_lists : list of lists
            List of protein identifier sets from enrichment analysis (e.g., set of proteins linked to specific GO term)
        x_fe : array-like
            Array of fold enrichments for each protein set.
        x_pvals : array-like
            Array of p-values for each protein set.
        adjust_log : bool, default=True
            Whether to automatically test and conduct.
            - log2 transformormation for fold enrichment values
            - -log10 transformation for p-values.
        verbose : bool, default = True
            Verbosity mode on or off

        Returns
        -------
        e_scores : numpy.ndarray
            Array of enrichment ranking scores (E scores) for each protein.

        Examples
        --------
        >>> e_score(ids=['protein1', 'protein2'], id_lists=[['protein1', 'protein2'], ['protein2']],
        ... x_fe=[2, 1.5], x_pvals=[0.05, 0.1])
        [0., 1.0]
        """
        # Checking functions
        # TODO check for duplicated Term
        _, x_fe, x_pvals = check_input_scoring_match(id_lists, x_fe, x_pvals)
        check_numeric_elements(x_fe, name="x_fe")
        check_numeric_elements(x_pvals, name="x_pvals")
        check_all_non_negative(x_fe, name="x_fe")
        ut.check_bool(name="verbose", val=verbose)
        # Get E-score
        e_scores = e_score(ids=ids, id_lists=id_lists, x_fe=x_fe, x_pvals=x_pvals,
                           adjust_log=adjust_log, verbose=verbose)
        return e_scores

    @staticmethod
    def c_score(ids=None, df_imp=None, col_id=None):
        """Obtain protein proteomics confidence score (C score) from cImpute output

        Parameters
        ----------
        ids : list or array-like
            List or array of protein identifiers.
        df_imp : pandas.DataFrame
            Data Frame from cImpute.
        col_id : str, default=None
            Name of id column from 'df_imp'. If None, index will be considered for ids.

        Returns
        -------
        c_scores : numpy.ndarray
            Array of confidence scores (C scores) from imputation for each protein.
        """
        c_scores = c_score(ids=ids, df_imp=df_imp, col_id=col_id)
        return c_scores

    @staticmethod
    def e_hits(ids=None, id_lists=None, terms=None, terms_sub_list=None, n_ids=None, n_terms=None, sort_alpha=False):
        """Get association matrix for protein ids and enrichment terms.

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
        _check_non_negative_number(name="n_ids", val=n_ids, min_val=1, accept_none=True, just_int=True)
        _check_non_negative_number(name="n_terms", val=n_terms, min_val=1, accept_none=True, just_int=True)
        # Obtain gene/protein associations with enrichment terms
        df_e_hits = e_hits(ids=ids, id_lists=id_lists, terms=terms, terms_sub_list=terms_sub_list, n_ids=n_ids,
                           n_terms=n_terms, sort_alpha=sort_alpha)
        return df_e_hits

