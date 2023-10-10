"""
This is a script for an enrichment association hit analysis
"""
import pandas as pd
import xomics._utils as ut


# I Helper Functions
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
def e_hits(ids=None, id_lists=None, terms=None, terms_sub_list=None, n_ids=None, n_terms=None, sort_alpha=False):
    """
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
    list_hits = []
    for ids_in_term in id_lists:
        _ids_in_term = ut.flatten_list(ids_in_term)
        list_hits.append([int(x in _ids_in_term) for x in ids])
    df_e_hits = pd.DataFrame(list_hits, columns=ids, index=terms)
    # Filter results
    if terms_sub_list is not None:
        if sort_alpha:
            terms_sub_list = sorted(terms_sub_list)
        df_e_hits = df_e_hits.loc[terms_sub_list]
    if n_ids is not None:
        # Filter the DataFrame by keeping only the 'n_ids' with the highest number of associations
        id_sums = df_e_hits.sum(axis=0).sort_values(ascending=False)
        top_ids = id_sums.nlargest(n_ids).index
        if sort_alpha:
            top_ids = sorted(top_ids)
        df_e_hits = df_e_hits[top_ids]
    if n_terms is not None and terms_sub_list is None:
        # Filter the DataFrame by keeping only the 'n_terms' with the highest number of associations
        term_sums = df_e_hits.sum(axis=1).sort_values(ascending=False)
        top_terms = term_sums.nlargest(n_terms).index
        if sort_alpha:
            top_terms = sorted(top_terms)
        df_e_hits = df_e_hits.loc[top_terms]
    return df_e_hits
