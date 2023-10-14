"""
This is a script for the backend of the e_hits (enrichment association hit) function of the pRank object.
"""
import pandas as pd
import xomics.utils as ut


# I Helper Functions


# II Main Functions
def e_hits(ids=None, id_lists=None, terms=None, terms_sub_list=None, n_ids=None, n_terms=None, sort_alpha=False):
    """
    Get matrix with associations between protein/gene ids and id sets representing protein/gene lists
    associated with specific biological terms obtained from an enrichment analysis (referred to as 'enrichment terms')
    such as GO or KEGG pathway terms.
    """
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
