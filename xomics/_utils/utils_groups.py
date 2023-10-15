"""
This is a script for functions to access group information from a provided dataframe.
"""


def get_dict_qcol_group(df=None, groups=None, str_quant=None):
    """
    Create a dictionary with groups from df based on lfq_str and given groups

    Parameters
    ----------
    df: DataFrame containing quantified features including missing values
    groups: List with group names
    str_quant: Substring in column indicating quantification
    """
    if str_quant is None:
        raise ValueError("'str_quant' must be given")
    dict_col_group = {}
    for col in list(df):
        if str_quant in col:
            col_wo_str_quant = col.replace(str_quant, "")
            for group in groups:
                if group in col_wo_str_quant:
                    dict_col_group[col] = group
    return dict_col_group


def get_dict_group_qcols(df=None, groups=None, str_quant=None):
    """
    Create a dictionary with groups from df based on lfq_str and given groups

    Parameters
    ----------
    df: DataFrame containing quantified features including missing values
    groups: List with group names
    str_quant: Substring in column indicating quantification
    """
    if str_quant is None:
        raise ValueError("'str_quant' must be given")
    dict_col_group = get_dict_qcol_group(df=df, groups=groups, str_quant=str_quant)
    dict_group_cols = {g: [k for k, v in dict_col_group.items() if v == g] for g in groups}
    return dict_group_cols


def get_qcols(df=None, groups=None, str_quant=None):
    """
    Create a list with groups from df based on lfq_str and given groups

    Parameters
    ----------
    df: DataFrame containing quantified features including missing values
    groups: List with group names
    str_quant: Substring in column indicating quantification
    """
    if str_quant is None:
        raise ValueError("'str_quant' must be given")
    dict_col_group = get_dict_qcol_group(df=df, groups=groups, str_quant=str_quant)
    dict_group_cols = {g: [k for k, v in dict_col_group.items() if v == g] for g in groups}
    list_group_cols = [col for group_cols in dict_group_cols.values() for col in group_cols]
    return list_group_cols
