"""
This is a script for data checking utility functions.
"""
import pandas as pd
import numpy as np
from sklearn.utils import check_array

from ._utils import add_str
from .utils_types import VALID_INT_TYPES, VALID_INT_FLOAT_TYPES
import xomics._utils.check_type as check_type


# Helper functions
def _convert_2d(val=None, name=None, str_add=None):
    """
    Convert array-like data to 2D array. Handles lists of arrays, lists of lists, and 1D lists.
    """
    str_error = add_str(str_error=f"'{name}' should be a 2D list or 2D array with rows having the same number of columns.",
                        str_add=str_add)
    if isinstance(val, list):
        # Check if List with arrays and return if yes
        if all(isinstance(i, np.ndarray) for i in val):
            try:
                val = np.asarray(val)
            except ValueError:
                raise ValueError(str_error)
        # Convert 1D list to 2D list
        elif all(not isinstance(i, list) for i in val):
            try:
                val = np.asarray([val])
            except ValueError:
                raise ValueError(str_error)
        # For nested lists, ensure they are 2D (list of lists with equal lengths)
        else:
            try:
                val = np.array(val)  # Convert nested list to numpy array
                if val.ndim != 2:
                    raise ValueError
            except ValueError:
                raise ValueError(str_error)
    elif hasattr(val, 'ndim') and val.ndim == 1:
        try:
            val = np.asarray([val])
        except ValueError:
            raise ValueError(str_error)
    return val


# Check array like
def check_array_like(name=None, val=None, dtype=None, ensure_2d=False, allow_nan=False, convert_2d=False,
                     accept_none=False, expected_dim=None, str_add=None):
    """Check if the provided value is array-like and matches the specified dtype."""
    if val is None:
        if accept_none:
            return None  # Skip tests
        else:
            raise ValueError(f"'{name}' should not be None.")
    # Extend dtype to handle a list of dtypes including bool
    dtype = check_type.check_str(name="dtype", val=dtype, accept_none=True)
    valid_dtypes = ["numeric", "int", "float", "bool", None]
    if dtype not in valid_dtypes:
        str_error = add_str(str_error=f"'dtype' should be one of the following: {valid_dtypes}", str_add=str_add)
        raise ValueError(str_error)
    dict_expected_dtype = {"numeric": "numeric", "int": "int64", "float": "float64", "bool": "bool"}
    expected_dtype = dict_expected_dtype[dtype] if dtype is not None else None
    # Specific check for boolean arrays
    if dtype == "bool":
        flattened_val = np.array(val).flatten()
        if not all(isinstance(item, (bool, np.bool_)) for item in flattened_val):
            str_error = add_str(str_error=f"All elements in '{name}' must be of type 'bool' (either Python native or NumPy bool).",
                                str_add=str_add)
            raise ValueError(str_error)
    # Convert a 1D list or array to a 2D array if needed
    if convert_2d:
        val = _convert_2d(val=val, name=name, str_add=str_add)
    # Utilize Scikit-learn's check_array for robust checking
    try:
        # Convert list to array
        val = check_array(val, dtype=expected_dtype, ensure_2d=ensure_2d, force_all_finite=not allow_nan)
    except Exception as e:
        dtype = "any type" if dtype is None else dtype
        raise ValueError(f"'{name}' should be array-like with '{dtype}' values."
                         f"\nScikit message:\n\t{e}")
    # Check dimensions if specified
    if expected_dim is not None and len(val.shape) != expected_dim:
        str_error = add_str(str_error=f"'{name}' should have {expected_dim} dimensions, but has {len(val.shape)}.",
                            str_add=str_add)
        raise ValueError(str_error)
    return val


# Check feature matrix and labels
def check_X(X, X_name="X", min_n_samples=3, min_n_features=2, min_n_unique_features=None,
            ensure_2d=True, allow_nan=False, accept_none=False, str_add=None):
    """Check the feature matrix X is valid."""
    if X is None:
        if not accept_none:
            raise ValueError(f"'{X_name}' should not be None")
        else:
            return None
    X = check_array_like(name="X", val=X, dtype="float", ensure_2d=ensure_2d, allow_nan=allow_nan)
    if np.isinf(X).any():
        str_error = add_str(str_error=f"'X' should not contain infinite values",
                            str_add=str_add)
        raise ValueError(str_error)
    n_samples, n_features = X.shape
    if n_samples < min_n_samples:
        raise ValueError(f"n_samples ({n_samples} in 'X') should be >= {min_n_samples}."
                         f"\nX = {X}")
    if n_features < min_n_features:
        raise ValueError(f"n_features ({n_features} in 'X') should be >= {min_n_features}."
                         f"\nX = {X}")
    if min_n_unique_features is not None:
        n_unique_features = sum([len(set(X[:, col])) > 1 for col in range(n_features)])
        if n_unique_features < min_n_unique_features:
            str_error = add_str(str_error=f"'n_unique_features' ({n_unique_features}) should be >= {min_n_unique_features}",
                                str_add=str_add)
            raise ValueError(str_error)
    return X


def check_X_unique_samples(X, min_n_unique_samples=3, str_add=None):
    """Check if the matrix X has a sufficient number of unique samples."""
    n_unique_samples = len(set(map(tuple, X)))
    if n_unique_samples == 1:
        str_error = add_str(str_error="Feature matrix 'X' should not have all identical samples.",
                            str_add=str_add)
        raise ValueError(str_error)
    if n_unique_samples < min_n_unique_samples:
        raise ValueError(f"n_unique_samples ({n_unique_samples}) should be >= {min_n_unique_samples}."
                         f"\nX = {X}")
    return X


def check_labels(labels=None, name="labels", vals_requiered=None, len_requiered=None, allow_other_vals=True,
                 n_per_group_requiered=None, accept_float=False, str_add=None):
    """Check the provided labels against various criteria like type, required values, and length."""
    if labels is None:
        raise ValueError(f"'{name}' should not be None.")
    labels = check_type.check_list_like(name=name, val=labels)
    # Convert labels to a numpy array if it's not already
    labels = np.asarray(labels)
    # Ensure labels is at least 1-dimensional
    if labels.ndim == 0:
        labels = np.array([labels.item()])  # Convert 0-d array to 1-d array
    unique_labels = set(labels)
    if len(unique_labels) == 1:
        str_error = add_str(str_error=f"'{name}' should contain more than one different value ({unique_labels}).",
                            str_add=str_add)
        raise ValueError(str_error)
    valid_types = VALID_INT_TYPES if not accept_float else VALID_INT_FLOAT_TYPES
    wrong_types = [l for l in unique_labels if not isinstance(l, valid_types)]
    if wrong_types:
        str_error = add_str(str_error=f"Labels in '{name}' should be type int, but contain: {set(map(type, wrong_types))}",
                            str_add=str_add)
        raise ValueError(str_error)
    if vals_requiered is not None:
        missing_vals = [x for x in vals_requiered if x not in labels]
        if len(missing_vals) > 0:
            str_error = add_str(str_error=f"'{name}' ({unique_labels}) does not contain requiered values: {missing_vals}",
                                str_add=str_add)
            raise ValueError(str_error)
        if not allow_other_vals:
            wrong_vals = [x for x in labels if x not in vals_requiered]
            if len(wrong_vals) > 0:
                str_error = add_str(str_error=f"'{name}' ({unique_labels}) does contain wrong values: {wrong_vals}",
                                    str_add=str_add)
                raise ValueError(str_error)
    if len_requiered is not None and len(labels) != len_requiered:
        str_error = add_str(str_error=f"'{name}' (n={len(labels)}) should contain {len_requiered} values.",
                            str_add=str_add)
        raise ValueError(str_error)
    # Check for minimum length per group
    if n_per_group_requiered is not None:
        label_counts = {label: np.sum(labels == label) for label in unique_labels}
        underrepresented_labels = {label: count for label, count in label_counts.items() if
                                   count < n_per_group_requiered}
        if underrepresented_labels:
            str_error = add_str(str_error=f"Each label should have at least {n_per_group_requiered} occurrences. "
                                          f"Underrepresented labels: {underrepresented_labels}",
                                str_add=str_add)
            raise ValueError(str_error)
    return labels


def check_match_X_labels(X=None, X_name="X", labels=None, labels_name="labels", check_variability=False,
                         str_add=None):
    """Check if the number of samples in X matches the number of labels."""
    n_samples, n_features = X.shape
    if n_samples != len(labels):
        str_error = add_str(str_error=f"n_samples does not match for '{X_name}' ({len(X)}) and '{labels_name}' ({len(labels)}).",
                            str_add=str_add)
        raise ValueError(str_error)
    if check_variability:
        unique_labels = np.unique(labels)
        for label in unique_labels:
            group_X = X[labels == label]
            if not np.all(np.var(group_X, axis=0) != 0):
                str_error = add_str(str_error=f"Variance in 'X' for label '{label}' from '{labels_name}' is too low.",
                                    str_add=str_add)
                raise ValueError(str_error)


def check_match_X_list_labels(X=None, list_labels=None, check_variability=False, vals_requiered=None, str_add=None):
    """Check if each label set is matching with X"""
    for i, labels in enumerate(list_labels):
        check_labels(labels=labels, vals_requiered=vals_requiered)
        check_match_X_labels(X=X, labels=labels, labels_name=f"list_labels (set {i+1})",
                             check_variability=check_variability, str_add=str_add)


def check_match_list_labels_names_datasets(list_labels=None, names_datasets=None, str_add=None):
    """Check if length of list_labels and names match"""
    if names_datasets is None:
        return None # Skip check
    if len(list_labels) != len(names_datasets):
        str_error = add_str(str_error=f"Length of 'list_labels' ({len(list_labels)}) and 'names_datasets'"
                                      f" ({len(names_datasets)} does not match)",
                            str_add=str_add)
        raise ValueError(str_error)


# Check sets
def check_superset_subset(subset=None, superset=None, name_subset=None, name_superset=None, str_add=None):
    """Check if all elements of the subset are contained in the superset."""
    wrong_elements = [x for x in subset if x not in superset]
    if len(wrong_elements) != 0:
        str_error = add_str(str_error=f"'{name_superset}' does not contain the following elements of '{name_subset}': {wrong_elements}",
                            str_add=str_add)
        raise ValueError(str_error)


# df checking functions
def check_df(name="df", df=None, accept_none=False, accept_nan=True, check_all_positive=False,
             cols_requiered=None, cols_forbidden=None, cols_nan_check=None, str_add=None):
    """Check if the provided DataFrame meets various criteria such as NaN values, required/forbidden columns, etc."""
    # Check DataFrame and values
    if df is None:
        if not accept_none:
            raise ValueError(f"'{name}' should not be None")
        else:
            return None
    if not isinstance(df, pd.DataFrame):
        str_error = add_str(str_error= f"'{name}' ({type(df)}) should be DataFrame",
                            str_add=str_add)
        raise ValueError(str_error)
    if not accept_nan and df.isna().any().any():
        str_error = add_str(str_error=f"'{name}' contains NaN values, which are not allowed",
                            str_add=str_add)
        raise ValueError(str_error)
    if check_all_positive:
        numeric_df = df.select_dtypes(include=['float', 'int'])
        if numeric_df.min().min() <= 0:
            str_error = add_str(str_error=f"'{name}' should not contain non-positive values.",
                                str_add=str_add)
            raise ValueError(str_error)
    # Check columns
    args = dict(accept_str=True, accept_none=True, str_add=str_add)
    cols_requiered = check_type.check_list_like(name='cols_requiered', val=cols_requiered, **args)
    cols_forbidden = check_type.check_list_like(name='cols_forbidden', val=cols_forbidden, **args)
    cols_nan_check = check_type.check_list_like(name='cols_nan_check', val=cols_nan_check, **args)
    if cols_requiered is not None:
        missing_cols = [col for col in cols_requiered if col not in df.columns]
        if len(missing_cols) > 0:
            str_error = add_str(str_error=f"'{name}' is missing required columns: {missing_cols}",
                                str_add=str_add)
            raise ValueError(str_error)
    if cols_forbidden is not None:
        forbidden_cols = [col for col in cols_forbidden if col in df.columns]
        if len(forbidden_cols) > 0:
            str_error = add_str(str_error=f"'{name}' is contains forbidden columns: {forbidden_cols}",
                                str_add=str_add)
            raise ValueError(str_error)
    if cols_nan_check is not None:
        if df[cols_nan_check].isna().sum().sum() > 0:
            str_error = add_str(str_error=f"NaN values are not allowed in '{cols_nan_check}'.",
                                str_add=str_add)
            raise ValueError(str_error)
    return df


def check_col_in_df(df=None, name_df=None, cols=None, name_cols=None, accept_nan=False, error_if_exists=False, accept_none=False):
    """
    Check if the column or columns exists in the DataFrame, if the values have the correct type, and if NaNs are allowed.
    """
    # Check if the column already exists and raise error if error_if_exists is True
    if error_if_exists and (cols in df.columns):
        raise ValueError(f"Column '{cols}' already exists in '{name_df}'")
    if cols is None:
        if not accept_none:
            raise ValueError(f"'{name_cols}' should not be None.")
        return None
    # Check if the column exists in the DataFrame
    if isinstance(cols, str):
        if cols not in df.columns:
            raise ValueError(f"'{cols}' must be a column in '{name_df}': {list(df.columns)}")
    else:
        wrong_columns = [c for c in cols if c not in df.columns]
        if len(wrong_columns) > 0:
            raise ValueError(f"Following columns are not in '{name_df}': {wrong_columns}")

    # Check if NaNs are present when they are not accepted
    if not accept_nan:
        if df[cols].isna().any().any():
            raise ValueError(f"NaN values are not allowed in '{cols}'.")