"""
This a script for general decorators used in xOmics.
# Dev: use runtime decorator only for internal methods since they destroy the signature for some IDEs
"""
import warnings
import traceback
from sklearn.exceptions import ConvergenceWarning, UndefinedMetricWarning
import functools
import re

# Helper functions


# Catch Runtime
class CatchRuntimeWarnings:
    """Context manager to catch RuntimeWarnings and store them in a list."""
    def __enter__(self):
        self._warn_list = []
        self._other_warnings = []
        self._showwarning_orig = warnings.showwarning
        warnings.showwarning = self._catch_warning
        return self

    def __exit__(self, exc_type, exc_value, tb):
        warnings.showwarning = self._showwarning_orig
        # Re-issue any other warnings that were caught but not RuntimeWarning
        for warn_message, warn_category, filename, lineno in self._other_warnings:
            warnings.warn_explicit(warn_message, warn_category, filename, lineno)

    def _catch_warning(self, message, category, filename, lineno, file=None, line=None):
        if category == RuntimeWarning:
            line_content = traceback.format_list([(filename, lineno, "", line)])[0].strip()
            warning_msg = f"{message}: {line_content.split(', in')[1]}"
            self._warn_list.append(warning_msg)
        else:
            # Store other warnings for re-issuing later
            self._other_warnings.append((message, category, filename, lineno))

    def get_warnings(self):
        return self._warn_list


def catch_runtime_warnings():
    """Decorator to catch RuntimeWarnings and store them in a list.

    Returns:
        decorated_func: The decorated function
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with CatchRuntimeWarnings() as crw:
                result = func(*args, **kwargs)
            if crw.get_warnings():
                list_warnings = crw.get_warnings()
                n = len(list_warnings)
                summary_msg = f"The following {n} 'RuntimeWarnings' were caught:\n" + "\nRuntimeWarning: ".join(crw.get_warnings())
                warnings.warn(summary_msg, RuntimeWarning)
            return result
        return wrapper

    return decorator


# Catch convergence
class ClusteringConvergenceException(Exception):
    def __init__(self, message, distinct_clusters):
        super().__init__(message)
        self.distinct_clusters = distinct_clusters


def catch_convergence_warning():
    """Decorator to catch ConvergenceWarnings and raise custom exceptions.

    Returns:
        decorated_func: The decorated function.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with warnings.catch_warnings(record=True) as w:
                # Trigger the "always" behavior for ConvergenceWarning
                warnings.simplefilter("always", ConvergenceWarning)
                result = func(*args, **kwargs)  # Call the original function

                # Check if the warning is the one we're interested in
                for warn in w:
                    if issubclass(warn.category, ConvergenceWarning):
                        message = str(warn.message)
                        if "Number of distinct clusters" in message:
                            distinct_clusters = int(message.split("(")[1].split(")")[0].split()[0])
                            raise ClusteringConvergenceException(f"Process stopped due to ConvergenceWarning.", distinct_clusters)
            return result
        return wrapper

    return decorator


# Catch invalid division (could be added to AAclust().comp_medoids())
class InvalidDivisionException(Exception):
    pass


def catch_invalid_divide_warning():
    """Decorator to catch specific RuntimeWarnings related to invalid division
       and raise custom exceptions.

    Returns:
        decorated_func: The decorated function.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with CatchRuntimeWarnings() as crw:
                result = func(*args, **kwargs)
            if crw.get_warnings():
                raise InvalidDivisionException(f"\nError due to 'RuntimeWarning': {crw.get_warnings()[0]}")
            return result
        return wrapper
    return decorator

# Catch UndefinedMetricWarnings
class CatchUndefinedMetricWarning:
    """Context manager to catch and aggregate UndefinedMetricWarnings."""
    def __enter__(self):
        self._warn_set = set()
        self._other_warnings = []
        self._showwarning_orig = warnings.showwarning
        warnings.showwarning = self._catch_warning
        return self

    def __exit__(self, exc_type, exc_value, tb):
        warnings.showwarning = self._showwarning_orig
        for warn_message, warn_category, filename, lineno in self._other_warnings:
            warnings.warn_explicit(warn_message, warn_category, filename, lineno)

    def _catch_warning(self, message, category, filename, lineno, file=None, line=None):
        if category == UndefinedMetricWarning:
            self._warn_set.add(str(message))  # Add message to set (duplicates are automatically handled)
        else:
            self._other_warnings.append((message, category, filename, lineno))

    def get_warnings(self):
        return list(self._warn_set)

def catch_undefined_metric_warning():
    """Decorator to catch and report UndefinedMetricWarnings once per unique message."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with CatchUndefinedMetricWarning() as cumw:
                result = func(*args, **kwargs)
            if cumw.get_warnings():
                summary_msg = "The following 'UndefinedMetricWarning' was caught:\n" + "\n".join(cumw.get_warnings())
                summary_msg += ("\n This warning was likely triggered due to 'precision' or 'f1' metrics and "
                                "an imbalanced and/or small dataset.")
                warnings.warn(summary_msg, UndefinedMetricWarning)
            return result
        return wrapper
    return decorator