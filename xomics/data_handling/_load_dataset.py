"""
This is a script for protein benchmark loading function. Please define new loading functions by their loaded data
by introducing a new data table in docs/source/index/tables_templates.rst.
"""
import os
from pandas import DataFrame
from typing import Optional, Literal
import xomics.utils as ut

# Constants
FOLDER_USE_CASES = ut.FOLDER_DATA + "use_cases" + ut.SEP


# I Helper Functions
# Check functions
def check_name_of_dataset(name="Overview", folder_in=None):
    """"""
    if name == "Overview":
        return
    list_datasets = [x.split(".")[0] for x in os.listdir(folder_in)
                     if "." in x and not x.startswith(".") and ".tsv" in x]
    if name not in list_datasets:
        raise ValueError(f"'name' ({name}) should be one of the following. {list_datasets}")


# II Main Functions
def load_dataset(name: str = "Overview",
                 n: Optional[int] = None,
                 random: bool = False,
                 drop_na: bool = False,
                 ) -> DataFrame:
    """
    Loads protein benchmarking datasets.

    Parameters
    ----------
    name
        The name of the loaded dataset, from the 'Dataset' column in the overview table.
    n
        Number of proteins, selected by index. If None, the whole dataset will be returned.
    random
        If True, ``n`` randomly selected proteins will be chosen.
    drop_na
        If True, rows containing any missing value will be dropped.

    Returns
    -------
    df_q
        DataFrame with quantifiaction values for n samples, typically proteins or genes.
        ``Rows`` correspond samples and ``columns`` to quantifications over different conditions

    Examples
    --------
    >>> import xomics as xo
    >>> df_lfq = xo.load_dataset(name="PROT_DEMYELINATION", n=100)

    Notes
    -----
    - The name of ``df_q`` can be specified for the respective quantification type such as 'df_lfq' for
      label free quantification (lfq) commonly used in mass-sepctromatery (MS)-based use_cases

    See Also
    --------
    * Overview of all benchmarks in :ref:`t1_overview_datasets`.
    * Step-by-step guide in the `Data Loading Tutorial <tutorial2a_data_loader.html>`_.
    """
    check_name_of_dataset(name=name, folder_in=FOLDER_USE_CASES)
    ut.check_number_range(name="n", val=n, min_val=1, accept_none=True, just_int=True)
    ut.check_bool(name="drop_na", val=drop_na)
    ut.check_bool(name="random", val=random)
    # Load overview table
    if name == "Overview":
        return ut.read_excel_cached(FOLDER_USE_CASES + "Overview.xlsx")
    df = ut.read_csv_cached(FOLDER_USE_CASES + name + ".tsv", sep="\t")
    df.columns = [c.lower().replace(" ", "_") for c in list(df)]
    for col in [ut.COL_PROT_ID, ut.COL_PROT_NAME, ut.COL_GENE_NAME]:
        if col in df:
            df = df.dropna(subset=col)
    if drop_na:
        df = df.dropna()
    if n is not None:
        if random:
            df = df.sample(n)
        else:
            df = df.head(n)
    # Adjust index
    df_q = df.reset_index(drop=True)
    return df_q
