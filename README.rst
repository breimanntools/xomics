Welcome to the xOmics documentation
===================================
.. Developer Notes:
    Please update badges in README.rst and vice versa
.. image:: https://github.com/breimanntools/xomics/workflows/Build/badge.svg
   :target: https://github.com/breimanntools/xomics/actions
   :alt: Build Status

.. image:: https://github.com/breimanntools/xomics/workflows/Python-check/badge.svg
   :target: https://github.com/breimanntools/xomics/actions
   :alt: Python-check

.. image:: https://img.shields.io/pypi/status/xomics.svg
   :target: https://pypi.org/project/xomics/
   :alt: PyPI - Status

.. image:: https://img.shields.io/pypi/pyversions/xomics.svg
   :target: https://pypi.python.org/pypi/xomics
   :alt: Supported Python Versions

.. image:: https://img.shields.io/pypi/v/xomics.svg
   :target: https://pypi.python.org/pypi/xomics
   :alt: PyPI - Package Version

.. image:: https://anaconda.org/conda-forge/xomics/badges/version.svg
   :target: https://anaconda.org/conda-forge/xomics
   :alt: Conda - Package Version

.. image:: https://readthedocs.org/projects/xomics/badge/?version=latest
   :target: https://xomics.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://img.shields.io/github/license/breimanntools/xomics.svg
   :target: https://github.com/breimanntools/xomics/blob/master/LICENSE
   :alt: License

.. image:: https://pepy.tech/badge/xomics
   :target: https://pepy.tech/project/xomics
   :alt: Downloads

**xOmics** (eXplainable Omics) is a Python framework developed for interpretable omics analysis,
focusing on differential proteomics expression data. It introduces the following key algorithms:

- **cImpute**: Conditional Imputation - A transparent method for hybrid missing value imputation.
- **pRank**: Protein-centric ranking of (prote)omic data, including integration with functional enrichment results.
- **pIntegrate**: Protein-centric integration of multiple (prote)omic datasets for differential analysis.

In addition, **xOmics** provides functional capabilities for efficiently loading benchmark proteomics datasets via
**load_datasets**, accompanied by corresponding enrichment data.A suite of supportive functions is also available to
facilitate a smooth and efficient (prote)omic analysis pipeline.

Install
=======
**xOmics** can be installed either from `PyPi <https://pypi.org/project/xomics>`_ or
`conda-forge <https://anaconda.org/conda-forge/xomics>`_:

.. code-block:: bash

   pip install -u xomics
   or
   conda install -c conda-forge xomics

Contributing
============
We appreciate bug reports, feature requests, or updates on documentation and code. For details, please refer to
`Contributing Guidelines <CONTRIBUTING.rst>`_. For further questions or suggestions, please email stephanbreimann@gmail.com.

Citations
=========
If you use xOmics in your work, please cite the respective publication as follows:

**xOmics**:
   [Citation details and link if available]

**cImpute**:
   [Citation details and link if available]

**pRank**:
   [Citation details and link if available]

**pIntegrate**:
   [Citation details and link if available]
