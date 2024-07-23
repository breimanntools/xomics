Welcome to the xOmics documentation
===================================
.. Developer Notes:
    Please update badges in README.rst and vice versa
.. Group 1: Package badges
.. |PyPI Status| image:: https://img.shields.io/pypi/status/xomics.svg
   :target: https://pypi.org/project/xomics/
   :alt: PyPI - Status

.. |PyPI Version| image:: https://img.shields.io/pypi/v/xomics.svg
   :target: https://pypi.python.org/pypi/xomics
   :alt: PyPI - Package Version

.. |Supported Python Versions| image:: https://img.shields.io/pypi/pyversions/xomics.svg
   :target: https://pypi.python.org/pypi/xomics
   :alt: Supported Python Versions

.. |Downloads| image:: https://pepy.tech/badge/xomics
   :target: https://pepy.tech/project/xomics
   :alt: Downloads

.. |License| image:: https://img.shields.io/github/license/breimanntools/xomics.svg
   :target: https://github.com/breimanntools/xomics/blob/master/LICENSE
   :alt: License

.. Group 2: Testing badges
.. |Unit Tests| image:: https://github.com/breimanntools/xomics/actions/workflows/main.yml/badge.svg
   :target: https://github.com/breimanntools/xomics/actions/workflows/main.yml
   :alt: CI/CD Pipeline

.. |CodeQL| image:: https://github.com/breimanntools/xomics/actions/workflows/codeql_analysis.yml/badge.svg
   :target: https://github.com/breimanntools/xomics/actions/workflows/codeql_analysis.yml
   :alt: CodeQL

.. |Codecov| image:: https://codecov.io/gh/breimanntools/xomics/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/breimanntools/xomics
   :alt: Codecov

.. |Documentation Status| image:: https://readthedocs.org/projects/xomics/badge/?version=latest
   :target: https://xomics.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


.. Group 3: Potential badges for future
.. |Conda Version| image:: https://anaconda.org/conda-forge/xomics/badges/version.svg
   :target: https://anaconda.org/conda-forge/xomics
   :alt: Conda - Package Version


..
    Missing badges
    |Conda Version|

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - **Package**
     - |PyPI Status| |PyPI Version| |Supported Python Versions| |Downloads| |License|
   * - **Testing**
     - |Unit Tests| |CodeQL| |Codecov| |Documentation Status|

**xOmics** (eXplainable Omics) is a Python framework developed for interpretable omics analysis,
focusing on differential proteomics expression data. It introduces the following key algorithms:

- **cImpute**: Conditional Imputation - A transparent method for hybrid missing value imputation.
- **pRank**: Protein-centric ranking of (prote)omic data, including integration with functional enrichment results.
- **pIntegrate**: Protein-centric integration of multiple (prote)omic datasets for differential analysis.
- **QARIP**: Quantitative proteomic analysis of regulated intramembrane proteolysis.

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

**QARIP**:
    *QARIP: a web server for quantitative proteomic analysis of regulated intramembrane proteolysis*
    `Nucleic Acids Research <https://academic.oup.com/nar/article/41/W1/W459/1105195>`__.
