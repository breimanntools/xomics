..
   Developer Notes:
   This is the landing page for the xOmics documentation using Sphinx, containing the root `toctree` directive.
   The documentation will be hosted on Read the docs.
..

Welcome to the xOmics documentation!
====================================
.. include:: index/badges.rst
.. include:: index/overview.rst

Install
=======
**xOmics** can be installed either from `PyPi <https://pypi.org/project/xomics>`_ or
`conda-forge <https://anaconda.org/conda-forge/xomics>`_:

.. code-block:: bash

   pip install -u xomics
   or
   conda install -c conda-forge xomics


.. toctree::
   :maxdepth: 1
   :caption: OVERVIEW

   index/introduction.rst
   index/usage_principles.rst
   index/CONTRIBUTING_COPY.rst

.. toctree::
   :maxdepth: 1
   :caption: EXAMPLES

   tutorials.rst

.. toctree::
   :maxdepth: 2
   :caption: REFERENCES

   api.rst

.. toctree::
   :maxdepth: 1

   index/tables.rst
   index/references.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citation
========

.. include:: index/citations.rst
