.. currentmodule:: xomics

.. _api:

===
API
===

This Application Programming Interface (API) is the public interface for the objects and functions of our xOmics
Python toolkit, which can be imported by:

.. code-block:: python

    import xomics as xo

You can then access all methods and objects via the `xo` alias, such as `xo.load_dataset`.

.. _data_api:

Data Handling
-------------
.. autosummary::
    :toctree: generated/


.. _imputation_api:


Imputation
----------
.. autosummary::
    :toctree: generated/

    xomics.cImpute

.. _ranking_api:


Ranking
-------
.. autosummary::
    :toctree: generated/



Integration
------------
.. autosummary::
    :toctree: generated/


.. _plot_api:


Plotting
--------
.. autosummary::
    :toctree: generated/

        xomics.plot_volcano

Plot Utilities
--------------
.. autosummary::
    :toctree: generated/

        xomics.plot_settings
        xomics.plot_legend
        xomics.plot_get_clist
        xomics.plot_gcfs

