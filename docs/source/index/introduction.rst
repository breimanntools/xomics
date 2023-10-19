Introduction
============

**xOmics** (eXplainable Omics) is a Python framework developed for streamlined and explainable omics analysis, with a
spotlight on differential proteomics expression data.

Workflow
--------
A typical workflow involves the following algorithms:

- **cImpute**: Conditional Imputation - A transparent method for hybrid missing value imputation.
- **pIntegrate**: Protein-centric integration of multiple (prote)omic datasets to find commonalities and differences.
- **pRank**: Protein-centric ranking of (prote)omic data, leveraging functional enrichment results.

Missing value imputation, data integration, and protein-centric ranking can all be performed independently. See examples
and practical usage in our :ref:`tutorials <tutorials>`.

Data Flow: Components of xOmics
-------------------------------

The xOmics toolkit uses different DataFrames. It starts with DataFrames containing quantifications for omics datasets
and (optionally) corresponding enrichment data (**df_q**, **df_enrich**). The omics quantification data can be first
imputed and then used integration and protein-centric ranking, where enrichment data can be included as well.
Various plotting functions are provided to enable data interpretation.

See the primary data flow within the xOmics toolkit in this diagram:

.. image:: /_artwork/diagrams/components.png

Entry Points: Bridges to External Software
------------------------------------------
