Data Flow and Enry Points
=========================

Data Flow: Components of xOmics
-------------------------------

The xOmics toolkit uses different DataFrames. It starts with DataFrames containing omics datasets and (optionally)
corresponding functional enrichment data (**df_omics**, **df_enrich**). Data imputation and ranking can be directly
applied to the omics data, where for the latter enrichment data can be utilized as well. Multiple omics datasets
can be integrated and subsequently ranked. Various plotting functions are provided to enable data interpreation.

See the primary data flow within the xOmics toolkit in this diagram:

.. image:: /_artwork/diagrams/components.png

Entry Points: Bridges to External Libraries
-------------------------------------------
