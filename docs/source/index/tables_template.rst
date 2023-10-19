..
   Developer Notes:
   This is the index file for all tables of the xOmics documentation.
   Tables should be saved in the /tables directory. This file serves as a template
   for tables.rst, which is automatically generated based on the information here and
   in the .csv tables from the /tables directory.

   Instructions for Adding a New Table:
   1. Store the table as a .csv file in the index/tables directory. Name it using the format tX,
      where X is incremented based on the last entry's number.
   2. Update the t0_mapper.xlsx with a corresponding entry for the new table.
   3. Create a new descriptive section here that elucidates the table's columns and any
      essential data types, such as categories.

   Note: Each table should include a 'Reference' column (include exceptions in create_tables_doc.py).

   # Key Annotations for Automated Table Generation via create_tables_doc.py:
   _XXX: A string to be stripped from the references. This prevents redundancies that may result
         in broken links.
   ADD-TABLE: Placeholder indicating where tables for the corresponding section should be inserted.
..

.. _tables_XXX:

Tables
======

.. contents::
    :local:
    :depth: 1

.. _t0_mapper_XXX:

Overview Table
--------------
All tables from the xOmics documentation are listed here:

ADD-TABLE

.. _t1_omics_fields_XXX:

Omics fields
------------
All omics field of which data can be analyzed by the xOmics toolkit are summarized in the following table:

ADD-TABLE

.. _t2_quantification_methods_XXX:

Quantification methods
----------------------
The different quantification methods used in these omic fields are described in this overview:

ADD-TABLE

.. _t3_overview_datasets_XXX:

Overview of Datasets
--------------------
Example proteomic datasets are given for different mouse models for neurodegenerative diseases.

Datasets are named according to the assessed disease (e.g., Alzheimer´s disease (AD)) and the name of the mouse models,
as described in `Alzforum <https://www.alzforum.org/research-models>`_ or defined in the respective publication. Datasets
were obtained by mass spectrometry (MS)-based proteomics.

ADD-TABLE

.. _t4_omics_analysis_tools_XXX:

Omics Analysis Tools
--------------------
Overview of different omics analysis software tools such as `MaxQuant <https://www.maxquant.org/>`_ or
`DIA-NN <https://www.nature.com/articles/s41592-019-0638-x>`_ for proteomics data are given.

ADD-TABLE

.. _t5_omics_post_analysis_tools_XXX:

Post-Analysis Tools
-------------------
Post-analysis tools for omics data are diverse software solutions that facilitate specialized types of data evaluations,
like differential gene expression analysis. These tools span from Graphical User Interface (GUI) applications
such as `Perseus <https://maxquant.net/perseus/>`_ to Python-based packages tailored for specific analyses, such as
`Scanpy <https://scanpy.readthedocs.io/en/stable/>`_ for single-cell RNAseq data analysis.

ADD-TABLE

.. _t6_enrichment_tools_XXX:

Enrichment Tools
----------------
Enrichment analysis for omics data (most often genes) is a computational method used to identify which predefined sets
of genes are statistically over-represented in a large set of genes. It helps in deciphering the biological significance
behind large-scale molecular data by linking genes to known pathways, functions, or other biological categories.
While proteins are analyzed based on their gene names using `Gene Ontology terms <https://geneontology.org/>`_ or
pathway terms of databases such as `Reactome <https://reactome.org/>`_, enrichment analysis tools for lipids are
improving with the annotation scope of the  `Lipid Ontology <https://lipidomicssociety.org/interest_groups/lipid-ontology/>`_.
See on overview of diverse enrichment tools here:

ADD-TABLE