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
All tables from the xOmics documentation are listed here, in chronological order based on the project history.

ADD-TABLE

.. _t1_overview_proteomics_XXX:

Proteomics Example Datasets
---------------------------
Example proteomic datasets are given for different mouse models for neurodegenerative diseases.

Datasets are named according to the assessed disease (e.g., Alzheimer´s disease (AD)) and the name of the mouse models,
as described in `Alzforum <https://www.alzforum.org/research-models>`_ or defined in the respective publication. Datasets
were obtained by mass spectrometry (MS)-based proteomics.

ADD-TABLE

.. _t2_omics_analysis_tools_XXX:

Omics Analysis Tools
--------------------
Overview of different omics analysis software tools such as `MaxQuant <https://www.maxquant.org/>`_ or
`DIA-NN <https://www.nature.com/articles/s41592-019-0638-x>`_ for proteomics data are given.

ADD-TABLE

.. _t3_omics_post_analysis_tools_XXX:

Post-Analysis Tools
-------------------
Post-analysis tools for omics data are diverse software solutions that facilitate specialized types of data evaluations,
like differential gene expression analysis. These tools span from Graphical User Interface (GUI) applications
such as `Perseus <https://maxquant.net/perseus/>`_ to Python-based packages tailored for specific analyses, such as
`Scanpy <https://scanpy.readthedocs.io/en/stable/>`_ for single-cell RNAseq data analysis.

ADD-TABLE

.. _t4_gene_enrichment_tools_XXX:

Gene Enrichment Tools
---------------------
Gene enrichment analysis for omics data is a computational method used to identify which predefined sets of genes
or proteins are statistically over-represented in a large set of genes or proteins. It helps in deciphering the
biological significance behind large-scale molecular data by linking genes to known pathways, functions, or other
biological categories.

ADD-TABLE