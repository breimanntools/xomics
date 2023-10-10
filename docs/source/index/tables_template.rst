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
as described in <https://www.alzforum.org/research-models> or defined in the respective publication. Datasets
were obtained by mass spectrometry (MS)-based proteomics.

ADD-TABLE
