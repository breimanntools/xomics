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

.. _tables:

Tables
======

.. contents::
    :local:
    :depth: 1

.. _t0_mapper:

Overview Table
--------------
All tables from the xOmics documentation are listed here, in chronological order based on the project history.


.. list-table::
   :header-rows: 1
   :widths: 8 8 8

   * - Table
     - Description
     - See Also
   * - t1_overview_proteomics
     - Proteomic example datasets
     - aa.load_dataset
   * - t2_omics_analysis_tools
     - Analysis tools for omics data
     - nan
   * - t3_omics_post-analysis_tools
     - Post-analysis tools for omics data
     - nan
   * - t4_gene_enrichment_tools
     - Gene enrichment analysis tools
     - nan


.. _t1_overview_proteomics:

Proteomics Example Datasets
---------------------------
Example proteomic datasets are given for different mouse models for neurodegenerative diseases.

Datasets are named according to the assessed disease (e.g., Alzheimer´s disease (AD)) and the name of the mouse models,
as described in `Alzforum <https://www.alzforum.org/research-models>`_ or defined in the respective publication. Datasets
were obtained by mass spectrometry (MS)-based proteomics.


.. list-table::
   :header-rows: 1
   :widths: 8 8 8

   * - Dataset
     - Description
     - Reference
   * - XX
     - XX
     - :ref:`Monasor20 <Monasor20>`


.. _t2_omics_analysis_tools:

Omics Analysis Tools
--------------------
Overview of different omics analysis software tools such as `MaxQuant <https://www.maxquant.org/>`_ or
`DIA-NN <https://www.nature.com/articles/s41592-019-0638-x>`_ for proteomics data are given.


.. list-table::
   :header-rows: 1
   :widths: 8 8 8 8 8 8 8 8 8 8 8

   * - Tool
     - Description
     - Usability
     - Open-source
     - Documentation
     - Community & Support
     - Integration
     - Programming Language/GUI
     - Advantages
     - Disadvantages
     - Publication
   * - MaxQuant
     - Proteomics data analysis, especially for label-free quantification
     - Proteomics
     - No
     - Excellent
     - Active
     - Moderate
     - Standalone / GUI
     - Robust algorithms, widely used
     - Requires high computational resources
     - [Link to paper]
   * - Spectronaut
     - Analysis of DIA (data-independent acquisition) mass spectrometry data
     - Proteomics
     - No
     - Good
     - Managed by Biognosys
     - Limited
     - Standalone / GUI
     - Optimized for DIA, high reproducibility
     - Proprietary software
     - [Link to paper]
   * - DIA-NN
     - Software suite for DIA data analysis
     - Proteomics
     - Yes
     - Good
     - Growing
     - Good
     - Command-line
     - Open-source, versatile
     - Command-line based
     - [Link to paper]
   * - Skyline
     - Targeted mass spec data analysis
     - Proteomics
     - Yes
     - Excellent
     - Active
     - Excellent
     - Standalone / GUI
     - Supports multiple instrument vendors, extensible
     - Mainly for targeted proteomics
     - [Link to paper]
   * - LipidSearch
     - Software for lipidomics data processing and identification
     - Lipidomics
     - No
     - Good
     - Managed by Thermo Fisher
     - Moderate
     - Standalone / GUI
     - Comprehensive lipid databases, integration with mass spec instruments
     - Proprietary software
     - [Link to paper]
   * - LipidHunter
     - Identification of lipids from LC-MS/MS data
     - Lipidomics
     - Yes
     - Good
     - Active
     - Good
     - Python
     - Open-source, comprehensive output
     - Requires good understanding of lipidomics
     - [Link to paper]
   * - MZmine
     - Framework for processing, visualization, and analysis of mass spectrometry data
     - Metabolomics
     - Yes
     - Good
     - Active
     - Good
     - Java
     - Modular, supports various data processing tasks
     - Java-centric, learning curve
     - [Link to paper]
   * - MetaboAnalyst
     - Comprehensive web-based tool for metabolomics data analysis
     - Metabolomics
     - Yes
     - Excellent
     - Active
     - Good
     - Web-based
     - Wide range of statistical methods, user-friendly interface
     - Web-based, can limit very large analyses
     - [Link to paper]
   * - XCMS
     - Processing and analysis of untargeted metabolomics data
     - Metabolomics
     - Yes
     - Good
     - Active
     - Good
     - R
     - Widely used in the community, high flexibility
     - Requires R programming knowledge
     - [Link to paper]
   * - Compound Discoverer
     - Software for metabolite identification and quantitative analysis
     - Metabolomics
     - No
     - Good
     - Managed by Thermo Fisher
     - Moderate
     - Standalone / GUI
     - Comprehensive workflow, integration with mass spec instruments
     - Proprietary software
     - [Link to paper]


.. _t3_omics_post-analysis_tools_XXX:

Post-Analysis Tools
-------------------
Post-analysis tools for omics data are diverse software solutions that facilitate specialized types of data evaluations,
like differential gene expression analysis. These tools span from Graphical User Interface (GUI) applications
such as `Perseus <https://maxquant.net/perseus/>`_ to Python-based packages tailored for specific analyses, such as
`Scanpy <https://scanpy.readthedocs.io/en/stable/>`_ for single-cell RNAseq data analysis.


.. list-table::
   :header-rows: 1
   :widths: 8 8 8 8 8 8 8 8 8 8 8

   * - Tool
     - Description
     - Usability
     - Open-source
     - Documentation
     - Community & Support
     - Integration
     - Programming Language/GUI
     - Advantages
     - Disadvantages
     - Publication
   * - MaxQuant
     - Proteomics data analysis, especially for label-free quantification
     - Proteomics
     - No
     - Excellent
     - Active
     - Moderate
     - Standalone / GUI
     - Robust algorithms, widely used
     - Requires high computational resources
     - [Link to paper]
   * - Spectronaut
     - Analysis of DIA (data-independent acquisition) mass spectrometry data
     - Proteomics
     - No
     - Good
     - Managed by Biognosys
     - Limited
     - Standalone / GUI
     - Optimized for DIA, high reproducibility
     - Proprietary software
     - [Link to paper]
   * - DIA-NN
     - Software suite for DIA data analysis
     - Proteomics
     - Yes
     - Good
     - Growing
     - Good
     - Command-line
     - Open-source, versatile
     - Command-line based
     - [Link to paper]
   * - Skyline
     - Targeted mass spec data analysis
     - Proteomics
     - Yes
     - Excellent
     - Active
     - Excellent
     - Standalone / GUI
     - Supports multiple instrument vendors, extensible
     - Mainly for targeted proteomics
     - [Link to paper]
   * - LipidSearch
     - Software for lipidomics data processing and identification
     - Lipidomics
     - No
     - Good
     - Managed by Thermo Fisher
     - Moderate
     - Standalone / GUI
     - Comprehensive lipid databases, integration with mass spec instruments
     - Proprietary software
     - [Link to paper]
   * - LipidHunter
     - Identification of lipids from LC-MS/MS data
     - Lipidomics
     - Yes
     - Good
     - Active
     - Good
     - Python
     - Open-source, comprehensive output
     - Requires good understanding of lipidomics
     - [Link to paper]
   * - MZmine
     - Framework for processing, visualization, and analysis of mass spectrometry data
     - Metabolomics
     - Yes
     - Good
     - Active
     - Good
     - Java
     - Modular, supports various data processing tasks
     - Java-centric, learning curve
     - [Link to paper]
   * - MetaboAnalyst
     - Comprehensive web-based tool for metabolomics data analysis
     - Metabolomics
     - Yes
     - Excellent
     - Active
     - Good
     - Web-based
     - Wide range of statistical methods, user-friendly interface
     - Web-based, can limit very large analyses
     - [Link to paper]
   * - XCMS
     - Processing and analysis of untargeted metabolomics data
     - Metabolomics
     - Yes
     - Good
     - Active
     - Good
     - R
     - Widely used in the community, high flexibility
     - Requires R programming knowledge
     - [Link to paper]
   * - Compound Discoverer
     - Software for metabolite identification and quantitative analysis
     - Metabolomics
     - No
     - Good
     - Managed by Thermo Fisher
     - Moderate
     - Standalone / GUI
     - Comprehensive workflow, integration with mass spec instruments
     - Proprietary software
     - [Link to paper]


.. _t4_gene_enrichment_toolsXXX:

Gene Enrichment Tools
---------------------
Gene enrichment analysis for omics data is a computational method used to identify which predefined sets of genes
or proteins are statistically over-represented in a large set of genes or proteins. It helps in deciphering the
biological significance behind large-scale molecular data by linking genes to known pathways, functions, or other
biological categories.

ADD-TABLE