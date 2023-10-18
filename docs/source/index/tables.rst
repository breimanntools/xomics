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
   * - t3_omics_post_analysis_tools
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


.. _t3_omics_post_analysis_tools:

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
   * - Perseus
     - Comprehensive platform for in-depth analysis of proteomics data
     - Proteomics
     - No
     - Good
     - Moderate
     - Limited
     - GUI
     - Comprehensive analysis for MaxQuant data
     - Limited to specific datasets
     - [Link to paper]
   * - PEPPI
     - Tool for analyzing protein-protein interactions and functional associations
     - Proteomics
     - Unknown
     - Moderate
     - Unknown
     - Unknown
     - Likely GUI
     - Protein interaction analysis
     - Unknown support and documentation
     - Unknown
   * - MSstats
     - Statistical relative quantification in mass spectrometry-based proteomics
     - Proteomics
     - Yes
     - Good
     - Active
     - Moderate
     - R / GUI
     - Robust statistical framework
     - R learning curve for some
     - [Link to paper]
   * - Pyteomics
     - Collection of tools for various tasks in proteomics data analysis
     - Proteomics
     - Yes
     - Good
     - Moderate
     - Good
     - Python
     - Python-based, flexible
     - Requires Python expertise
     - [Link to paper]
   * - AlphaPept
     - Peptide identification and quantification
     - Proteomics
     - Yes
     - Good
     - Growing
     - Limited
     - Python / GUI
     - Fast and accurate peptide identification
     - Still maturing
     - [Link to paper]
   * - Seurat
     - Toolkit for quality control, analysis, and exploration of single-cell RNA-seq data
     - scRNA-seq
     - Yes
     - Excellent
     - Very Active
     - Good
     - R / GUI
     - Comprehensive scRNA-seq toolkit
     - R learning curve for some
     - [Link to paper]
   * - Scanpy
     - Analyzing and visualizing single-cell RNA-seq data with emphasis on scalability and speed
     - scRNA-seq
     - Yes
     - Excellent
     - Very Active
     - Excellent
     - Python
     - Scalable, integration with other tools
     - Python-centric
     - [Link to paper]
   * - SCope
     - Fast, scalable, and user-friendly tool for visualizing and interpreting large datasets from scRNA-seq
     - scRNA-seq
     - Yes
     - Good
     - Active
     - Good
     - Web-based
     - User-friendly, web-based
     - Limited to visualization
     - [Link to paper]
   * - AnnData
     - Handling matrix data with annotations
     - scRNA-seq
     - Yes
     - Good
     - Associated with Scanpy
     - Good
     - Python
     - Efficient data structure for large datasets
     - Primarily a data structure, not a full toolkit
     - [Link to paper]
   * - MetaboAnalyst
     - Comprehensive platform for metabolomics data analysis and interpretation
     - Metabolomics
     - Yes
     - Excellent
     - Active
     - Good
     - Web-based
     - Comprehensive, user-friendly
     - Web-based might limit large-scale analyses
     - [Link to paper]
   * - XCMS
     - LC/MS and GC/MS data preprocessing
     - Metabolomics
     - Yes
     - Excellent
     - Very Active
     - Excellent
     - R / GUI
     - Industry standard for LC/MS data
     - R learning curve for some
     - [Link to paper]
   * - MZmine
     - MS-based molecular profile data processing and analysis
     - Metabolomics, Lipidomics
     - Yes
     - Good
     - Active
     - Good
     - Java / GUI
     - Versatile and supports various data formats
     - Java-based, might be slower on large data
     - [Link to paper]
   * - LipidSearch
     - Accurate identification and quantification of lipids from LC-MS/MS data
     - Lipidomics
     - No
     - Good
     - Managed by Thermo
     - Limited
     - GUI
     - Accurate lipid identification
     - Proprietary and expensive
     - [Link to paper]
   * - LipidHunter
     - Direct annotation of lipid species from LC-MS datasets
     - Lipidomics
     - Yes
     - Moderate
     - Moderate
     - Moderate
     - Python
     - Direct lipid species annotation
     - Requires command-line experience
     - [Link to paper]


.. _t4_gene_enrichment_tools:

Gene Enrichment Tools
---------------------
Gene enrichment analysis for omics data is a computational method used to identify which predefined sets of genes
or proteins are statistically over-represented in a large set of genes or proteins. It helps in deciphering the
biological significance behind large-scale molecular data by linking genes to known pathways, functions, or other
biological categories.


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
   * - GSEA
     - Tool for gene set enrichment analysis
     - Genomics
     - Yes
     - Excellent
     - Active
     - Good
     - Java / GUI
     - Benchmark for GSEA, widely used
     - Java-centric, may be slower on huge datasets
     - [Link to paper]
   * - Enrichr
     - Web-based tool for gene set enrichment analysis
     - Genomics
     - Yes
     - Excellent
     - Active
     - Excellent
     - Web-based
     - Comprehensive databases, user-friendly interface
     - Web-based might limit very large analyses
     - [Link to paper]
   * - DAVID
     - Bioinformatics resources for gene functional classification
     - Genomics
     - No
     - Good
     - Moderate
     - Limited
     - Web-based
     - Multiple annotation tools, widely recognized
     - Outdated interface, limited updates
     - [Link to paper]
   * - WebGestalt
     - Web-based gene set analysis toolkit
     - Genomics
     - Unknown
     - Good
     - Active
     - Good
     - Web-based
     - Multiple enrichment methods, integrated databases
     - Limited by web-interface constraints
     - [Link to paper]
   * - g:Profiler
     - Functional profiling of gene lists from large-scale experiments
     - Genomics
     - Yes
     - Good
     - Active
     - Good
     - Web-based
     - Multi-level annotation, user-friendly interface
     - Web-based, can have slow response times
     - [Link to paper]
   * - PANTHER
     - Protein ANalysis THrough Evolutionary Relationships
     - Genomics
     - No
     - Excellent
     - Managed by PANTHER
     - Limited
     - Web-based
     - Classification system, evolutionary data
     - Mainly for protein-centric analysis
     - [Link to paper]
   * - Metascape
     - Tool for gene annotation and analysis resource
     - Genomics
     - Unknown
     - Good
     - Active
     - Good
     - Web-based
     - Multiple methods and databases combined
     - Limited to predefined gene sets
     - [Link to paper]
   * - LION/web
     - Lipidome isotope labeling-based ontology
     - Lipidomics
     - Unknown
     - Good
     - Growing
     - Moderate
     - Web-based
     - Comprehensive lipid databases
     - Web-based constraints
     - [Link to paper]
   * - ClueGO
     - Cytoscape plug-in to decipher functionally grouped gene ontology networks
     - Genomics
     - Unknown
     - Good
     - Active
     - Excellent
     - Cytoscape plug-in
     - Visual representation, integrates multiple data
     - Requires Cytoscape
     - [Link to paper]
   * - FAST
     - Functional Annotation of the Mammalian Genome
     - Genomics
     - Unknown
     - Good
     - Managed by FANTOM
     - Limited
     - Web-based
     - Broad mammalian genome annotation
     - Focused on mammalian genomes
     - [Link to paper]

