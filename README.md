cImpute (conditional Imputation) is a hybrid imputation algorithm for missing values (MVs) in (prote)omics data.
Missing values can be distinguished into three categories as described by Lazar et al., 2016 and Wei et al., 2018
for proteomic data sets as follows:

    a) Missing Completely At Random (MCAR): MVs due to random errors and stochastic fluctuations during process of 
        data acquisition. Since MCAR MVs can not be explained by measured intensities, they are uniformly distributed.
    b) Missing At Random (MAR): MVs due to suboptimal data processing and conditional dependencies. MAR is a more 
        general class than MCAR, where all MCAR MVs are MAR MVs. The distribution of MAR MVs can just be speculated 
        and likely differs highly between experiments. 
    c) Missing Not At Random (MNAR): MVs due to experimental bias (i.e., the detection limit in mass-spectrometry
        experiments). Commonly, MNAR MVs are described by a left-censored Gaussian distribution (i.e., the Gaussian
        distribution is truncated on the region of lower abundances, which is the left side of the distribution).

cImpute aims to impute only MVs matching well-defined confidence criteria and consists of following four steps:

1. Definition of the upper bound for MNAR MVs called upMNAR to distinguish between MNAR and MCAR.
    Let the detection range (DR) be defined as

        DR = Dmax - Dmin

    where Dmax is the smallest detected value (i.e., detection limit) and Dmax is the largest detected value over
    the a whole data set. Such that the upper bound for MNAR MVs is defined as

        upMNAR = Dmin + DR * l

    where l is the location factor for the upMNAR given as relative proportion of the DR (0-1, by default 0.1).

2. For each detected protein in an experimental group, cImpute distinguishes between following four MVs scenarios:
    
    a) MCAR if all values > upMNAR
    b) MNAR if all values <= upMNAR
    c) NM if no MVs occur
    d) MAR otherwise.

3. A confidence score (CS) [0-1] will be assigned to each detected protein in an experimental group based on
    the fraction of MVs and the MVs category:

    a) For MCAR, CS is equal to the fraction of quantified values (non MVs)
    b) For MNAR, CS is equal to the fraction of MVs
    d) For NM, CS is set to 1
    c) For MAR, CS is set to 0 (due to the lack of certainty about the source of the MVs)

4. Imputation methods are applied on MCAR and MNAR MVs (MinProb and KNN, respectively)
    with a CS higher than a given threshold for each experimental group (0.5 by default).

The upMNAR (step 1) is computed over the whole dataset and steps 2-4 are preformed for each experimental group
    (c.f., Liu and Dongre, 2020).

As a result, just MVs are imputed fulfilling a certain confidence level and a high transparency is guaranteed
due to the missing value classification and the given confidence score.

References
----------
Lazar et al., 2016, Accounting for the Multiple Nature of Missing Values in Label-Free Quantitative Proteomics
    Data Sets to Compare Imputation Strategies (Journal of Proteomics Research)
Liang et al., 2021, A comparative study of evaluating missing value imputation methods in label-free proteomics
    (scientific reports)
Palstrom et al., 2020, Data Imputation in Merged Isobaric Labeling-Based Relative Quantification Datasets
    (Mass Spectrometry Data Analysis in Proteomics, Springer Protocols)
Liu and Dongre, 2020, Proper imputation of missing values in proteomics datasets for differential expression analysis
    (Briefings in Bioinformatics)

See also
--------
    https://www.rdocumentation.org/packages/imputeLCMD/versions/2.0/topics/impute.MinProb
    https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/MissingValues.html
