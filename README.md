# GENOMIC_PREDICTION
Project Overview: Implements a genomic prediction pipeline using GBLUP, combining BLUE-adjusted phenotypes and SNP marker data to estimate Genomic Estimated Breeding Values (GEBVs) for selection in plant breeding.

Note: The dataset is hypothetical and includes intentional inconsistencies

Folder Structure: 
Gene_Prediction_Script/   : R scripts for genomic prediction (GBLUP)
Gene_Prediction_Data/     : Input phenotypic and genotypic datasets
Gene_Prediction_Results/  : Output results (BLUEs, GEBV rankings)

Steps:
1. Phenotypic data cleaning and quality assessment
2. Experimental design validation (Location, Replication)
3. BLUE estimation using mixed models
4. Genotype–phenotype alignment
5. SNP missing data imputation
6. Genomic Relationship Matrix (VanRaden method)
7. GBLUP model fitting
8. Extraction and ranking of GEBVs

Notes:
Phenotypic data are adjusted using BLUEs prior to genomic prediction.
SNP markers are coded as 0/1/2 with mean imputation for missing values.
The dataset is hypothetical and includes intentional inconsistencies to demonstrate robust data handling.
All scripts are ready to run, and outputs are saved in Gene_Prediction_Results/.
