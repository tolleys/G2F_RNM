# G2F_RNM
Within this repository, we are sharing the R script used in the analysis for Investigating the Genomic Background and Predictive Ability of Genotype-by-environment Interactions in Maize Grain Yield Based on Reaction Norm Models.  

Phenotype, genotype, weather, and metadata are provided at Genomes 2 Fields https://www.genomes2fields.org/resources/

Script Descriptions

1.G2F_PHENO_GENO_Prep_sat.R was used to prepare the data from the Genomes 2 Fields experiment for BLUPF90.

2.G2F_NASAPowerEnvironmentalGradientPrediction_sat.R was used to use weather and soil data to estimate the merit of a given environment in a machine learning model.

3.G2F_GeneticParameters_sat.R used the output from BLUPF90 to evaluate heritability, genetic correlation, SNP effect, and linkage disequilibrium.

4.G2F_GenomicPrediciton_sat.R used the output from BLUPF90 to evaluate genomic prediction of various validation scenarios. 
