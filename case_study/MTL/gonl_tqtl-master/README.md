# gonl_tqtl
Telomere length inheritance in the Dutch population

This repository contains scripts for the following types of analysis. 

## Linear regression to assess correlation of telomere length with age and sex 

Script name: linear_regression_age_sex.R
Input data: total.df.RData # phenotype data of all individuals 

The script was used to generate Figure 1. 

## Pairwise correlations between factors: parental ages, ages at conception, telomere lengths and sex 

Script name: pairwise_cors_supplementary_fig1.R
Input data: child.df.RData # phenotype data for children and their parents 248 families

The script was used to generate Supplementary Figure 1. 


## Multiple linear regression to assess correlation of telomere length also with parental telomeres and ages at conception

Script name: MLR_telomere_inheritance.R
Input data: child.df.RData # phenotype data for children and their parents 248 families

The script was used to generate Table 1 and Supplementary table 1. 

## Multiple linear regression also including single nucleotide variations

Script name: MLR_SNV.R
Input data: child.df.RData # phenotype data for children and their parents 248 families
Input data: SNV map and ped files # available upon request


The script was used to perform regression analysis for each SNP and generate the MLR inheritance model datasheet in Supplementary table 2. 
