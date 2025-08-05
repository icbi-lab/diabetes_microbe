# diabetes_microbe

![Version Badge](https://img.shields.io/badge/Version-1.0.2-brightgreen?style=for-the-badge)

## The  gut microbiota in diabetes mellitus (T3cDM vs T1DM)

Metagenomic study of the  gut microbiota between patients with pancreoprive diabetes (T3cDM) and type 1 diabetes (T1DM)

This study consists of 48 patients ... 

## Repository structure

01-tables: Contains initial metadata about samples and data files that are used for downstream analysis. This includes locations of comparative data.

02-scripts: Contains all scripts and code notebooks used in the day-to-day analysis during the project. Can optionally include sub-directories for each language (e.g., R, python).

05-results: Contains figures and final figures

Raw data is uploaded to zenodo here. 

## Dependencies

Key software and packages used include:

- Nextflow (workflow manager): nf-core/ampliseq v2.12.0
- cutadapt v4.6
- FastQC v0.12.1
- DADA2 v1.30.0
- QIIME2 v2023.7.0
- SILVA v138.1
- R v4.3.0
- python=3.9.4 

Python packages mainly used: pandas==1.5.3, seaborn==0.13.2, matplotlib-base=3.8.4, upsetplot==0.9.0. For more details on Python dependencies see diabetes_microbe.yml

For more details on R packages see R_session_info.txt


## Data preprocessing
- 16S rRNA amplicon data from DNA. The merged fastq files are stored in zenodo. 
- The fastq files are analysed using nf-core/ampliseq v2.12.0

```
sbatch 02-scripts/nf-core_ampliseq/run_nf_core_ampliseq.slurm
```

## Data anaylsis

- Compositional analysis

```
02-scripts/01_compositional.R
```

```
05-results/figures/compositional_plots_familiy.png
```

- Alpha diversity

```
02-scripts/04_alpha diversity.R
```
- Beta diversity

```
02-scripts/05_beta_diversity.ipynb
```

- Microbial logistic regression model 
```
02-scripts/06_logistic_regression_H_T1DM.R
02-scripts/06_logistic_regression_H_T3cDM.R
02-scripts/06_logistic_regression_T1DM_T3cDM.R
```

- Microbial community metabolic modeling
```
02-scripts/materials/08_micom.ipynb
```
- Correlation Analyses
```
02-scripts/012_zscore_heatmap.R
02-scripts/013_taxa_clinical_correlation.R
```

- Linear discriminant analysis for microbiome data using LEfSe

```
02-scripts/011_lefse.R
02-scripts/016_lefse_sex.R
```

- Linear Regression

```
02-scripts/014_linear_regression.R
```

- Longitudinal linear mixed models (LMM)
```
02-scripts/015_metadata_analysis.R
```

## Authors

Erika Kvalem


