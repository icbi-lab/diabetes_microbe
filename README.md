# diabetes_microbe

![Version Badge](https://img.shields.io/badge/Version-1.0.2-brightgreen?style=for-the-badge)

## The  gut microbiota in diabetes mellitus (T3cDM vs T1DM)

Metagenomic study of the  gut microbiota between patients with pancreoprive diabetes (T3cDM) and type 1 diabetes (T1DM)

This study consists of 48 patients ... 

## Repository structure

01-tables: Contains initial input files for preprocessing, main downstream analysis and also supplementary analysis. Contains also statistical results.

02-scripts: Contains all the scripts for downstream analysis and producing the figures. 

05-results: Contains final figures

Raw data is uploaded to zenodo here. https://doi.org/10.5281/zenodo.16794434

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

Figure 1

- Alpha diversity

```
02-scripts/Figure1a_alpha_diversity.ipynb
02-scripts/Figure1a_simpson.R
```
- Beta diversity

```
02-scripts/Figure1b_beta_diversity.ipynb
```
- Upset plot 
```
02-scripts/Figure1c_upsetplot.ipynb
```
- PCA
```
02-scripts/Figure1d_pca.ipynb
``` 

Figure 2

- Compositional analysis

```
02-scripts/Figure2a_compositional.R
```
- Core microbiome

```
02-scripts/Figure2b_core_microbiome_heatmap.R
02-scripts/Figure2c_violin_relative_abundance.R
```

Figure 3

- Microbial logistic regression model 

```
02-scripts/Figure3abc_logistic_regression_T1DM_T3cDM_CV.R
```

- Venn diagram 

```
02-scripts/Figure3d_venn_diagram.R
```

- Genera of interest 

```
02-scripts/Figure3e_violin_relative_abundance_genus.R
```


Figure 4

- Microbial community metabolic modeling
```
02-scripts/materials/micom.ipynb
02-scripts/materials/Figure4a_growthrate_boxplot.R
02-scripts/materials/Figure4b_metabolome_barplot_all_taxa.R
02-scripts/materials/Figure4c_metabolome_barplot_eshigella.R
```



## Authors

Erika Kvalem


