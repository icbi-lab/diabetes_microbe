# Preprocessing & analysis

##  Preprocessing: nf-core/ampliseq
The preprocessing was performed using nf-core/ampliseq: Amplicon sequencing analysis workflow using DADA2 and QIIME2. 

Check ./nf-core_ampliseq/run_nf_core_ampliseq.slurm for more information. 

Samplesheet available at diabetes_microbe/01-tables/nf-core_ampliseq

##  Analysis

01_compositional.R  - The stacked bar plot depicts the relative abundance for the 10 most abundant families in all samples

02_upsetplot.ipynb - UpSetPlot visualizing the set overlap for all the recorded families 

03_pca.ipynb - Principal component analysis (PCA) for all the samples where the two principal components account for most of the variance 

04_alpha_diversity.ipynb - Community α-diversity analysis described by Shannon, Simpson and Chao1 index measurements

05_beta_diversity.ipynb - Community Ꞵ-diversity analysis described by Bray-Curtis and Jaccard distances.

06_logistic_regression_H_T1DM.R - Microbial logistic regression model comparison H vs T1DM confusuion matrix, ROC AUC curve and log odd coefficients

06_logistic_regression_H_T3cDM.R - Microbial logistic regression model comparison H vs T3cDM , ROC AUC curve and log odd coefficients

06_logistic_regression_T1DM_T3cDM.R - Microbial logistic regression model comparison T1D vs T3cDM, ROC AUC curve and log odd coefficients

07_violin_relative_abundance_top_features.R - Relative abundance boxplots top features from  microbial logistic regression model pairwise comparison

08_venn_diagram.R - Venn diagram using top features from  microbial logistic regression model pairwise comparison

./materials/08_micom.ipynb - Microbial community modelling 

./materials/09_growthrate_boxplot.R - Microbial community modelling  growth rate boxplot

./materials/010_metabolite_flux_boxplot_volcano.R - Microbial community modelling metabolite flux boxplot 

011_lefse.R -   Linear discriminant analysis Effect Size (LefSe) comparing disease groups  pairwise comparison

013_taxa_clinical_correlation.R   - Pearson correlation analyses between relative abundances of top genera identified from logistic models and selected clinical variables relevant to T1D and T3cDM.

014_linear_regression.R - Linear regression to evaluate the relationship between Escherichia-Shigella relative abundance and the binary clinical variable "Pancreatectomy"

015_metadata_analysis.R  - Linear mixed-effects models (GLMM) on clinical and laboratory variables collected at baseline and follow-up. Forest plot. 

016_lefse_sex.R  - Linear discriminant analysis Effect Size (LefSe) separated by disease comparing sex (female vs male)