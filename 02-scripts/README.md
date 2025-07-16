# Preprocessing & analysis

##  Preprocessing: nf-core/ampliseq
The preprocessing was performed using nf-core/ampliseq: Amplicon sequencing analysis workflow using DADA2 and QIIME2. 

Check ./nf-core_ampliseq/run_nf_core_ampliseq.slurm for more information. 

Samplesheet available at diabetes_microbe/01-tables/nf-core_ampliseq

##  Analysis

01_compositional.R  - The stacked bar plot depicts the relative abundance for the 10 most abundant families in all samples
02_upsetplot.ipynb - UpSetPlot visualizing the set overlap for all the recorded families 

