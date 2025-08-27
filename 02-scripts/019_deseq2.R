library(tidyverse)

library(phyloseq)

library(microbiomeMarker)

library(rsample)
library(janitor)
library(dplyr)
library(microbiomeMarker)



library(dplyr)
library(tidyr)
library(stringr)
library(viridis)

library(readr)
library(themis)
library(MiscMetabar)



library(dplyr)
library(tidyr)
library(purrr)
library("DESeq2")
library(ggplot2)


# In this code the disease categories used are the following
# PDM refers to T3cDM
# DM refers to T1DM 
# K refers to H
library("mia")
phy <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")
print(phy)

phy

phy <- subset_taxa(phy, !is.na(Phylum) & Phylum != "")

phy <- filter_taxa(phy, function(x) sum(x > 0) > 3, TRUE)

phy <- transform_sample_counts(phy, function(x) x / sum(x))
phy <- filter_taxa(phy, function(x) mean(x > 0.03) > 0.05, TRUE) # adjust 0.03 as needed

# Create "Type" column based on pattern in "sample_information"
sample_data(phy)$Type <- ifelse(grepl("PDM", sample_data(phy)$sample_information), "PDM",
                                ifelse(grepl("K", sample_data(phy)$sample_information), "K", "DM"))


phy <- subset_samples(phy, Type %in% c("DM", "PDM"))

phy_fam <- aggregate_taxa(phy, level = "Family")
