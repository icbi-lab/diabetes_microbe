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
library(ggplot2)
library(networkD3)


# In this code the disease categories used are the following
# PDM refers to T3cDM
# DM refers to T1DM 
# K refers to H

phy <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")
print(phy)

phy

phy <- subset_taxa(phy, !is.na(Phylum) & Phylum != "")

phy <- filter_taxa(phy, function(x) sum(x > 0) > 3, TRUE)

phy <- transform_sample_counts(phy, function(x) x / sum(x))
phy <- filter_taxa(phy, function(x) mean(x > 0.03) > 0.05, TRUE) 


sample_data(phy)$Type <- ifelse(grepl("PDM", sample_data(phy)$sample_information), "PDM",
                                ifelse(grepl("K", sample_data(phy)$sample_information), "K", "DM"))

tax_datatable(phy)

gp <- subset_taxa(phy, phy@tax_table[, 1] == "Bacteria")

p <- sankey_pq(gp, taxa = 1:6, fontSize = 20, nodeWidth = 20)

nodes <- p$x$nodes
links <- p$x$links

links$group <- nodes$name[links$source + 1]   
nodes$group <- nodes$name

colourJS <- 'd3.scaleOrdinal(d3.schemeCategory20);'


p_col <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value  = "value",
  NodeID = "name",
  NodeGroup = "group",
  LinkGroup = "group",
  fontSize  = 20,
  nodeWidth = 20,
  colourScale = colourJS
)

htmlwidgets::saveWidget(p_col, "/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/sankey_genus.html", selfcontained = TRUE)
