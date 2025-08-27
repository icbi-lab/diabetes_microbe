
#Load libraries 
library(tidyverse)
library(phyloseq)
library(microbiomeMarker)
library(janitor)
library(dplyr)
library(microbiomeMarker)
#library(knitr)
library(dplyr)
library(tidyr)
library(stringr)
library(viridis)
library(readr)
library(themis)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)


phy <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")
print(phy)
phy


phy <- subset_taxa(phy, !is.na(Phylum) & Phylum != "")




phy <- filter_taxa(phy, function(x) sum(x > 0) > 3, TRUE)


phy <- transform_sample_counts(phy, function(x) x / sum(x))
phy <- filter_taxa(phy, function(x) mean(x > 0.03) > 0.05, TRUE) 

sample_data(phy)$Type <- ifelse(grepl("PDM", sample_data(phy)$sample_information), "PDM",
                                ifelse(grepl("K", sample_data(phy)$sample_information), "K", "DM"))

sample_data(phy)$Type <- ifelse(grepl("PDM", sample_data(phy)$sample_information), "T3cDM",
                                ifelse(grepl("K", sample_data(phy)$sample_information), "H", "T1DM"))
phy <- microbiome::transform(phy, "compositional")

tax <- as.data.frame(tax_table(phy))

tax <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_exact")]

colnames(tax)[colnames(tax) == "Species_exact"] <- "Species"

tax_fixed <- tax_table(as.matrix(tax))
rownames(tax_fixed) <- taxa_names(phy)

###############
#HEATMAP prevalence 0.3

library(phyloseq)
library(dplyr)
library(ggplot2)

# Your phyloseq object: phy
stopifnot("Type" %in% colnames(sample_data(phy)))

# Prevalence across all samples
otu_mat <- as(otu_table(phy), "matrix")
if (!taxa_are_rows(phy)) otu_mat <- t(otu_mat)

prev <- rowMeans(otu_mat > 0)                 # proportion of samples with nonzero abundance
core_taxa <- names(prev[prev >= 0.5])

phy_core <- prune_taxa(core_taxa, phy)

# Change "Genus" to your desired rank (e.g., "Family", "Phylum")
phy_core <- tax_glom(phy_core, taxrank = "Family", NArm = TRUE)

phy_core_rel <- transform_sample_counts(phy_core, function(x) x / sum(x))


df <- psmelt(phy_core_rel) %>%
  mutate(
    Taxon = if ("Family" %in% colnames(.)) as.character(Family) else as.character(OTU),
    logAbundance = log1p(Abundance)  
  )


# Order taxa by decreasing overall prevalence (or mean abundance)
taxa_order <- df %>%
  group_by(Taxon) %>%
  summarize(prev = mean(Abundance > 0), mean_abund = mean(Abundance)) %>%
  arrange(desc(prev), desc(mean_abund)) %>%
  pull(Taxon)

df$Taxon <- factor(df$Taxon, levels = taxa_order)

# Order samples within each Type by their total abundance of core taxa (nice for readability)
sample_order <- df %>%
  group_by(Sample, Type) %>%
  summarize(total = sum(Abundance), .groups = "drop") %>%
  arrange(Type, desc(total)) %>%
  pull(Sample)

df$Sample <- factor(df$Sample, levels = sample_order)



p <- ggplot(df, aes(x = Sample, y = Family, fill = logAbundance)) +
  geom_tile() +
  facet_wrap(~ Type, scales = "free_x") +
  scale_fill_gradient2(
    low = "#ffffbf", mid = "#ffffbf", high = "#fc8d59",
    midpoint = median(df$logAbundance, na.rm = TRUE),  # where white sits
    name = "log(R.abundance+1)"
  ) +
  labs(x = NULL, y = NULL, title = "") +
  theme_minimal(base_size = 
                  25) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.1, "lines")
  )


p

ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/heatmap_core_micrbiome.svg", height = 5, width = 18,dpi=300)
ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/heatmap_core_micriobiome.png", height =5, width = 18,dpi=300)

