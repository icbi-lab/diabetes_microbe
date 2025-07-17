
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
library(tidymodels)
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

tax <- as.data.frame(tax_table(phy))

tax <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_exact")]

colnames(tax)[colnames(tax) == "Species_exact"] <- "Species"

tax_fixed <- tax_table(as.matrix(tax))
rownames(tax_fixed) <- taxa_names(phy)

tax_table(phy) <- tax_fixed



# Step 1: Extract relevant genus-level data
target_genera <- c("Fusicatenibacter", "Akkermansia", "Subdoligranulum",
                   "Blautia", "CAG-352", "Escherichia-Shigella",
                   "Faecalibacterium", "Streptococcus")

# Filter phyloseq object to only include those genera
genus_taxa <- taxa_names(phy)[tax_table(phy)[, "Genus"] %in% target_genera]
phy_target <- prune_taxa(genus_taxa, phy)

# Step 2: Melt to long format with sample, genus, and relative abundance
df_genus <- psmelt(phy_target) %>%
  #filter(Type %in% c("DM", "PDM")) %>%  # exclude K group if present
  mutate(Genus = as.character(Genus))


 


#############################################


# Recreate pairwise p-values with FDR adjustment
stat_df <- df_genus %>%
  filter(Type %in% c("DM", "PDM", "K")) %>%
  group_by(Genus) %>%
  summarise(
    DM_PDM = wilcox.test(Abundance[Type == "DM"], Abundance[Type == "PDM"])$p.value,
    DM_K   = wilcox.test(Abundance[Type == "DM"], Abundance[Type == "K"])$p.value,
    PDM_K  = wilcox.test(Abundance[Type == "PDM"], Abundance[Type == "K"])$p.value,
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("DM_") | starts_with("PDM_"),
    names_to = "comparison",
    values_to = "p.adj"
  ) %>%
  mutate(
    group1 = case_when(
      comparison == "DM_PDM" ~ "DM",
      comparison == "DM_K" ~ "DM",
      comparison == "PDM_K" ~ "PDM"
    ),
    group2 = case_when(
      comparison == "DM_PDM" ~ "PDM",
      comparison == "DM_K" ~ "K",
      comparison == "PDM_K" ~ "K"
    ),
    p.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Dynamically space brackets above each genus
stat_df <- stat_df %>%
  group_by(Genus) %>%
  mutate(y.position = max(df_genus$Abundance[df_genus$Genus == unique(Genus)], na.rm = TRUE) + 0.03 * row_number()) %>%
  ungroup()

p <- ggplot(df_genus, aes(x = Type, y = Abundance)) +
  geom_violin(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A", "K" = "blue")) +
  stat_pvalue_manual(
    stat_df,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 3
  ) +
  theme_minimal() +
  labs(
    title = "",
    y = "Relative Abundance",
    x = NULL
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

p



