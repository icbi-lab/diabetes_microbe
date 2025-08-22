
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



tax_table(phy) <- tax_fixed
ps_fam <- aggregate_taxa(phy, level = "Family")

top_families <- names(sort(taxa_sums(ps_fam), decreasing = TRUE))[1:20]



target_genera <- c("Enterobacteriaceae","Bacteroidaceae", "Ruminococcaceae", "Streptococcaceae","Rikenellaceae","Lachnospiraceae","Bifidobacteriaceae","Akkermansiaceae")

genus_taxa <- taxa_names(phy)[tax_table(phy)[, "Family"] %in% target_genera]
phy_target <- prune_taxa(genus_taxa, phy)

df_genus <- psmelt(phy_target) %>%
  #filter(Type %in% c("DM", "PDM")) %>%  # exclude K group if present
  mutate(Family = as.character(Family))


#####################################

stat_df <- df_genus %>%
  filter(Type %in% c("T1DM", "T3cDM", "H")) %>%
  group_by(Family) %>%
  summarise(
    T1DM_T3cDM = wilcox.test(Abundance[Type == "T1DM"], Abundance[Type == "T3cDM"])$p.value,
    T1DM_H     = wilcox.test(Abundance[Type == "T1DM"], Abundance[Type == "H"])$p.value,
    T3cDM_H    = wilcox.test(Abundance[Type == "T3cDM"], Abundance[Type == "H"])$p.value,
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("T1DM_") | starts_with("T3cDM_"),
    names_to = "comparison",
    values_to = "p.adj"
  ) %>%
  mutate(
    group1 = case_when(
      comparison == "T1DM_T3cDM" ~ "T1DM",
      comparison == "T1DM_H" ~ "T1DM",
      comparison == "T3cDM_H" ~ "T3cDM"
    ),
    group2 = case_when(
      comparison == "T1DM_T3cDM" ~ "T3cDM",
      comparison == "T1DM_H" ~ "H",
      comparison == "T3cDM_H" ~ "H"
    ),
    p.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

stat_df <- stat_df %>%
  group_by(Family) %>%
  mutate(y.position = max(df_genus$Abundance[df_genus$Family == unique(Family)], na.rm = TRUE) + 0.003 * row_number()) %>%
  ungroup()

stat_df$y.position <- pmin(stat_df$y.position, 0.049)
#write_csv(stat_df,"/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/stat_df_top_8_features_logreg.csv")
library(dplyr)

df_genus_clean <- df_genus %>%
  group_by(Family, Type) %>%
  mutate(
    Q1 = quantile(Abundance, 0.25, na.rm = TRUE),
    Q3 = quantile(Abundance, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 1.5 * IQR,
    upper = Q3 + 1.5 * IQR
  ) %>%
  filter(Abundance >= lower & Abundance <= upper) %>%
  ungroup()

df_genus_clean <- df_genus_clean %>%
  mutate(
    Type = case_when(
      str_starts(sample_information, "PDM") ~ "T3cDM",
      str_starts(sample_information, "DM")  ~ "T1DM",
      TRUE ~ "H"
    )
  )

# Plot
p <- ggplot(df_genus_clean, aes(x = Type, y = Abundance)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +

  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  facet_wrap(~ Family, scales = "free_y", ncol = 4) +
 # scale_fill_manual(values = c("T1DM" = "#E1812C", "T3cDM" = "#3A923A", "H" = "blue")) +
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
  )+
  coord_cartesian(ylim = c(0, 0.05))  


p
#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/boxplot_core_families.svg", height = 8, width = 18,dpi=300)
#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/boxplot_core_families.png", height = 8, width = 18,dpi=300)

#############################
