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
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)



phy <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")

df <- read_csv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/PDM merged 3.0_modified.csv")
library(dplyr)

df <- rename(df, sample_information = Probennummer)

samdf <- as.data.frame(sample_data(phy))
samdf$sex <- NA


for (i in 1:nrow(samdf)) {
  sample_id <- samdf$sample_information[i]
  

  match_index <- which(df$sample_information == sample_id)
  
  
  if (length(match_index) == 1) {
    samdf$sex[i] <- df$sex[match_index]
  }
}


samdf$disease <- ifelse(grepl("PDM", samdf$sample_information), "T3cDM",
                        ifelse(grepl("K", samdf$sample_information), "H", "T1D"))


sample_data(phy) <- samdf


phy_T1D <- subset_samples(phy, disease == "T1D")


phy_T3cDM <- subset_samples(phy, disease == "T3cDM")

phy_H <- subset_samples(phy, disease == "H")





phy <- phy_T1D

colnames(tax_table(phy))[colnames(tax_table(phy)) == "Species_exact"] <- "Species"


tax <- tax_table(phy)

tax_clean <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]


tax_table(phy) <- tax_clean
mm_lefse_T1D <- run_lefse(
  phy,
  wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  group = "sex",
  multigrp_strat = TRUE,  taxa_rank = "Genus"
)


get_lefse_df <- function(lef_out, comp_name, group1, group2) {
  df <- marker_table(lef_out) %>%
    data.frame() %>%
    mutate(
      feature_mod = stringr::str_extract(feature, "[^|]+$"),
      signed_lda = ifelse(enrich_group == group1, -ef_lda, ef_lda),
      comparison = comp_name
    )
}


dat_T1D <- get_lefse_df(mm_lefse_T1D, "f vs m", "f", "m")


#


phy <- phy_T3cDM

colnames(tax_table(phy))[colnames(tax_table(phy)) == "Species_exact"] <- "Species"


tax <- tax_table(phy)


tax_clean <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

tax_table(phy) <- tax_clean
mm_lefse_T3cDM <- run_lefse(
  phy,
  wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  group = "sex",
  multigrp_strat = TRUE,  taxa_rank = "Genus"
)


#

phy <- phy_H

colnames(tax_table(phy))[colnames(tax_table(phy)) == "Species_exact"] <- "Species"

tax <- tax_table(phy)

tax_clean <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]


tax_table(phy) <- tax_clean


mm_lefse_H <- run_lefse(
  phy,
  wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  group = "sex",
  multigrp_strat = TRUE,  taxa_rank = "Genus"
)

####################################################

library(dplyr)
library(ggplot2)

# Add a 'condition' column to each
dat_T3cDM <- get_lefse_df(mm_lefse_T3cDM, "f vs m", "f", "m") %>%
  mutate(condition = "T3cDM")

dat_T1D <- get_lefse_df(mm_lefse_T1D, "f vs m", "f", "m") %>%
  mutate(condition = "T1D")

dat_H <- get_lefse_df(mm_lefse_H, "f vs m", "f", "m") %>%
  mutate(condition = "H")


# Combine both
dat_all <- bind_rows(dat_T3cDM, dat_T1D,dat_H)

# Create the combined plot
p_all <- ggplot(dat_all, aes(x = feature, y = signed_lda, fill = enrich_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_y_continuous(name = "LDA SCORE (log10)                  p.adj < 0.05") +
  scale_fill_manual(values = c("f" = "#ef8a62", "m" = "#67a9cf")) +
  scale_x_discrete(labels = dat_all$feature_mod) +
  facet_wrap(~ condition, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank()
  )


p_all


#write_csv(dat_all, "/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/dat_all_barplot_lefse_sex_comparison.csv")

#ggsave(plot=p_all,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/barplot_lefse_sex_comparison.svg", height = 10, width = 10, dpi=300)
#ggsave(plot=p_all,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/barplot_lefse_sex_comparison.png", height = 10, width = 10, dpi=300)


