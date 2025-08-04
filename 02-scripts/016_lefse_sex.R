library(tidyverse)

library(phyloseq)

library(microbiomeMarker)

library(rsample)
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

#library(tidymodels)


library(dplyr)
library(tidyr)
library(purrr)

library(ggplot2)


# In this code the disease categories used are the following
# PDM refers to T3cDM
# DM refers to T1DM 
# K refers to H

phy <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")

df <- read_csv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/PDM merged 3.0_modified.csv")
df <- df %>%
  rename(sample_information = Probennummer)

samdf <- as.data.frame(sample_data(phy))
samdf$sex <- NA


for (i in 1:nrow(samdf)) {
  sample_id <- samdf$sample_information[i]
  

  match_index <- which(df$sample_information == sample_id)
  
  
  if (length(match_index) == 1) {
    samdf$sex[i] <- df$sex[match_index]
  }
}




# Create the 'disease' column based on 'sample_information'
samdf$disease <- ifelse(grepl("PDM", samdf$sample_information), "T3cDM",
                        ifelse(grepl("K", samdf$sample_information), "H", "T1D"))


sample_data(phy) <- samdf

# Filter for T1D
phy_T1D <- subset_samples(phy, disease == "T1D")

# Filter for T3cDM 
phy_T3cDM <- subset_samples(phy, disease == "T3cDM")

############################################################################## LEFSE M vs F T1D 
phy <- phy_T1D

print(phy)

phy

phy <- subset_taxa(phy, !is.na(Phylum) & Phylum != "")

phy <- filter_taxa(phy, function(x) sum(x > 0) > 3, TRUE)

phy <- transform_sample_counts(phy, function(x) x / sum(x))
phy <- filter_taxa(phy, function(x) mean(x > 0.03) > 0.05, TRUE) # adjust 0.03 as needed




tax <- as.data.frame(tax_table(phy))

tax <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_exact")]

colnames(tax)[colnames(tax) == "Species_exact"] <- "Species"

tax_fixed <- tax_table(as.matrix(tax))
rownames(tax_fixed) <- taxa_names(phy)

tax_table(phy) <- tax_fixed


lef<- run_lefse(phy, group = "sex", norm = "CPM", kw_cutoff = 0.1, taxa_rank = "Genus",lda_cutoff = 1)

get_lefse_df <- function(lef_out, comp_name, group1, group2) {
  df <- marker_table(lef_out) %>%
    data.frame()
  
  # figure out what the correct column name is
  feature_col <- names(df)[grep("feature|taxa|taxon", names(df), ignore.case = TRUE)][1]
  
  df <- df %>%
    mutate(
      feature_mod = stringr::str_extract(.[[feature_col]], "[^|]+$"),
      signed_lda = ifelse(enrich_group == group1, -ef_lda, ef_lda),
      comparison = comp_name
    )
  
  return(df)
}

dat <- get_lefse_df(lef, "m vs f", "m", "f")



dat_all <- dat_all %>%
  group_by(comparison) %>%
  arrange(enrich_group, desc(signed_lda), .by_group = TRUE) %>%
  mutate(feature_mod = factor(feature_mod, levels = rev(unique(feature_mod))))

dat_all$enrich_group <- factor(dat_all$enrich_group, levels = c("DM", "PDM", "K"))
dat_all$enrich_group <- factor(dat_all$enrich_group, levels = c("DM", "PDM", "K"))

dat_all <- dat_all %>%
  group_by(comparison) %>%
  arrange(enrich_group, desc(signed_lda), .by_group = TRUE) %>%
  mutate(feature_mod = factor(feature_mod, levels = rev(unique(feature_mod)))) %>%
  ungroup()

dat_all$comparison <- factor(
  dat_all$comparison,
  levels = c("PDM vs DM", "DM vs K", "PDM vs K")
)

dat_all <- dat_all %>%
  mutate(feature_id = paste(comparison, feature_mod, sep = " | "))

dat_all <- dat_all %>%
  group_by(comparison) %>%
  arrange(enrich_group, desc(signed_lda), .by_group = TRUE) %>%
  mutate(feature_id = factor(feature_id, levels = rev(unique(feature_id)))) %>%
  ungroup()

dat_all$comparison <- factor(
  dat_all$comparison,
  levels = c("PDM vs DM", "DM vs K", "PDM vs K")
)

dat_all <- dat_all %>%
  mutate(feature_id = paste(comparison, feature_mod, sep = " | "))

dat_all <- dat_all %>%
  group_by(comparison) %>%
  arrange(enrich_group, desc(signed_lda), .by_group = TRUE) %>%
  mutate(feature_id = factor(feature_id, levels = rev(unique(feature_id)))) %>%
  ungroup()

legend_df <- data.frame(
  feature_id = NA,
  signed_lda = NA,
  enrich_group = NA,
  padj_label = "p.adj < 0.05"
)

# Plot
p_all <- ggplot(dat_all, aes(x = feature_id, y = signed_lda, fill = enrich_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_y_continuous(name = "LDA SCORE (log10)                  p.adj < 0.05") +
  scale_fill_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A", "K" = "#3274A1")) +
  scale_x_discrete(labels = dat_all$feature_mod) +   # 👈 custom axis labels
  facet_wrap(~ comparison, scales = "free_y", ncol = 1) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank()
  )

p_all



#ggsave(plot=p_all,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/barplot_lefse_all_comparison.svg", height = 10, width = 10, dpi=300)
#ggsave(plot=p_all,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/barplot_lefse_all_comparison.png", height = 10, width = 10, dpi=300)

