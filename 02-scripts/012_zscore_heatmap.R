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
library(tidymodels)
library(phyloseq)
library(pheatmap)
library(dplyr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)


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

tax <- as.data.frame(tax_table(phy))

tax <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_exact")]

colnames(tax)[colnames(tax) == "Species_exact"] <- "Species"

tax_fixed <- tax_table(as.matrix(tax))
rownames(tax_fixed) <- taxa_names(phy)

tax_table(phy) <- tax_fixed


selected_genera <- c("Blautia", "Subdoligranulum", "CAG-352", "Akkermansia", "Fusicatenibacter","Subdoligranulum", "Streptococcus", "Akkermansia", "Blautia", "CAG-352","Streptococcus", "Blautia", "Escherichia-Shigella", "Faecalibacterium", "CAG-352")


phy_genus <- tax_glom(phy, taxrank = "Genus")

genus_table <- tax_table(phy_genus)[, "Genus"]
taxa_to_keep <- taxa_names(phy_genus)[as.character(genus_table) %in% selected_genera]
phy_subset <- prune_taxa(taxa_to_keep, phy_genus)

otu_mat <- as(otu_table(phy_subset), "matrix")
if (!taxa_are_rows(phy_subset)) {
  otu_mat <- t(otu_mat)
}
otu_mat <- sweep(otu_mat, 2, colSums(otu_mat), FUN = "/")

genus_names <- as.character(tax_table(phy_subset)[, "Genus"])
rownames(otu_mat) <- genus_names

log_otu_mat <- log10(otu_mat + 1e-6)

sample_ids <- as.character(sample_data(phy_subset)$sample_information)
names(sample_ids) <- sample_names(phy_subset)
colnames(log_otu_mat) <- sample_ids[colnames(log_otu_mat)]

sample_info <- data.frame(sample_data(phy_subset))
sample_info$sample_id <- sample_info$sample_information
sample_info$Type <- ifelse(grepl("PDM", sample_info$sample_information), "PDM",
                           ifelse(grepl("K", sample_info$sample_information), "K", "DM"))

custom_order <- c("K1",  "K2",  "K3", "K4", "K5", "K6", "K7", "K8", "K9", "K10",
                  "DM1", "DM2", "DM3", "DM4", "DM5", "DM6", "DM7", "DM8", "DM9", "DM10", 
                  "DM11", "DM12", "DM13", "DM14", "DM15", "DM16", "DM17", "DM18", "DM19", "DM20", "DM21",
                  "PDM1", "PDM2", "PDM4", "PDM5", "PDM6", "PDM7", "PDM8", "PDM9", "PDM10", 
                  "PDM11", "PDM12", "PDM13", "PDM14", "PDM15", "PDM16", "PDM17", "PDM18")

sample_id_map <- setNames(sample_names(phy_subset), sample_info$sample_information)
custom_order_filtered <- custom_order[custom_order %in% colnames(log_otu_mat)]
log_otu_mat <- log_otu_mat[, custom_order_filtered]
annotation_col <- sample_info[sample_id_map[custom_order_filtered], , drop = FALSE]


############################################### Heatmap zscore per type
selected_genera <- unique(c(
  "Blautia", "Subdoligranulum", "CAG-352", "Akkermansia",
  "Fusicatenibacter", "Streptococcus", "Escherichia-Shigella",
  "Faecalibacterium"
))

otu_df <- as.data.frame(t(log_otu_mat))
otu_df$sample_information <- rownames(otu_df)
rownames(annotation_col) <- annotation_col$sample_information

log_otu_mat_t <- t(log_otu_mat) 

annotation_clean <- annotation_col[rownames(log_otu_mat_t), , drop = FALSE]
annotation_clean$sample_information <- rownames(annotation_clean)


otu_annotated <- left_join(otu_df, annotation_clean, by = "sample_information") %>%
  filter(!is.na(Type))

numeric_cols <- sapply(otu_annotated, is.numeric)
otu_clean <- otu_annotated[, numeric_cols]
otu_clean$Type <- otu_annotated$Type

group_means_df <- otu_clean %>%
  group_by(Type) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames("Type")

group_means <- t(group_means_df)
colnames(group_means) <- rownames(group_means_df)

group_means_selected <- group_means[rownames(group_means) %in% selected_genera, ]

group_means_z <- t(scale(t(group_means_selected)))

group_means_z <- as.matrix(group_means_z)
mode(group_means_z) <- "numeric"



ph <- pheatmap(t(group_means_z), 
               color = colorRampPalette(c("#998ec3", "white", "#f1a340"))(100),
               cluster_rows = TRUE,     
               cluster_cols = TRUE,     
               clustering_distance_rows = "euclidean",
               clustering_method = "complete", angle_col = 45,
               main = "")

otu_df <- as.data.frame(t(log_otu_mat))
otu_df$sample <- rownames(otu_df)

annotation_clean <- annotation_col
annotation_clean$sample <- rownames(annotation_clean)

otu_annotated <- left_join(otu_df, annotation_clean, by = "sample")

otu_annotated$Type <- factor(otu_annotated$Type, levels = c("K", "DM", "PDM"))

kruskal_results <- sapply(colnames(log_otu_mat_t), function(genus) {
  tryCatch({
    kruskal.test(otu_annotated[[genus]] ~ otu_annotated$Type)$p.value
  }, error = function(e) NA)
})

p_adj <- p.adjust(kruskal_results, method = "fdr")

kw_df <- data.frame(
  Genus = names(kruskal_results),
  P_value = kruskal_results,
  Adj_P_value = p_adj
)

significant_genera <- kw_df[kw_df$Adj_P_value < 0.05, ]
print(significant_genera)

get_stars <- function(p) {
  if (is.na(p)) return("")
  else if (p <= 0.001) return("***")
  else if (p <= 0.01) return("**")
  else if (p <= 0.05) return("*")
  else return("")
}

signif_genus <- kw_df$Genus
signif_stars <- sapply(kw_df$Adj_P_value, get_stars)
names(signif_stars) <- signif_genus

stars_matrix <- matrix("", nrow = nrow(group_means_z), ncol = ncol(group_means_z),
                       dimnames = dimnames(group_means_z))

for (genus in names(signif_stars)) {
  if (genus %in% rownames(stars_matrix)) {
    stars_matrix[genus, ] <- signif_stars[genus]
  }
}


ph <- pheatmap(t(group_means_z),                   
         color = colorRampPalette(c("#998ec3", "white", "#f1a340"))(100),
         cluster_rows = TRUE,                    
         cluster_cols = TRUE,                    
         display_numbers = t(stars_matrix),       
         number_color = "black",
         fontsize_number = 12,
         angle_col = 45)

#grid::grid.text(
#  "Z-score Normalized Heatmap Kruskal Wallis *p.adj<0.05 **p.adj<0.01",
#  x = 0.5, y = 0.97, gp = grid::gpar(fontface = "plain", fontsize = 14)
#)

ph
dev.off()

ggsave(plot=ph,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/zscore_heatmap.svg", height = 5, width = 10, dpi=300)
ggsave(plot=ph,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/zscore_heatmap.png", height = 5, width = 10, dpi=300)

