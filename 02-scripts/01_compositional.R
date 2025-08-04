# Compositional barplot

# Load libraries
library(microbiome) 
library(phyloseq)
library(RColorBrewer) 
library(ggpubr) 
library(dplyr) 


#Load data
#ps1 <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")
ps1 <- readRDS("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/dada2_phyloseq.rds")
print(ps1)
ps1.com <- ps1

#Rename taya as ASVs and extract taxonomy table
taxa_names(ps1.com) <- paste0("ASV_", rownames(tax_table(ps1.com)))
taxic <- as.data.frame(ps1.com@tax_table)
taxic$OTU <- rownames(taxic) 
colnames(taxic)

#Update and clean tayonomy table
taxmat <- as.matrix(taxic) 
new.tax <- tax_table(taxmat)
tax_table(ps1.com) <- new.tax 
tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "Unclassified family"

# set legend formatting
guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))

#Remove phylogenetic tree if it exists
ps1.com@phy_tree <- NULL

# Transform to relative abundance 
ps1.com <- microbiome::transform(ps1.com, "compositional")


# Filter out taxa with mean relative abundance < 0.1% 
ps1.com <- filter_taxa(ps1.com, function(x) mean(x) > 0.001, prune = TRUE)

#Agregate taxa at Family level and select top 10 families
ps_fam <- aggregate_taxa(ps1.com, level = "Family")
top_families <- names(sort(taxa_sums(ps_fam), decreasing = TRUE))[1:10]
ps1.com.fam <- prune_taxa(top_families, ps_fam)

#Relative abundance transformation
ps1.com.fam.rel <- microbiome::transform(ps1.com.fam, "compositional")

#Plot
plot.composition.relAbun <- plot_composition(ps1.com.fam.rel,
                                             sample.sort = "scientific_name",
                                             x.label = "env_material") 

plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
plot.composition.relAbun <- plot.composition.relAbun + scale_fill_brewer("Family", palette = "Paired") + theme_bw() 
plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))
print(plot.composition.relAbun)



# Reorder samples based on clnical information
ps_fam <- aggregate_taxa(ps1.com, level = "Family")
top_families <- names(sort(taxa_sums(ps_fam), decreasing = TRUE))[1:10]
ps1.com.fam <- prune_taxa(top_families, ps_fam)

ps1.com.fam.rel <- microbiome::transform(ps1.com.fam, "compositional")
custom_order <- c("K1",  "K2",  "K3", "K4", "K5",   "K6",  "K7",  "K8",   "K9", "K10",
                  "DM1", "DM2",  "DM3", "DM4", "DM5",  "DM6", "DM7", "DM8", "DM9", "DM10", 
                  "DM11", "DM12", "DM13", "DM14", "DM15", "DM16", "DM17", "DM18", "DM19", "DM20", "DM21",
                  "PDM1", "PDM2", "PDM4", "PDM5", "PDM6", "PDM7", "PDM8", "PDM9", "PDM10", 
                  "PDM11", "PDM12", "PDM13", "PDM14", "PDM15", "PDM16", "PDM17", "PDM18")

sample_data(ps1.com.fam.rel)$sample_information <- factor(
  sample_data(ps1.com)$sample_information,
  levels = custom_order
)

# Plot
plot.composition.relAbun <- plot_composition(ps1.com.fam.rel,
                                             sample.sort = "sample_information", 
                                             x.label = "sample_information")      

# Add aesthetics
plot.composition.relAbun <- plot.composition.relAbun +
  theme(legend.position = "bottom") +
  scale_fill_brewer("Family", palette = "Paired") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,size = 15),    axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +

  labs(y = "Relative abundance") + 
  guide_italics +
  theme(legend.title = element_text(size = 18))

# Print
print(plot.composition.relAbun)


#write.csv(plot.composition.relAbun[["data"]], "/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/compositional_plots_family.csv", row.names = FALSE)


#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/compositional_plots_family.svg", height = 10, width = 10, dpi=300)
#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/compositional_plots_family.png", height = 10, width = 10, dpi=300)



