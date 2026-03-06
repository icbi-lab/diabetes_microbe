# Statistical Power calculation

# Load libraries
library(microbiome) 
library(phyloseq)
library(RColorBrewer) 
library(ggpubr) 
library(dplyr)
library(devtools)
library(micropower)
library(knitr)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(parallel)
library(vegan)
cores = detectCores()
set.seed(515087345)
PATH='data'

#Load data

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

# Set legend formatting
guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))

# Remove phylogenetic tree if it exists
ps1.com@phy_tree <- NULL

# Transform to relative abundance 
ps1.com <- microbiome::transform(ps1.com, "compositional")


# Build Group_new in metadata and store it in ps1
md <- data.frame(sample_data(ps1)) %>%
  mutate(Group_new = case_when(
    grepl("PDM", sample_information) ~ "T3cDM",
    grepl("DM",  sample_information) ~ "T1DM",
    grepl("K",   sample_information) ~ "H",
    TRUE ~ NA_character_
  ))

# Keep only samples with a group
keep <- rownames(md)[!is.na(md$Group_new)]
ps1g <- prune_samples(keep, ps1)

# Write Group_new back into phyloseq object
sample_data(ps1g)$Group_new <- md[keep, "Group_new"]

table(sample_data(ps1g)$Group_new)


# Define the simulation-based power function

power_permanova_ps <- function(ps,
                               group_var = "Group_new",
                               groups_keep = c("T1DM", "T3cDM"),
                               method = "bray",
                               n_per_group = 15,
                               R = 500,
                               alpha = 0.05,
                               permutations = 999,
                               seed = 1) {
  # Sets seed for reproducibility.
  
  set.seed(seed)
  
  # Keep only the groups of interest T1DM and T3cDM
  
  meta <- data.frame(sample_data(ps))
  meta <- meta[meta[[group_var]] %in% groups_keep, , drop = FALSE]
  
  # Check feasibility enough samples in each group to subsample
  
  tab <- table(meta[[group_var]])
  if (any(tab < n_per_group)) {
    stop(paste0("Not enough samples for n_per_group=", n_per_group,
                ". Available: ", paste(names(tab), tab, sep="=", collapse=", ")))
  }
  
 # TRUE/FALSE vector length R to store whether each simulation was significant.
  sig <- logical(R)
  
  # Loop over simulations
  
  for (i in seq_len(R)) {
    keep <- unlist(lapply(groups_keep, function(g) {
      samp <- rownames(meta)[meta[[group_var]] == g]
      sample(samp, n_per_group, replace = FALSE)
    }))
    
    ps_sub <- prune_samples(keep, ps)
    meta_sub <- data.frame(sample_data(ps_sub))
    
    # Subset phyloseq object containing only the sampled individuals
    d <- phyloseq::distance(ps_sub, method = method)
    
    # Distance matrix (Bray–Curtis by default) from the subset
    
    fit <- vegan::adonis2(as.dist(d) ~ meta_sub[[group_var]],
                          permutations = permutations)
    
    sig[i] <- fit$`Pr(>F)`[1] < alpha
  }
  # Group effect was significant at p < 0.05
  mean(sig)
}

# Run the power simulation for multiple sample sizes
ns <- c(6, 8, 10, 12, 14, 16, 17)


# Using 500 simulations
pow <- sapply(ns, function(n)
  power_permanova_ps(ps1g,
                     groups_keep = c("T1DM","T3cDM"),
                     n_per_group = n,
                     R = 500,          # increase to 1000 for smoother results
                     seed = 100 + n)
)


# Table of sample size per group vs estimated power
power_table <- data.frame(n_per_group = ns, power = pow)
power_table

