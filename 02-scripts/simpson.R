library(tidyverse)

library(phyloseq)

library(ggplot2)
library(ggpmisc)

library(dplyr)

library(microbiomeMarker)




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


phy <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")
print(phy)

phy

# Remove NA phylum
phy <- subset_taxa(phy, !is.na(Phylum) & Phylum != "")

# Minimum prevalence filter: keep features in >3 samples
phy <- filter_taxa(phy, function(x) sum(x > 0) > 3, TRUE)

# Abundance threshold filter (e.g., 3%)
phy <- transform_sample_counts(phy, function(x) x / sum(x))
phy <- filter_taxa(phy, function(x) mean(x > 0.03) > 0.05, TRUE) # adjust 0.03 as needed

# Create "Type" column based on pattern in "sample_information"
sample_data(phy)$Type <- ifelse(grepl("PDM", sample_data(phy)$sample_information), "PDM",
                                ifelse(grepl("K", sample_data(phy)$sample_information), "K", "DM"))


tax <- as.data.frame(tax_table(phy))

# Keep only the 7 standard taxonomic ranks
tax <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_exact")]

# Rename 'Species_exact' to 'Species'
colnames(tax)[colnames(tax) == "Species_exact"] <- "Species"

# Reassign as tax_table and preserve rownames
tax_fixed <- tax_table(as.matrix(tax))
rownames(tax_fixed) <- taxa_names(phy)

# Set it back into phyloseq object
tax_table(phy) <- tax_fixed


df <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/PDM merged 3.0_modified.csv")%>%clean_names()


name_map_full <- c(
  "probennummer" = "sample_information",
  "verstorben" = "Deceased",
  "pankreatektomie" = "Pancreatectomy",
  "c2" = "C-Peptide (FU)",
  "nikotin" = "Smoking Status",
  "sex" = "Sex",
  "age" = "Age",
  
  # Baseline (BS)
  "hyperlipid_mie_20181" = "Hyperlipidemia (BS)",
  "arterielle_hyperotnie1" = "Arterial Hypertension (BS)",
  "bmi1" = "Body Mass Index (BMI) (BS)",
  "gr_e1" = "Height (BS)",
  "gewicht1" = "Weight (BS)",
  "lipidsenker1" = "Lipid-lowering Medication (BS)",
  "art_von_lipidsenker1" = "Type of Lipid-lowering Medication (BS)",
  "rr_medikation1" = "Blood Pressure Medication (BS)",
  "b_blocker1" = "Beta Blockers (BS)",
  "ace_hemmer1" = "ACE Inhibitors (BS)",
  "diuretika1" = "Diuretics (BS)",
  "insulin1" = "Insulin Therapy (BS)",
  "langzeit_insulin1" = "Long-acting Insulin (BS)",
  "kurzzeit_insulin1" = "Short-acting Insulin (BS)",
  "misch_inuslin1" = "Mixed Insulin (BS)",
  "masld1" = "MASLD (BS)",
  "khk1" = "Coronary Heart Disease (CHD) (BS)",
  "ca1" = "Cancer (unspecified) (BS)",
  
  # Follow-up (FU)
  "hyperlipid_mie_20182" = "Hyperlipidemia (FU)",
  "arterielle_hyperotnie2" = "Arterial Hypertension (FU)",
  "bmi2" = "Body Mass Index (BMI) (FU)",
  "gr_e2" = "Height (FU)",
  "gewicht2" = "Weight (FU)",
  "lipidsenker2" = "Lipid-lowering Medication (FU)",
  "art_von_lipidsenker2" = "Type of Lipid-lowering Medication (FU)",
  "rr_medikation2" = "Blood Pressure Medication (FU)",
  "b_blocker2" = "Beta Blockers (FU)",
  "ace_hemmer2" = "ACE Inhibitors (FU)",
  "diuretika2" = "Diuretics (FU)",
  "insulin2" = "Insulin Therapy (FU)",
  "langzeit_insulin2" = "Long-acting Insulin (FU)",
  "kurzzeit_insulin2" = "Short-acting Insulin (FU)",
  "misch_inuslin2" = "Mixed Insulin (FU)",
  "masld2" = "MASLD (FU)",
  "khk2" = "Coronary Heart Disease (CHD) (FU)",
  "ca2" = "Cancer (unspecified) (FU)",
  
  # Lab values (BS)
  "leukozyten1" = "Leukocytes (BS)",
  "h_moglobin1" = "Hemoglobin (BS)",
  "h_matokrit1" = "Hematocrit (BS)",
  "thrombozyten1" = "Platelets (BS)",
  "harnstoff1" = "Urea (BS)",
  "creatinin_enzym_idms_1" = "Creatinine (IDMS) (BS)",
  "glomerul_re_filtrationsrate1" = "Glomerular Filtration Rate (GFR) (BS)",
  "bilirubin_gesamt1" = "Total Bilirubin (BS)",
  "natrium1" = "Sodium (BS)",
  "got_asat_1" = "AST (BS)",
  "gpt_alat_1" = "ALT (BS)",
  "gamma_gt1" = "GGT (BS)",
  "alkalische_phosphatase1" = "Alkaline Phosphatase (BS)",
  "lactat_dehydrogenase_ldh_1" = "LDH (BS)",
  "c_reaktives_prot_crp_1" = "CRP (BS)",
  "quicktest_pt_1" = "Prothrombin Time (BS)",
  "inr_pt_1" = "INR (BS)",
  "part_thrombopl_zeit_a_ptt_1" = "aPTT (BS)",
  "fibrinogen_funkt_n_clauss1" = "Fibrinogen (Clauss) (BS)",
  "albumin1" = "Albumin (BS)",
  "glukose1" = "Glucose (BS)",
  "hb_a1c_dcct_ngsp_1" = "HbA1c (DCCT/NGSP) (BS)",
  "hb_a1c_ifcc_1" = "HbA1c (IFCC) (BS)",
  "cholesterin1" = "Total Cholesterol (BS)",
  "non_hdl_cholesterin1" = "Non-HDL Cholesterol (BS)",
  "triglyceride1" = "Triglycerides (BS)",
  "hdl_cholesterin1" = "HDL Cholesterol (BS)",
  "ldl_cholesterin1" = "LDL Cholesterol (BS)",
  "eisen1" = "Iron (BS)",
  "ferritin1" = "Ferritin (BS)",
  "transferrin1" = "Transferrin (BS)",
  "transferrins_ttigung1" = "Transferrin Saturation (BS)",
  "troponin_t_hoch_sens_1" = "Troponin T (High Sensitivity) (BS)",
  "nt_pro_bnp1" = "NT-proBNP (BS)",
  
  # Lab values (FU)
  "leukozyten2" = "Leukocytes (FU)",
  "h_moglobin2" = "Hemoglobin (FU)",
  "h_matokrit2" = "Hematocrit (FU)",
  "thrombozyten2" = "Platelets (FU)",
  "harnstoff2" = "Urea (FU)",
  "creatinin_enzym_idms_2" = "Creatinine (IDMS) (FU)",
  "glomerul_re_filtrationsrate2" = "Glomerular Filtration Rate (GFR) (FU)",
  "bilirubin_gesamt2" = "Total Bilirubin (FU)",
  "natrium2" = "Sodium (FU)",
  "got_asat_2" = "AST (FU)",
  "gpt_alat_2" = "ALT (FU)",
  "gamma_gt2" = "GGT (FU)",
  "alkalische_phosphatase2" = "Alkaline Phosphatase (FU)",
  "lactat_dehydrogenase_ldh_2" = "LDH (FU)",
  "c_reaktives_prot_crp_2" = "CRP (FU)",
  "quicktest_pt_2" = "Prothrombin Time (FU)",
  "inr_pt_2" = "INR (FU)",
  "part_thrombopl_zeit_a_ptt_2" = "aPTT (FU)",
  "fibrinogen_funkt_n_clauss2" = "Fibrinogen (Clauss) (FU)",
  "albumin2" = "Albumin (FU)",
  "glukose2" = "Glucose (FU)",
  "hb_a1c_dcct_ngsp_2" = "HbA1c (DCCT/NGSP) (FU)",
  "hb_a1c_ifcc_2" = "HbA1c (IFCC) (FU)",
  "cholesterin2" = "Total Cholesterol (FU)",
  "non_hdl_cholesterin2" = "Non-HDL Cholesterol (FU)",
  "triglyceride2" = "Triglycerides (FU)",
  "hdl_cholesterin2" = "HDL Cholesterol (FU)",
  "ldl_cholesterin2" = "LDL Cholesterol (FU)",
  "eisen2" = "Iron (FU)",
  "ferritin2" = "Ferritin (FU)",
  "transferrin2" = "Transferrin (FU)",
  "transferrins_ttigung2" = "Transferrin Saturation (FU)",
  "troponin_t_hoch_sens_2" = "Troponin T (High Sensitivity) (FU)",
  "nt_pro_bnp2" = "NT-proBNP (FU)"
)


df <- df %>%
  dplyr::rename_with(.cols = all_of(names(name_map_full)), .fn = ~ name_map_full[.x])

sample_data(phy) <- sample_data(phy)[, "sample_information", drop = FALSE]

sam <- as(sample_data(phy), "data.frame")

sam_new <- dplyr::left_join(sam, df, by = "sample_information")


rownames(sam_new) <- rownames(sam)  


sample_data(phy) <- phyloseq::sample_data(sam_new)


otu_df <- as.data.frame(otu_table(phy))

# Create "Type" column based on pattern in "sample_information"
sample_data(phy)$Type <- ifelse(grepl("PDM", sample_data(phy)$sample_information), "PDM",
                                ifelse(grepl("K", sample_data(phy)$sample_information), "K", "DM"))


# Packages
library(phyloseq)
library(vegan)   # for diversity()
library(dplyr)

# 1) Get the abundance matrix from your phyloseq object
otu_mat <- as(otu_table(phy), "matrix")

# Make samples = rows, taxa = columns (as vegan expects)
if (taxa_are_rows(phy)) {
  otu_mat <- t(otu_mat)
}

# 2) Compute Simpson diversity per sample
# Note: vegan::diversity(index = "simpson") returns the Gini–Simpson index = 1 - D
#simpson <- vegan::diversity(otu_mat, index = "simpson")
inv_simpson <- vegan::diversity(otu_mat, index = "invsimpson")


# 3) Build a data frame with sample Type
sdat <- as.data.frame(sample_data(phy))
# Ensure order aligns with otu_mat rownames
sdat <- sdat[rownames(otu_mat), , drop = FALSE]

alpha_df <- data.frame(
  sample  = rownames(otu_mat),
  Type    = sdat$Type,
  inv_simpson = inv_simpson,
  row.names = NULL
)

# 4) Summarize by Type (mean, sd, median, IQR, n)
by_type <- alpha_df %>%
  group_by(Type) %>%
  summarise(
    n      = dplyr::n(),
    mean   = mean(inv_simpson, na.rm = TRUE),
    sd     = sd(inv_simpson, na.rm = TRUE),
    median = median(inv_simpson, na.rm = TRUE),
    IQR    = IQR(inv_simpson, na.rm = TRUE),
    .groups = "drop"
  )

by_type


alpha_df <- alpha_df %>%
  mutate(Type = case_when(
    Type == "PDM" ~ "T3cDM",
    Type == "DM"  ~ "T1DM",
    Type == "K"   ~ "H",
    TRUE ~ Type
  ))

ggplot(alpha_df, aes(x = Type, y = inv_simpson, fill = Type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Inverse Simpson Diversity Index by Group",
    x = "",
    y = "Inverse Simpson (1 / D)"
  ) +
  scale_fill_manual(values = c("T1DM" = "#E1812C", "T3cDM" = "#3A923A", "H" = "#3C91E6")) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

write_csv(alpha_df, "/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/simpson_calculated_alpha_div.csv")

