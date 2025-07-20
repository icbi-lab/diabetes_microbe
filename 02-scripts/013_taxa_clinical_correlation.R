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


##################################### NEW 
# Sample metadata (clinical features)
clinical_df <- data.frame(sample_data(phy), check.names = FALSE)
clinical_df$Sex_ <- ifelse(clinical_df$Sex == "f", 0, 1)

clinical_vars <- c("Smoking Status","Insulin Therapy (BS)","Short-acting Insulin (BS)","Mixed Insulin (BS)","Long-acting Insulin (BS)","Pancreatectomy","HbA1c (DCCT/NGSP) (BS)", "MASLD (BS)","Age","Sex_"

)

abund <- as.data.frame(otu_table(phy))
if (taxa_are_rows(phy)) {
  abund <- t(abund)
}
all(rownames(clinical_df) == rownames(abund))  # should be TRUE

shared_samples <- intersect(rownames(clinical_df), rownames(abund))
clinical_df <- clinical_df[shared_samples, , drop = FALSE]
abund <- abund[shared_samples, , drop = FALSE]

tax <- as.data.frame(tax_table(phy))
tax$tax_id <- rownames(tax)

tax_long <- tax %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = Genus, names_to = "rank", values_to = "name") %>%
  filter(!is.na(name), name != "") %>%
  mutate(feature_mod = paste0(tolower(substr(rank, 1, 1)), "__", name)) %>%
  dplyr::select(tax_id, feature_mod)



abund <- as.data.frame(otu_table(phy))
if (taxa_are_rows(phy)) {
  abund <- t(abund)
}
abund$tax_id <- colnames(otu_table(phy))
abund_tidy <- as.data.frame(t(otu_table(phy)))
abund_tidy$sample_id <- rownames(abund_tidy)


abund_long <- abund_tidy %>%
  pivot_longer(-sample_id, names_to = "tax_id", values_to = "abundance")

abund_long <- abund_long %>%
  left_join(tax_long, by = "tax_id") %>%
  filter(!is.na(feature_mod))

abund_grouped <- abund_long %>%
  group_by(sample_id, feature_mod) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

abund_wide <- abund_grouped %>%
  pivot_wider(names_from = feature_mod, values_from = abundance)

abund_matrix <- as.data.frame(abund_wide)
rownames(abund_matrix) <- abund_matrix$sample_id
abund_matrix$sample_id <- NULL

#features_to_use <- unique(dat_all$feature_mod)
#features_to_use <- paste0("g__", features_to_use)

features_to_use <- c(
  "g__Akkermansia",
  "g__CAG-352",
  "g__Escherichia-Shigella",
  "g__Faecalibacterium",
  "g__Blautia",
  "g__Streptococcus",
  "g__Subdoligranulum",
  "g__Fusicatenibacter"
)

abund_sel <- abund_matrix[, features_to_use]


abund_sel <- abund_matrix[, features_to_use, drop = FALSE]

clinical_df <- data.frame(sample_data(phy), check.names = FALSE)

shared_samples <- intersect(rownames(clinical_df), rownames(abund_sel))
clinical_df <- clinical_df[shared_samples, ]
abund_sel <- abund_sel[shared_samples, ]

clinical_vars <- clinical_vars[clinical_vars %in% colnames(clinical_df)]

clinical_vars <- clinical_vars[clinical_vars %in% names(clinical_df[sapply(clinical_df, is.numeric)])]

desired_order <- c(  "Age",
                     "Sex_",
  "Smoking Status",
  "HbA1c (DCCT/NGSP) (BS)",
  "Insulin Therapy (BS)",
  "Short-acting Insulin (BS)",
  "Mixed Insulin (BS)",
  "Long-acting Insulin (BS)",
  "Pancreatectomy",
  "MASLD (BS)"

)


cor_results <- expand.grid(
  microbe = colnames(abund_sel),
  clinical_var = clinical_vars,
  stringsAsFactors = FALSE
) %>%
  mutate(
    cor = purrr::map2_dbl(microbe, clinical_var, ~{
      x <- abund_sel[[.x]]
      y <- clinical_df[[.y]]
      if (length(unique(na.omit(x))) < 2 || length(unique(na.omit(y))) < 2) return(NA)
      cor.test(x, y, method = "pearson")$estimate
    }),
    p_value = purrr::map2_dbl(microbe, clinical_var, ~{
      x <- abund_sel[[.x]]
      y <- clinical_df[[.y]]
      if (length(unique(na.omit(x))) < 2 || length(unique(na.omit(y))) < 2) return(NA)
      cor.test(x, y, method = "pearson")$p.value
    }),
    p_adj = p.adjust(p_value, method = "fdr")
  )

cor_results <- cor_results %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

cor_results$clinical_var <- factor(cor_results$clinical_var, levels = desired_order)


pp <- ggplot(cor_results, aes(x = clinical_var, y = microbe, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif), size = 4, color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Pearson R"
  ) +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

pp



ggsave(plot=pp,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/correlation_micro_clinical_doctor.svg", height = 5, width = 9, dpi=300)
ggsave(plot=pp,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/correlation_micro_clinical_doctor.png", height = 5, width = 9, dpi=300)

