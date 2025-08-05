library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggpmisc)
library(dplyr)
library(microbiomeMarker)
library(ggplot2)
library(ggpmisc)
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

library(phyloseq)
library(dplyr)

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

############################################# LINEAR REGRESSION



phy_genus <- tax_glom(phy, taxrank = "Genus")



phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))

abund <- as.data.frame(otu_table(phy_genus_rel))
taxa <- as.data.frame(tax_table(phy_genus_rel))


blautia_row <- rownames(taxa)[taxa$Genus == "Escherichia-Shigella"]
blautia_abund <- as.numeric(abund[blautia_row, ])
names(blautia_abund) <- colnames(abund)
# Make sure sample names align
common_samples <- intersect(rownames(sam_new), names(blautia_abund))



plot_df <- sam_new[common_samples, ]
plot_df$Blautia <- blautia_abund[common_samples]



plot_df$Type <- ifelse(grepl("PDM", plot_df$sample_information), "PDM",
                                ifelse(grepl("K", plot_df$sample_information), "K", "DM"))

# Predict PDM vs DM 
plot_df <- plot_df %>%
  filter(Type %in% c("PDM")) %>%
  mutate(Type = factor(Type))


pdm_df <- subset(plot_df, Type == "PDM")
pdm_df$Type <- "T3cDM"



p <- ggplot(pdm_df, aes(x = Blautia, y = `Pancreatectomy`, color = Type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_poly_eq(
    aes(label = paste(..p.value.label..)),
    formula = y ~ x, parse = TRUE,
    label.x.npc = "left", label.y.npc = 0.95
  ) +
  labs(title = "", x = "Escherichia-Shigella", y = "Pancreatectomy") +
  scale_color_manual(values = c("DM" = "#E1812C", "T3cDM" = "#3A923A", "K" = "blue")) +
  theme_minimal()
p


#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/linear_regression_Escherichia-Shigella_Pancreatectomy.svg", height = 3, width = 5, dpi=300)
#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/linear_regression_Escherichia-Shigella_Pancreatectomy.png", height = 3, width = 5, dpi=300)



############################################## new 

pancreatectomy_map <- c(
  "PDM19" = "Teilresektion links",
  "PDM12" = "Teilresektion rechts",
  "PDM20" = "Teilresektion links",
  "PDM5"  = "Teilresektion links",
  "PDM16" = "Teilresektion links",
  "PDM23" = "Teilresektion rechts",
  "PDM6"  = "Resektion",
  "PDM8"  = "Teilresektion rechts",
  "PDM3"  = "Teilresektion links",
  "PDM9"  = "Resektion",
  "PDM18" = "Resektion",
  "PDM4"  = "Teilresektion links"
)

# Create a new column by mapping from rownames (or another identifier column)
sam_new$Pancreatectomy_cat <- pancreatectomy_map[sam_new$sample_information]
sam_new$Pancreatectomy_cat <- recode(
  sam_new$Pancreatectomy_cat,
  "Teilresektion links" = "Partial resection left",
  "Teilresektion rechts" = "Partial resection right",
  "Resektion" = "Resection"
)
sam_new$bacteria_ab <- blautia_abund[rownames(sam_new)]

pan <- ggplot(sam_new, aes(x = Pancreatectomy_cat, y = bacteria_ab)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgreen") +
  geom_jitter(width = 0.2, color = "darkgreen", size = 1.5) +
  labs(
    title = "",
    x = "Pancreatectomy Category",
    y = "Escherichia-Shigella relative abundance"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


my_comparisons <- list(
  c("Resection", "Partial resection left"),
  c("Resection", "Partial resection right"),
  c("Partial resection left", "Partial resection right")
)


summary_stats <- sam_new %>%
  group_by(Pancreatectomy_cat) %>%
  summarise(
    n = sum(!is.na(bacteria_ab)),
    median = median(bacteria_ab, na.rm = TRUE),
    IQR = IQR(bacteria_ab, na.rm = TRUE)
  )

print(summary_stats)

# Convert to wide named vector for easy access
medians <- summary_stats$median
names(medians) <- summary_stats$Pancreatectomy_cat

# Compute fold changes (median1 / median2)
fold_change_1 <- medians["Resection"] / medians["Partial resection left"]
fold_change_2 <- medians["Resection"] / medians["Partial resection right"]
fold_change_3 <- medians["Partial resection right"] / medians["Partial resection left"]

cat("Fold change (Resection vs Left Partial):", round(fold_change_1, 2), "\n")
cat("Fold change (Resection vs Right Partial):", round(fold_change_2, 2), "\n")
cat("Fold change (Right Partial vs Left Partial):", round(fold_change_3, 2), "\n")



pan
#ggsave(plot=pan,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/boxplot_Escherichia-Shigella_Pancreatectomy.svg", height = 5, width = 3, dpi=300)
#ggsave(plot=pan,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/boxplot_Escherichia-Shigella_Pancreatectomy.png", height = 5, width = 3, dpi=300)

