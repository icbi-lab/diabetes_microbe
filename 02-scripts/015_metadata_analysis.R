library(tidyverse)
library(readxl)
library(janitor)
library(GGally)
library(ggpubr)
library(dplyr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom.mixed)
library(lmerTest)


df_all <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/PDM merged 3.0_modified.csv")%>%
  clean_names() 

df_lab <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/laboratory_variables.csv")%>%
  clean_names() 

var_names_lab <- colnames(df_lab)
base_vars_lab <- unique(str_remove(var_names_lab, "[12]$"))

df_epi <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/epidemiological_variables.csv")%>%
  clean_names() 

var_names_epi <- colnames(df_epi)
base_vars_epi <- unique(str_remove(var_names_epi, "[12]$"))



df_all <- df_all %>%
  mutate(Type = case_when(
    grepl("^PDM", probennummer) ~ "PDM",
    grepl("^K", probennummer) ~ "K",
    TRUE ~ "DM"
  ))


var_names <- colnames(df_all)
base_vars <- unique(str_remove(var_names, "[12]$"))

base_vars <- setdiff(base_vars, c("probennummer", "Type","verstorben","pankreatektomie"))


var_names <- unlist(lapply(base_vars, function(v) paste0(v, c("1", "2"))))
var_names <- var_names[var_names %in% colnames(df_all)]  # keep existing

df_long <- df_all %>%
  mutate(ID = probennummer) %>%
  pivot_longer(
    cols = all_of(var_names),
    names_to = c(".value", "time"),
    names_pattern = "(.*)([12])"
  ) %>%
  mutate(time = as.integer(time))

results <- list()

for (var in base_vars) {

  if (!all(c(var, "time", "probennummer") %in% colnames(df_long))) next
  
  df_sub <- df_long %>%
    dplyr::select(probennummer, time,age, sex,Type, !!sym(var)) %>%
    dplyr::filter(!is.na(!!sym(var)), !is.na(Type))
  
  
  if (nrow(df_sub) < 10) next
  
  
  model <- NULL
  result <- NULL
  
  # Numeric -> LMM
  if (is.numeric(df_sub[[var]])) {
    model <- tryCatch(
      lmer(as.formula(paste(var, "~ time + Type + age + sex + (1 | probennummer)")), data = df_sub),
      error = function(e) NULL
    )
  }
  # Binary -> GLMM
  else if (all(df_sub[[var]] %in% c(0, 1))) {
    model <- tryCatch(
      glmer(as.formula(paste(var, "~ time + Type + age + sex + (1 | probennummer)")), data = df_sub, family = binomial),
      error = function(e) NULL
    )
  }
  
  if (!is.null(model)) {
    est <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
      dplyr::filter(term == "time") %>%
      dplyr::mutate(
        variable = var,
        p.value = 2 * (1 - pnorm(abs(estimate / std.error)))  # add this line
      ) %>%
      dplyr::select(variable, estimate, std.error, conf.low, conf.high, p.value)
    
    results[[var]] <- est
  }
}

forest_df <- bind_rows(results)
forest_df <- forest_df %>%
  dplyr::mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )



forest_df <- forest_df %>%
  dplyr::mutate(
    source = case_when(
      variable %in% base_vars_epi ~ "epi",
      variable %in% base_vars_lab ~ "lab",
      TRUE ~ "other"
    )
  )


run_models_by_group <- function(df_group, group_label, base_vars) {
  results <- list()
  
  for (var in base_vars) {
    if (!all(c(var, "time", "probennummer", "age", "sex") %in% colnames(df_group))) next
    
    df_sub <- df_group %>%
      dplyr::select(probennummer, time, age, sex, !!sym(var)) %>%
      dplyr::filter(!is.na(!!sym(var)), !is.na(age), !is.na(sex))
    
    if (nrow(df_sub) < 10) next
    
    model <- NULL
    
    if (is.numeric(df_sub[[var]])) {
      model <- tryCatch(
        lmer(as.formula(paste(var, "~ time + age + sex + (1 | probennummer)")), data = df_sub),
        error = function(e) NULL
      )
    }
    else if (all(df_sub[[var]] %in% c(0, 1), na.rm = TRUE)) {
      model <- tryCatch(
        glmer(as.formula(paste(var, "~ time + age + sex + (1 | probennummer)")), data = df_sub, family = binomial),
        error = function(e) NULL
      )
    }
    
    if (!is.null(model)) {
      est <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
        dplyr::filter(term == "time") %>%
        dplyr::mutate(
          variable = var,
          group = group_label,
          p.value = 2 * (1 - pnorm(abs(estimate / std.error)))
        ) %>%
        dplyr::select(variable, estimate, std.error, conf.low, conf.high, p.value, group)
      
      if ("glmerMod" %in% class(model)) {
        est <- est %>%
          dplyr::mutate(
            estimate = exp(estimate),
            conf.low = exp(conf.low),
            conf.high = exp(conf.high)
          )
      }
      
      results[[var]] <- est
    }
  }
  
  bind_rows(results)
}

df_pdm <- df_long %>% filter(Type == "PDM")
df_dm  <- df_long %>% filter(Type == "DM")

all_vars <- colnames(df_long)
base_vars <- setdiff(all_vars, c("probennummer", "time", "Type", "sex", "age"))

forest_pdm <- run_models_by_group(df_pdm, "PDM", base_vars)
forest_dm  <- run_models_by_group(df_dm,  "DM",  base_vars)

forest_df <- bind_rows(forest_pdm, forest_dm)

forest_df <- forest_df %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )
name_map <- c(
  "age" = "Age",
  "hyperlipid_mie_2018" = "Hyperlipidemia (2018)",
  "arterielle_hyperotnie" = "Arterial Hypertension",
  "bmi" = "Body Mass Index (BMI)",
  "gr_e" = "Height",
  "gewicht" = "Weight",
  "lipidsenker" = "Lipid-lowering Medication",
  "rr_medikation" = "Blood Pressure Medication",
  "langzeit_insulin" = "Long-acting Insulin",
  "kurzzeit_insulin" = "Short-acting Insulin",
  "misch_inuslin" = "Mixed Insulin",
  "masld" = "MASLD",
  "khk" = "Coronary Heart Disease (CHD)",
  "ca" = "Cancer (unspecified)",
  "leukozyten" = "Leukocytes",
  "h_moglobin" = "Hemoglobin",
  "h_matokrit" = "Hematocrit",
  "thrombozyten" = "Platelets",
  "harnstoff" = "Urea",
  "creatinin_enzym_idms_" = "Creatinine (IDMS)",
  "glomerul_re_filtrationsrate" = "Glomerular Filtration Rate (GFR)",
  "bilirubin_gesamt" = "Total Bilirubin",
  "natrium" = "Sodium",
  "got_asat_" = "AST",
  "gpt_alat_" = "ALT",
  "gamma_gt" = "GGT",
  "alkalische_phosphatase" = "Alkaline Phosphatase",
  "lactat_dehydrogenase_ldh_" = "LDH",
  "c_reaktives_prot_crp_" = "CRP",
  "quicktest_pt_" = "Prothrombin Time",
  "inr_pt_" = "INR",
  "part_thrombopl_zeit_a_ptt_" = "aPTT",
  "fibrinogen_funkt_n_clauss" = "Fibrinogen (Clauss)",
  "albumin" = "Albumin",
  "glukose" = "Glucose",
  "hb_a1c_dcct_ngsp_" = "HbA1c (DCCT/NGSP)",
  "hb_a1c_ifcc_" = "HbA1c (IFCC)",
  "cholesterin" = "Total Cholesterol",
  "non_hdl_cholesterin" = "Non-HDL Cholesterol",
  "triglyceride" = "Triglycerides",
  "hdl_cholesterin" = "HDL Cholesterol",
  "ldl_cholesterin" = "LDL Cholesterol",
  "eisen" = "Iron",
  "ferritin" = "Ferritin",
  "transferrin" = "Transferrin",
  "transferrins_ttigung" = "Transferrin Saturation",
  "troponin_t_hoch_sens_" = "Troponin T (High Sensitivity)",
  "nt_pro_bnp" = "NT-proBNP"
)

forest_df <- forest_df %>%
  mutate(variable = recode(variable, !!!name_map))


forest_df <- forest_df %>%
  group_by(group) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    sig = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )


forest_df_sig <- forest_df %>%
  filter(p.adj < 0.05) %>%
  mutate(variable = recode(variable, !!!name_map))



p <- ggplot(forest_df_sig, aes(x = estimate, y = reorder(variable, estimate), color = group)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = sig), hjust = -0.5, size = 4, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c("DM" = "#E1812C", "PDM" = "#3A923A")
  ) +
  labs(
    x = "Mean effect (Follow up vs Baseline)",
    y = "",
    title = "",
    color = ""
  ) +
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size = 12),  
    axis.text.x = element_text(size = 10), 
    legend.text = element_text(size = 11)  
  )

p <- p + theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA)
)

p <- p + 
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -1.5,
           label = "* p.adj < 0.05   ** p.adj < 0.01   *** p.adj < 0.001",
           size = 6, color = "black", fontface = "italic")



p

#ggsave(plot = p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/forest_plot_numerical_groups_padjsig.svg", height = 5, width = 7)
#ggsave(plot = p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/forest_plot_numerical_groups_padjsig.png", height = 5, width = 7)

write.csv(forest_df_sig, file = "/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/forest_df_sig_forest_plot_numerical_groups_padjsig.csv", row.names = FALSE)
