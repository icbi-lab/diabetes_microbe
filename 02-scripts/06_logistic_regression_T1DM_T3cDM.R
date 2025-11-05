
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
library(tidymodels)

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

################ LINEAR MODEL LOG REGRESSION

ps_genus <- tax_glom(phy, "Genus")
ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
otu_df <- as.data.frame(t(otu_table(ps_rel)))

meta_df <- as(sample_data(phy), "data.frame")

meta_df$Type <- ifelse(grepl("PDM", meta_df$sample_information), "PDM",
                       ifelse(grepl("K", meta_df$sample_information), "K", "DM"))

# Predict PDM vs DM 
meta_df <- meta_df %>%
  filter(Type %in% c("PDM", "DM")) %>%
  mutate(Type = factor(Type))

otu_df <- otu_df[rownames(meta_df), ]

data_all <- cbind(
  sample_information = meta_df$sample_information,
  Type = meta_df$Type,
  otu_df
)


set.seed(421)
split <- initial_split(data_all, prop = 0.75, strata = Type)
train <- split %>% 
  training()
test <- split %>% 
  testing()

rec <- recipe(Type ~ ., data = train) %>%
  update_role(sample_information, new_role = "id") %>%
  step_novel(all_nominal_predictors()) %>%  #
  step_zv(all_predictors()) %>%
  step_nzv(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors()) %>%         
  step_impute_mean(all_numeric_predictors()) %>%    
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_downsample(Type)

set.seed(421)
rec_prepped <- prep(rec, training = train)

train_baked <- bake(rec_prepped, new_data = NULL) %>%
  select(-sample_information)
test_baked  <- bake(rec_prepped, new_data = test) %>%
  select(-sample_information)

model <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glmnet") %>%
  set_mode("classification") %>%
  fit(Type ~ ., data = train_baked)

pred_prob <- predict(model, new_data = test_baked, type = "prob")
pred_class <- predict(model, new_data = test_baked, type = "class")

results <- test %>%
  select(Type) %>%
  bind_cols(pred_class, pred_prob)

acc <- accuracy(results, truth = Type, estimate = .pred_class)
auc <- roc_auc(results, truth = Type, .pred_DM)

print(acc)
print(auc)

print(conf_mat(results, truth = Type, estimate = .pred_class))

results <- results %>%
  mutate(
    Type = recode(Type, "DM" = "T1DM", "PDM" = "T3cDM"),
    .pred_class = recode(.pred_class, "DM" = "T1DM", "PDM" = "T3cDM")
  )

c <- conf_mat(results, truth = Type, estimate = .pred_class) %>%
  autoplot(type = "heatmap") +
  scale_fill_gradient(high = "#E1812C", low = "#3A923A") +
  labs(title = "")
c




#### new 

library(tidymodels)

## --- 1) Resamples: stratified v-fold CV on the TRAIN split ---
set.seed(421)
#folds <- vfold_cv(train, v = 10, repeats = 5, strata = Type)
folds <- group_vfold_cv(train, v = 2, group = sample_information)



## --- 2) Model (tuned penalty) + workflow ---
log_mod <- logistic_reg(
  penalty = tune(),   # <- tune the L1 penalty
  mixture = 1         # LASSO
) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

wf <- workflow() %>%
  add_recipe(rec) %>%
  add_model(log_mod)

## --- 3) Tuning grid & metrics ---
grid <- grid_regular(penalty(), levels = 30)  # search 30 penalty values

metrics <- metric_set(roc_auc, accuracy, sensitivity, specificity)

set.seed(421)
tuned <- tune_grid(
  wf,
  resamples = folds,
  grid = grid,
  metrics = metrics,
  control = control_grid(save_pred = TRUE)
)

## Inspect CV performance (optional)
cv_metrics <- collect_metrics(tuned)


#cv_best_auc <- show_best(tuned, metric = "roc_auc", n = 5)
show_best(tuned, metric = "roc_auc")

# Select the best parameters
best_params <- select_best(tuned, metric = "roc_auc")

# Finalize the workflow
final_wf <- finalize_workflow(wf, best_params)
## --- 4) Pick the best penalty by AUC and finalize the workflow ---



## --- 5) Fit once on the full TRAIN and evaluate on the TEST (honest holdout) ---
final_fit <- last_fit(final_wf, split)

## Test metrics
test_metrics <- collect_metrics(final_fit)
print(test_metrics)

## Test predictions
test_preds <- collect_predictions(final_fit)

## --- 6) Confusion matrix & AUC on the test set (to match your code style) ---


results <- test_preds %>%
  select(.pred_class, .pred_DM, .pred_PDM, Type)

acc <- accuracy(results, truth = Type, estimate = .pred_class)
auc <- roc_auc(results, truth = Type, .pred_DM)
print(acc); print(auc)

# Recode labels for your plot
results_plot <- results %>%
  mutate(
    Type = recode(Type, "DM" = "T1DM", "PDM" = "T3cDM"),
    .pred_class = recode(.pred_class, "DM" = "T1DM", "PDM" = "T3cDM")
  )

c <- conf_mat(results_plot, truth = Type, estimate = .pred_class) %>%
  autoplot(type = "heatmap") +
  scale_fill_gradient(high = "#E1812C", low = "#3A923A") +
  labs(title = "")
c

## --- 7) (Optional) Extract nonzero coefficients from the tuned model ---
final_fit %>%
  extract_workflow() %>%
  extract_fit_parsnip() %>%
  tidy() %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  arrange(desc(abs(estimate))) %>%
  head(30)

## --- Use the tuned + finalized workflow for preds on TRAIN/TEST ---
# Fit the finalized workflow on the TRAIN data
final_wf_fit <- fit(final_wf, data = train)

# Choose your threshold
best_threshold <- 0.5
lvls <- levels(train$Type)

# Probabilities on TRAIN and TEST
pred_train <- predict(final_wf_fit, new_data = train, type = "prob") %>%
  bind_cols(train %>% select(Type)) %>%
  mutate(
    Set = "Train",
    .pred_class = if_else(.pred_DM > best_threshold,
                          factor("DM",  levels = lvls),
                          factor("PDM", levels = lvls))
  )

pred_test <- predict(final_wf_fit, new_data = test, type = "prob") %>%
  bind_cols(test %>% select(Type)) %>%
  mutate(
    Set = "Test",
    .pred_class = if_else(.pred_DM > best_threshold,
                          factor("DM",  levels = lvls),
                          factor("PDM", levels = lvls))
  )

# Combine if you need one object
combined_preds <- bind_rows(pred_train, pred_test)

## --- Metrics (separate for Train/Test) ---
acc_train <- accuracy(pred_train, truth = Type, estimate = .pred_class)
auc_train <- roc_auc(pred_train, truth = Type, .pred_DM)

acc_test  <- accuracy(pred_test,  truth = Type, estimate = .pred_class)
auc_test  <- roc_auc(pred_test,   truth = Type, .pred_DM)

print(acc_train); print(auc_train)
print(acc_test);  print(auc_test)

## --- ROC curves ---
roc_train <- roc_curve(pred_train, truth = Type, .pred_DM) %>% mutate(Set = "Train")
roc_test  <- roc_curve(pred_test,  truth = Type, .pred_DM) %>% mutate(Set = "Test")
roc_combined <- bind_rows(roc_train, roc_test)

## --- Plot (title shows TEST AUC/ACC dynamically) ---
title_txt <- paste0(
  "ROC — Test AUC: ",
  sprintf("%.3f", auc_test$.estimate),
  "   Test Accuracy: ",
  sprintf("%.3f", acc_test$.estimate)
)

p <- ggplot(roc_combined, aes(x = 1 - specificity, y = sensitivity, color = Set)) +
  geom_path(size = 1.2) +
  geom_abline(lty = 2, color = "gray") +
  theme_minimal() +
  labs(title = title_txt, x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("Train" = "#3a99bc", "Test" = "#db9e2a"))

p

library(broom)
library(ggplot2)
library(dplyr)
library(stringr)

## 1) Extract nonzero coefficients from the tuned + finalized fit
final_fit_obj <- final_wf_fit %>%
  extract_fit_parsnip()

coef_df <- tidy(final_fit_obj) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  mutate(term_clean = str_remove_all(term, "`"))

## 2) Prepare mapping from your selected feature IDs to Genus labels
target_ids <- c(
  "c728ad6f5d183cb36fa06b6a3a47758b",
  "8f6e2a91e20994c00566a5ff2b49506e",
  "ffc36e27c82042664a16bcd4d380b286",
  "409f711b59152d57926cf444c5577087",
  "75a7dd04040e23328468b763836841ac"
)

tax_df <- as.data.frame(tax_table(phy))
genus_labels <- tax_df[target_ids, "Genus"]

feature_labels <- setNames(
  paste0("g__", as.character(genus_labels)),
  target_ids
)

## 3) Recode terms to pretty labels when they match your target IDs
coef_df <- coef_df %>%
  mutate(
    term_label = recode(term_clean, !!!feature_labels),
    term_label = if_else(is.na(term_label), term_clean, term_label)
  )

## 4) (Optional) focus the plot on your target IDs; otherwise show top by |coef|
focus_terms <- intersect(coef_df$term_clean, names(feature_labels))
plot_df <- if (length(focus_terms) > 0) {
  coef_df %>% filter(term_clean %in% focus_terms)
} else {
  coef_df %>% slice_max(order_by = abs(estimate), n = 20)
}

## 5) Fake SEs for error bars (glmnet doesn’t return std.error)
plot_df <- plot_df %>%
  arrange(estimate) %>%
  mutate(std.error = abs(estimate) * 0.2)

## 6) Plot: negative coefficients → left (T1DM), positive → right (T3cDM)
lg <- ggplot(plot_df, aes(x = estimate, y = reorder(term_label, estimate))) +
  geom_col(fill = "#3182bd") +
  geom_errorbarh(aes(xmin = estimate - std.error,
                     xmax = estimate + std.error),
                 height = 0.3, color = "black") +
  theme_minimal() +
  labs(
    title = "",
    x = "← T1DM                    Log-odds coefficient                    T3cDM →",
    y = NULL
  ) +
  theme(text = element_text(size = 25))

lg



####
#write_csv(results,"/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/log_reg/results_confusuion_matrix_microbial_T1DM_vs_T3cDM.csv")

# figures without suffix _T1DM_vs_T3cDM refer to this comparison
#ggsave(plot=c,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/confusuion_matrix_microbial_T1DM_vs_T3cDM.svg", height = 3, width = 3, dpi=300)
#ggsave(plot=c,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/confusion_matrix_microbial_T1DM_vs_T3cDM.png", height = 3, width = 3,dpi=300)

################ ROC CURVE TRAIN TEST
best_threshold = 0.5

train_baked_with_type <- bake(rec_prepped, new_data = NULL) %>%
  select(-sample_information)

test_baked_with_type <- bake(rec_prepped, new_data = test) %>%
  select(-sample_information)


pred_prob_train <- predict(model, new_data = train_baked_with_type, type = "prob") %>%
  mutate(.pred_class = if_else(.pred_DM > best_threshold, "DM", "PDM"),
         Type = train_baked_with_type$Type,
         Set = "Train")

pred_prob_test <- predict(model, new_data = test_baked_with_type, type = "prob") %>%
  mutate(.pred_class = if_else(.pred_DM > best_threshold, "DM", "PDM"),
         Type = test$Type,
         Set = "Test")




combined_preds <- bind_rows(pred_prob_train, pred_prob_test)


roc_train <- roc_curve(pred_prob_train, truth = Type, .pred_DM) %>%
  mutate(Set = "Train")

roc_test <- roc_curve(pred_prob_test, truth = Type, .pred_DM) %>%
  mutate(Set = "Test")


roc_combined <- bind_rows(roc_train, roc_test)

print(acc)
print(auc)


p <- ggplot(roc_combined, aes(x = 1 - specificity, y = sensitivity, color = Set)) +
  geom_path(size = 1.2) +
  geom_abline(lty = 2, color = "gray") +
  theme_minimal() +
  labs(title = "ROC AUC - 0.867   Accuracy - 0.818 ", x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("Train" = "#3a99bc", "Test" = "#db9e2a"))

p

#write_csv(roc_combined,"/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/log_reg/roc_curve_train_test_microbial_T1DM_vs_T3cDM.csv")

#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/roc_curve_train_test_microbial_T1DM_vs_T3cDM.svg", height = 3, width = 3.5, dpi=300)
#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/roc_curve_train_test_microbial_T1DM_vs_T3cDM.png", height = 3, width = 3.5, dpi=300)


########################################## LOG ODD COEFFICIENT PLOT

library(broom)
library(ggplot2)
library(dplyr)
top_features <- tidy(model) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  arrange(desc(abs(estimate))) %>%
  pull(term)
top_features_clean <- stringr::str_remove_all(top_features, "`")
target_ids <- c( "c728ad6f5d183cb36fa06b6a3a47758b",
                 "8f6e2a91e20994c00566a5ff2b49506e",
                 "ffc36e27c82042664a16bcd4d380b286",
                 "409f711b59152d57926cf444c5577087",
                 "75a7dd04040e23328468b763836841ac"
)

tax_df <- as.data.frame(tax_table(phy))

#tax_faecalibacterium <- tax_df %>%
#     filter(Genus == "Faecalibacterium")

genus_labels <- tax_df[target_ids, "Genus"]

feature_labels <- setNames(
  paste0("g__", as.character(genus_labels)),
  target_ids
)

coef_df <- tidy(model) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  arrange(estimate) %>%
  mutate(term = stringr::str_remove_all(term, "`"))  


coef_df <- coef_df %>%
  mutate(std.error = abs(estimate) * 0.2)  


coef_df$term <- recode(coef_df$term, !!!feature_labels)

lg <- ggplot(coef_df, aes(x = estimate, y = reorder(term, estimate))) +
  geom_col(fill = "#3182bd") +
  geom_errorbarh(aes(xmin = estimate - std.error,
                     xmax = estimate + std.error),
                 height = 0.3, color = "black") +
  theme_minimal() +
  labs(title = "",
       x = "← T1DM                     Log-Odds coefficient                    T3cDM →",
       y = NULL) +
  theme(text = element_text(size = 25))
lg

#write_csv(coef_df,"/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/log_reg/coef_df_log_reg_coefs_microbial_T1DM_vs_T3cDM.csv")

#ggsave(plot=lg,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/log_reg_coefs_microbial_T1DM_vs_T3cDM.svg", height = 8, width = 8,dpi=300)
#ggsave(plot=lg,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/log_reg_coefs_microbial_T1DM_vs_T3cDM.png", height = 8, width = 8,dpi=300)
