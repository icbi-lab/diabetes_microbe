
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
  filter(Type %in% c("PDM", "K")) %>%
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

model <- logistic_reg(penalty = 0.05, mixture = 1) %>%
  set_engine("glmnet") %>%
  set_mode("classification") %>%
  fit(Type ~ ., data = train_baked)

pred_prob <- predict(model, new_data = test_baked, type = "prob")
pred_class <- predict(model, new_data = test_baked, type = "class")

results <- test %>%
  select(Type) %>%
  bind_cols(pred_class, pred_prob)
acc <- accuracy(results, truth = Type, estimate = .pred_class)
auc <- roc_auc(results, truth = Type, .pred_K)

print(acc)
print(auc)

print(conf_mat(results, truth = Type, estimate = .pred_class))

results <- results %>%
  mutate(
    Type = recode(Type, "PDM" = "T3cDM", "K" = "H"),
    .pred_class = recode(.pred_class, "PDM" = "T3cDM", "K" = "H")
  )

c <- conf_mat(results, truth = Type, estimate = .pred_class) %>%
  autoplot(type = "heatmap") +
  scale_fill_gradient(high = "#E1812C", low = "#3A923A") +
  labs(title = "")
c

write_csv(results,"/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/log_reg/results_confusuion_matrix_microbial_T3cDM_vs_H.csv")

#ggsave(plot=c,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/confusuion_matrix_microbial_T3cDM_vs_H.svg", height = 3, width = 3, dpi=300)
#ggsave(plot=c,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/confusion_matrix_microbial_T3cDM_vs_H.png", height = 3, width = 3,dpi=300)

################ ROC CURVE TRAIN TEST
best_threshold = 0.1

train_baked_with_type <- bake(rec_prepped, new_data = NULL) %>%
  select(-sample_information)

test_baked_with_type <- bake(rec_prepped, new_data = test) %>%
  select(-sample_information)


pred_prob_train <- predict(model, new_data = train_baked_with_type, type = "prob") %>%
  mutate(.pred_class = if_else(.pred_K > best_threshold, "K", "PDM"),
         Type = train_baked_with_type$Type,
         Set = "Train")

pred_prob_test <- predict(model, new_data = test_baked_with_type, type = "prob") %>%
  mutate(.pred_class = if_else(.pred_K > best_threshold, "K", "PDM"),
         Type = test$Type,
         Set = "Test")




combined_preds <- bind_rows(pred_prob_train, pred_prob_test)


roc_train <- roc_curve(pred_prob_train, truth = Type, .pred_K) %>%
  mutate(Set = "Train")

roc_test <- roc_curve(pred_prob_test, truth = Type, .pred_K) %>%
  mutate(Set = "Test")


roc_combined <- bind_rows(roc_train, roc_test)

print(acc)
print(auc)


p <- ggplot(roc_combined, aes(x = 1 - specificity, y = sensitivity, color = Set)) +
  geom_path(size = 1.2) +
  geom_abline(lty = 2, color = "gray") +
  theme_minimal() +
  labs(title = "ROC AUC -  0.867 Accuracy - 0.75 ", x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("Train" = "#3a99bc", "Test" = "#db9e2a"))

p

#write_csv(roc_combined,"/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/log_reg/roc_curve_train_test_microbial_T3cDM_vs_H.csv")



#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/roc_curve_train_test_microbial_T3cDM_vs_H.svg", height = 3, width = 3.5, dpi=300)
#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/roc_curve_train_test_microbial_T3cDM_vs_H.png", height = 3, width = 3.5, dpi=300)


########################################## LOG ODD COEFFICIENT PLOT

library(broom)
library(ggplot2)
library(dplyr)
top_features <- tidy(model) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  arrange(desc(abs(estimate))) %>%
  pull(term)


top_features <- top_features[1:5]

target_ids <- c( "c6fdf631c16099e62a83c66db8a7a6f1",
                 "8f6e2a91e20994c00566a5ff2b49506e",
                 "409f711b59152d57926cf444c5577087",
                 "dfa833b266bd2993b86feab3617b34c3",
                 "75a7dd04040e23328468b763836841ac"
)

tax_df <- as.data.frame(tax_table(phy))

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
  arrange(desc(abs(estimate))) %>%
  slice(1:5)


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
       x = "← T3cDM                     Log-Odds coefficient                    H →",
       y = NULL) +
  theme(text = element_text(size = 18))
lg

#write_csv(coef_df,"/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/log_reg/coef_df_log_reg_coefs_microbial_T3cDM_vs_H.csv")



#ggsave(plot=lg,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/log_reg_coefs_microbial_T3cDM_vs_H.svg", height = 8, width = 8,dpi=300)
#ggsave(plot=lg,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/log_reg_coefs_microbial_T3cDM_vs_H.png", height = 8, width = 8,dpi=300)
