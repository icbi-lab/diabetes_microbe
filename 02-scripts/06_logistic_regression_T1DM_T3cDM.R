
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
library(tidymodels)
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

#ggsave(plot=c,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/confusuion_matrix_microbial.svg", height = 3, width = 3, dpi=300)
#ggsave(plot=c,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/confusion_matrix_microbial.png", height = 3, width = 3,dpi=300)

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


#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/roc_curve_train_test_microbial.svg", height = 3, width = 3.5, dpi=300)
#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/roc_curve_train_test_microbial.png", height = 3, width = 3.5, dpi=300)


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
  theme(text = element_text(size = 12))
lg

#ggsave(plot=lg,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/log_reg_coefs_microbial.svg", height = 8, width = 8,dpi=300)
#ggsave(plot=lg,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/log_reg_coefs_microbial.png", height = 8, width = 8,dpi=300)
###




train_baked_labeled <- bake(rec_prepped, new_data = NULL, composition = "tibble") %>%
  select(-sample_information)

# Then select and plot
plot_data <- train_baked_labeled %>%
  select(Type, all_of(top_features_clean)) %>%
  pivot_longer(-Type, names_to = "Feature", values_to = "Value")



# Retrieve tax_table as a data frame
tax_df <- as.data.frame(tax_table(phy))

# Look up the Genus for the target taxa
genus_labels <- tax_df[target_ids, "Genus", drop = FALSE]
genus_labels

# Create mapping of hash IDs to genus
label_map <- setNames(genus_labels$Genus, rownames(genus_labels))

# Relabel Feature column
plot_data_ordered <- plot_data%>%
  mutate(Feature = recode(Feature, !!!label_map))


plot_data_ordered <- plot_data %>%
  mutate(
    Feature = factor(
      Feature,
      levels = c(
        setdiff(unique(Feature), c("c728ad6f5d183cb36fa06b6a3a47758b","8f6e2a91e20994c00566a5ff2b49506e","ffc36e27c82042664a16bcd4d380b286","c728ad6f5d183cb36fa06b6a3a47758b","75a7dd04040e23328468b763836841ac")),"c728ad6f5d183cb36fa06b6a3a47758b","8f6e2a91e20994c00566a5ff2b49506e","ffc36e27c82042664a16bcd4d380b286","c728ad6f5d183cb36fa06b6a3a47758b","75a7dd04040e23328468b763836841ac" )
    )
  )

library(ggpubr)
###
padj_df <- plot_data_ordered %>%
  group_by(Feature) %>%
  summarise(
    p = wilcox.test(Value[Type == "DM"], Value[Type == "PDM"])$p.value,
    .groups = "drop"
  ) %>%
  mutate(padj = p.adjust(p, method = "fdr"))
# Estimate y-position for label above max value per feature
padj_df <- padj_df %>%
  left_join(plot_data_ordered %>% group_by(Feature) %>% summarise(y_pos = max(Value) + 0.5), by = "Feature") %>%
  mutate(label = paste0("padj = ", signif(padj, 3)))
library(ggplot2)

p <- ggplot(plot_data_ordered, aes(x = Type, y = Value, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  facet_wrap(~ Feature, scales = "free", ncol = 5) +
  scale_fill_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A")) +
  geom_text(
    data = padj_df,
    aes(x = 1.5, y = y_pos, label = label),  # x = 1.5 centers between DM and PDM
    inherit.aes = FALSE,
    size = 3.5
  ) +
  theme_minimal() +
  labs(title = "Top Predictive Features", 
       x = "", y = "Z-score") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

p


p <- ggplot(plot_data_ordered, aes(x = Type, y = Value, fill = Type)) +
  geom_boxplot(trim = FALSE, alpha = 0.7) +
  geom_jitter(
    aes(color = Type),
    width = 0.15,
    size = 1,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  facet_wrap(~ Feature, scales = "free", ncol = 5) +
  scale_fill_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A")) +
  scale_color_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A")) +
  geom_text(
    data = padj_df,
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  theme_minimal() +
  labs(
    title = "Top Predictive Features", 
    x = "", y = "Z-score"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )


p
library(ggplot2)

library(tidyr)
library(dplyr)
library(pheatmap)
group_means <- plot_data_ordered %>%
  group_by(Feature, Type) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Type, values_from = mean_value) %>%
  column_to_rownames("Feature")

pheatmap(
  mat = as.matrix(group_means),
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Mean Z-scores by Group"
)


#### sample heatmap 
train_baked_labeled <- bake(rec_prepped, new_data = NULL, composition = "tibble")

top_features_clean <- stringr::str_remove_all(top_features, "`")

plot_data <- train_baked_labeled %>%
  select(sample_information, Type, all_of(top_features_clean)) %>%
  pivot_longer(cols = -c(sample_information, Type), names_to = "Feature", values_to = "Value")

heatmap_df <- plot_data %>%
  pivot_wider(names_from = Feature, values_from = Value) %>%
  column_to_rownames("sample_information") %>%
  mutate(across(everything(), as.numeric))

heatmap_df <- plot_data %>%
  select(sample_information, Feature, Value) %>%   # Don't select 'Type'
  pivot_wider(names_from = Feature, values_from = Value) %>%
  column_to_rownames("sample_information") %>%
  mutate(across(everything(), as.numeric))


heatmap_df <- heatmap_df %>%
  select(-Type)

sample_ann <- plot_data %>%
  distinct(sample_information, Type) %>%
  column_to_rownames("sample_information")

# Extract genera for current heatmap features
feature_ids <- colnames(heatmap_df)

# Get matching genus names from tax_df
genus_labels <- tax_df[feature_ids, "Genus"]

# Create named vector for renaming
feature_labels <- setNames(
  paste0("g__", as.character(genus_labels)),
  feature_ids
)
colnames(heatmap_df) <- feature_labels[colnames(heatmap_df)]


pheatmap(
  mat = as.matrix(heatmap_df),
  annotation_row = sample_ann,
  cluster_rows = FALSE,  # disables clustering of samples
  cluster_cols = TRUE,
  scale = "column",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  main = "Top Predictive Genera - Heatmap"
)
