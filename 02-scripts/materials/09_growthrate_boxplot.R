
library(readr)
library(ggplot2)
library(rstatix)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggpubr)


data <- read_csv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/growth_rates.csv")
exchanges <-  read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/growthresults_exchanges.tsv")
exchanges <- exchanges[, -1]

data <- data %>%
  left_join(exchanges %>% select(sample_id, condition), by = "sample_id")

filtered_data <- data %>%
  mutate(condition = recode(condition,
                            "Kontrolle" = "H",
                            "Diabetes mellitus Typ1" = "T1DM",
                            "pankreopriver Diabetes" = "T3cDM"))



my_palette <- c("T3cDM" = "#6ABC6A", "T1DM" = "#FFA555", "H" = "#619FCA")

# Plotting all taxa 
ggplot(filtered_data, aes(x = taxon, y = growth_rate, fill = condition)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = my_palette) +
  theme_minimal() +
  labs(x = "Taxon", y = "Growth Rate", fill = "Condition", color = "Condition") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


taxa_with_all_conditions <- filtered_data %>%
  group_by(taxon) %>%
  summarise(n_conditions = n_distinct(condition)) %>%
  filter(n_conditions == 3) %>%
  pull(taxon)

clean_data <- filtered_data %>%
  filter(taxon %in% taxa_with_all_conditions)

clean_data <- clean_data %>%
  mutate(taxon = factor(taxon, levels = unique(taxon)))

clean_data %>%
  group_by(taxon, condition) %>%
  summarise(n = n()) %>%
  arrange(n)

comparisons <- list(c("H", "T1DM"), c("H", "T3cDM"), c("T1DM", "T3cDM"))

valid_taxa <- clean_data %>%
  group_by(taxon, condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 2) %>%
  group_by(taxon) %>%
  summarise(n_conditions = n_distinct(condition), .groups = "drop") %>%
  filter(n_conditions == 3) %>%
  pull(taxon)

filtered_for_wilcox <- clean_data %>%
  filter(taxon %in% valid_taxa)

filtered_ok <- filtered_for_wilcox %>%
  filter(!taxon %in% c("Anaerostipes", "Coprococcus"))

wilcox_results <- filtered_ok %>%
  group_by(taxon) %>%
  pairwise_wilcox_test(
    growth_rate ~ condition,
    comparisons = list(c("H", "T1DM"), c("H", "T3cDM"), c("T1DM", "T3cDM")),
    p.adjust.method = "BH"
  ) %>%
  ungroup()

effect_sizes <- filtered_ok %>%
  group_by(taxon) %>%
  wilcox_effsize(growth_rate ~ condition, comparisons = list(
    c("H", "T1DM"), c("H", "T3cDM"), c("T1DM", "T3cDM")
  )) %>%
  ungroup()

final_results <- wilcox_results %>%
  left_join(effect_sizes, by = c("taxon", "group1", "group2"))


taxa_large_effect <- final_results %>%
  filter(magnitude %in% c("moderate", "large")) %>%
  pull(taxon) %>%
  unique()

# Plotting only taxa with moderate and large effect size 
q <- ggplot(clean_data %>% filter(taxon %in% taxa_large_effect), aes(x = taxon, y = growth_rate, fill = condition)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = my_palette) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    position = position_dodge(width = 0.8),
    label.y.npc = "top"
  ) +
  labs(x = "Effect size > 0.3 P.adj < 0.05", y = "Growth Rate [1/h]", fill = "Condition") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),axis.text.y = element_text( size = 16),
    panel.grid.major.x = element_blank()
  )

q
ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/growthrate_boxplot.svg", plot = q,
       width = 10, height = 6, units = "in", dpi = 300)
ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/growthrate_boxplot.png", plot = q,
       width = 10, height = 6, units = "in", dpi = 300)
