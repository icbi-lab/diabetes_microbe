
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
path <- "/data/scratch/kvalem/projects/2024/diabetes_microbe/02-scripts/materials/escherichia"


files <- c("growth_rates_healthy.csv", 
           "growth_rates_t1dm.csv", 
           "growth_rates_t3cdm.csv")

data <- files %>%
  map_df(~ read_csv(file.path(path, .x)) %>%
           mutate(group = case_when(
             .x == "growth_rates_healthy.csv" ~ "H",
             .x == "growth_rates_t1dm.csv"    ~ "T1DM",
             .x == "growth_rates_t3cdm.csv"   ~ "T3cDM"
           )))


data <- data %>% 
  rename(condition = group)

exchanges <-  read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/growthresults_exchanges.tsv")
exchanges <- exchanges[, -1]


filtered_data <- data %>%
  mutate(condition = recode(condition,
                            "Kontrolle" = "H",
                            "Diabetes mellitus Typ1" = "T1DM",
                            "pankreopriver Diabetes" = "T3cDM"))

taxa_keep <- c("Akkermansia",
  "Blautia",
  "Escherichia",
"Streptococcus",
"Faecalibacterium",
"Subdoligranulum")  # Fusobacterium omitted

filtered_data <- filtered_data %>%
  filter(taxon %in% taxa_keep)

my_palette <- c("T3cDM" = "#6ABC6A", "T1DM" = "#FFA555", "H" = "#619FCA")


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
  filter(
    magnitude %in% c("moderate", "large"),
    p.adj < 0.05
  ) %>%
  pull(taxon) %>%
  unique()

q <- ggplot(filtered_data, aes(x = taxon, y = growth_rate, fill = condition)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = my_palette) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    position = position_dodge(width = 0.8),
    label.y.npc = "top"
  ) +  scale_y_continuous(limits = c(0, 0.1)) +
  labs(x = "* p.adj < 0.05", y = "Growth Rate [1/h]", fill = "Condition") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),axis.text.y = element_text( size = 16),
    panel.grid.major.x = element_blank()
  )

q


#write_csv(wilcox_results, "/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/microbial_community/wilcox_results_growthrate_boxplot.csv")
#write_csv(effect_sizes, "/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/microbial_community/effect_sizes_growthrate_boxplot.csv")

#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/growthrate_boxplot.svg", plot = q,
#       width = 10, height = 6, units = "in", dpi = 300)
#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/growthrate_boxplot.png", plot = q,
#       width = 10, height = 6, units = "in", dpi = 300)
