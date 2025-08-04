
library(readr)
library(dplyr)
library(purrr)
library(rstatix)
data <- read_tsv("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/materials/growthresults_exchanges.tsv")
data <- data[, -1]


data <- data %>%
  mutate(condition = recode(condition,
                            "Kontrolle" = "H",
                            "Diabetes mellitus Typ1" = "T1DM",
                            "pankreopriver Diabetes" = "T3cDM"))

comparisons <- list(c("H", "T1DM"), c("H", "T3cDM"), c("T1DM", "T3cDM"))

valid_descriptions <- data %>%
  filter(!is.na(flux)) %>%
  group_by(description, condition) %>%
  summarise(n = n(), unique_vals = n_distinct(flux), .groups = "drop") %>%
  filter(n >= 2, unique_vals > 1) %>%
  group_by(description) %>%
  summarise(n_conditions = n_distinct(condition), .groups = "drop") %>%
  filter(n_conditions == 3) %>%
  pull(description)

stat_df <- data %>%
  filter(description %in% valid_descriptions) %>%
  group_by(description) %>%
  pairwise_wilcox_test(flux ~ condition, p.adjust.method = "fdr") %>%
  filter((group1 == "H" & group2 == "T1DM") |
           (group1 == "H" & group2 == "T3cDM") |
           (group1 == "T1DM" & group2 == "T3cDM")) %>%
  ungroup()

eff_df <- data %>%
  filter(description %in% valid_descriptions) %>%
  group_by(description) %>%
  wilcox_effsize(flux ~ condition, paired = FALSE, ci = TRUE) %>%
  filter((group1 == "H" & group2 == "T1DM") |
           (group1 == "H" & group2 == "T3cDM") |
           (group1 == "T1DM" & group2 == "T3cDM")) %>%
  ungroup()

stat_df <- stat_df %>%
  left_join(eff_df %>% select(description, group1, group2, effsize), 
            by = c("description", "group1", "group2"))


stat_df <- stat_df %>%
  mutate(result_strength = case_when(
    p.adj < 0.05 & abs(effsize) >= 0.8 ~ "strong",
    p.adj < 0.05 & abs(effsize) >= 0.5 ~ "moderate",
    p.adj < 0.05 & abs(effsize) >= 0.2 ~ "weak",
    TRUE ~ "not significant"
  ))


stat_df_filtered <- stat_df %>%
  filter(result_strength %in% c("moderate", "strong","weak"))


##############################################################################
my_palette <- c("T3cDM" = "#6ABC6A", "T1DM" = "#FFA555", "H" = "#619FCA")

metabolites <- unique(stat_df_filtered$description)



data <- data %>%
  filter(description %in% metabolites) %>%
  mutate(description = factor(description, levels = metabolites))



valid_descriptions <- data %>%
  filter(!is.na(flux)) %>%
  group_by(description, condition) %>%
  summarise(n = n(), sd = sd(flux), .groups = "drop") %>%
  filter(n >= 2, sd > 0) %>%
  group_by(description) %>%
  summarise(n_conditions = n_distinct(condition), .groups = "drop") %>%
  filter(n_conditions == 3) %>%
  pull(description)


stat_df <- data %>%
  filter(description %in% valid_descriptions) %>%
  group_by(description) %>%
  pairwise_wilcox_test(flux ~ condition, p.adjust.method = "fdr") %>%
  filter((group1 == "H" & group2 == "T1DM") |
           (group1 == "H" & group2 == "T3cDM") |
           (group1 == "T1DM" & group2 == "T3cDM")) %>%
  mutate(
    comparison = paste(group1, "vs", group2),
    p.signif = symnum(
      p.adj,
      corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns"))
  ) %>%
  ungroup()

max_flux <- data %>%
  group_by(description) %>%
  summarise(y.max = max(flux, na.rm = TRUE), .groups = "drop")



stat_df <- stat_df %>%
  left_join(max_flux, by = "description") %>%
  group_by(description) %>%
  mutate(
    y.position = y.max + row_number() * (y.max * 0.15),  # stack brackets
    x.start = case_when(
      group1 == "H" & group2 == "T1DM" ~ 1,
      group1 == "H" & group2 == "T3cDM" ~ 1,
      group1 == "T1DM" & group2 == "T3cDM" ~ 2
    ),
    x.end = case_when(
      group1 == "H" & group2 == "T1DM" ~ 2,
      group1 == "H" & group2 == "T3cDM" ~ 3,
      group1 == "T1DM" & group2 == "T3cDM" ~ 3
    )
  ) %>%
  ungroup()

stat_df <- stat_df %>%
  left_join(eff_df %>% select(description, group1, group2, effsize), 
            by = c("description", "group1", "group2"))

stat_df <- stat_df %>%
  mutate(result_strength = case_when(
    p.adj < 0.05 & abs(effsize) >= 0.5 ~ "strong",
    p.adj < 0.05 & abs(effsize) <= 0.5  & abs(effsize) >= 0.3 ~ "moderate",
    p.adj < 0.05 & abs(effsize) >= 0.3 ~ "weak",
    TRUE ~ "not significant"
  ))

stat_df_filtered <- stat_df %>%
  filter(result_strength %in% c("moderate", "strong","weak"))

metabolites_filtered <- unique(stat_df_filtered$description)

data <- data %>%
  filter(description %in% metabolites_filtered) %>%
  mutate(description = factor(description, levels = metabolites_filtered))


p <- ggplot(data, aes(x = condition, y = flux, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_jitter(aes(color = condition), width = 0.2, size = 1.2, alpha = 0.6) +
  facet_wrap(~ description, scales = "free_y", ncol = 6) +
  
  
  geom_segment(
    data = stat_df,
    aes(x = x.start, xend = x.end, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  
  
  geom_text(
    data = stat_df,
    aes(x = (x.start + x.end) / 2, y = y.position + 0.03 * y.position, label = p.signif),
    inherit.aes = FALSE,
    size = 5  
  ) +
  
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = my_palette) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "",
    x = "Condition",
    y = "Flux"
  )

p

#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/metabolite_flux_all.svg", plot = p,
#       width = 20, height = 20, units = "in", dpi = 300)


#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/metabolite_flux_all.png", plot = p,
#       width = 20, height = 20, units = "in", dpi = 300)



####################### Volcano plot 


stat_df <- stat_df %>%
  mutate(comparison = paste(group1, "vs", group2))
library(ggrepel)

stat_df_unique <- stat_df[!duplicated(stat_df$description), ]

# Create volcano plot
q <- ggplot(stat_df_unique, aes(x = effsize, y = -log10(p.adj))) +
  geom_point(aes(color = result_strength), alpha = 0.8, size = 3) +
  geom_vline(xintercept = c(0, 0.5), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_text_repel(
    data = filter(stat_df,abs(effsize) >= 0.3),
    aes(label = description),
    size = 3.5, max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "strong" = "#B2182B",
    "moderate" = "#EF8A62",
    "weak" = "#FDB863",
    "not significant" = "grey70"
  )) +
  theme_minimal(base_size = 13) +
  labs(title = "",
       x = "Effect size",
       y = "-log10(P.adj)",
       color = "Result strength")

#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/volcano_plot.svg", plot = q,
#       width = 10, height = 6, units = "in", dpi = 300)
#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/volcano_plot.png", plot = q,
#       width = 10, height = 6, units = "in", dpi = 300)


###################### Boxplot with bioligical sig. effect size >0.3 and p.adj >0.05
#Selected bioligically meaningful metabolites 

filtered_metabolites <- c("Pyridoxine", "Fe3+", "Ubiquinone-8", "Zinc")
filtered_data <- data %>% filter(description %in% filtered_metabolites)
filtered_stat_df <- stat_df %>% filter(description %in% filtered_metabolites)


p <- ggplot(filtered_data, aes(x = condition, y = flux, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = condition), 
              width = 0.2, size = 1.2, alpha = 0.6) +  # Jitter layer
  facet_wrap(~ description, scales = "free_y", ncol = 6) +
  
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = my_palette) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.position = "none"
  ) +
  labs(
    title = "",
    x = "Effect size > 0.3 P.adj < 0.05",
    y = "Flux"
  )

p

#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/metabolite_flux_biological_sig.svg", plot = p,
#       width = 10, height = 6, units = "in", dpi = 300)
#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/metabolite_flux_biological_sig.png", plot = p,
#       width = 10, height = 6, units = "in", dpi = 300)


# Stats to cs v 
filtered_metabolites <- c("Pyridoxine", "Fe3+", "Ubiquinone-8", "Zinc")

# Filter the dataframe
filtered_stat_df <- stat_df %>%
  filter(description %in% filtered_metabolites)

#write_csv(filtered_stat_df, "/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/supplementary_tables/microbial_communitystat_df_metabolite_flux_biological_sig.csv")

