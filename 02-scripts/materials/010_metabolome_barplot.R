library(readr)
library(rstatix)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(broom)
library(httr)
library(jsonlite)
library(purrr)
library(ggplot2)
library(forcats)
library(viridis)   


df <- read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/tables_01/materials/exchange_fluxes_all_taxa.csv")
#df <-  read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/growthresults_exchanges.tsv")
df <- df[, -1]

ef <- read_csv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/tables_01/materials/exchange_fluxes.csv")

levels_keep <- c("T1DM", "T3cDM")
dat <- df %>%
  filter(type %in% levels_keep) 


tvals <- dat %>%
  group_by(metabolite) %>%
  filter(n_distinct(type) == 2) %>%
  group_split() %>%
  map_dfr(function(d) {
    fit <- try(lm(abundance ~ type, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(tibble())       
    tidy(fit) %>%
      filter(term == "conditionpankreopriver Diabetes") %>%
      transmute(
        metabolite = d$metabolite[1],
        t_value = statistic,
        estimate = estimate,
        p_value = p.value
      )
  }) %>%
  arrange(t_value)


get_bigg_metabolite <- function(id) {
  url <- paste0("http://bigg.ucsd.edu/api/v2/universal/metabolites/", id)
  res <- try(GET(url), silent = TRUE)
  if (inherits(res, "try-error") || res$status_code != 200) return(NULL)
  cont <- content(res, as = "text", encoding = "UTF-8")
  parsed <- fromJSON(cont, flatten = TRUE)
  tibble(
    metabolite = id,
    description = parsed$name,
    formula = parsed$formula,
    kegg_id = parsed$database_links$KEGG,
    category = parsed$major_bio_category  
  )
}


ids <- gsub("\\[e\\]$", "", tvals$metabolite)


anno <- map_dfr(ids, get_bigg_metabolite)

data_keyed <- dat %>%
  mutate(bigg_id = str_remove(metabolite, "\\[[^\\]]+\\]$"))


tvals_annot <- tvals %>%
  mutate(metabolite_clean = gsub("\\[e\\]$", "", metabolite)) %>%
  left_join(anno, by = c("metabolite_clean" = "metabolite"))


tvals_annot <- tvals_annot %>% drop_na(description)

tvals_plot <- tvals_annot %>%
  filter(!is.na(metabolite), !is.na(t_value)) %>%     
  mutate(metabolite = fct_reorder(metabolite, t_value, 
                                  .fun = median, .na_rm = TRUE))

top30 <- tvals_plot %>%
  arrange(desc(abs(t_value))) %>%  slice_head(n = 23)


top30_unique <- top30 %>%
  distinct(description, .keep_all = TRUE)


lookup <- c(
  # existing mappings
  "Oxalate" = "Carboxylic acids",
  "Octadecenoate (n-C18:1)" = "Lipids",
  "Acetoacetate" = "Carboxylic acids",
  "Hydrogen" = "Other",
  "Pullulan (n=1200 repeat units, alpha-1,4 and alph-1,6 bounds)" = "Carbohydrates",
  "12-Dehydrocholic acid; 12-Oxodeoxycholic acid; 12oxo-3alpha,7alpha-Dihydroxy-5beta-cholan-24-oic acid" = "Bile acids",
  "Cytosine" = "Nucleotides",
  "D-Mannose" = "Carbohydrates",
  "N-Acetylneuraminate" = "Carbohydrates",
  "D-Fructose" = "Carbohydrates",
  "Fumarate" = "Carboxylic acids",
  "Reduced glutathione" = "Amino acid related",
  "Deoxyguanosine" = "Nucleotides",
  "NMN C11H14N2O8P" = "Nucleotides",
  "Putrescine" = "Amino acid related",
  "Glycolate C2H3O3" = "Carboxylic acids",
  "D-Galactose" = "Carbohydrates",
  "Pectins" = "Carbohydrates",
  "Starch" = "Carbohydrates",
  "1,2-Diacyl-sn-glycerol (dioctadecanoyl, n-C18:0)" = "Lipids",
  "Pullulan 1200" = "Carbohydrates",
  "12-Dehydrocholic acid" = "Bile acids",
  "2-Oxoglutarate" = "Carboxylic acids",
  "Maltose C12H22O11" = "Carbohydrates",
  "1,2-Diacyl-sn-glycerol" = "Lipids",
  "Methanol" = "Other",
  "Nicotinamide" = "Nucleotides",     
  "Xanthine" = "Nucleotides"     ,  
  "Propionate (n-C3:0)" = "SCFA",
  "Acetate" = "SCFA",
  "Butyrate (n-C4:0)" = "SCFA",
  "Formate" = "SCFA"
)

top30_unique$class <- lookup[top30_unique$description]

top30_unique <- top30_unique %>%
  mutate(description = case_when(
    description == "12-Dehydrocholic acid; 12-Oxodeoxycholic acid; 12oxo-3alpha,7alpha-Dihydroxy-5beta-cholan-24-oic acid" ~ "12-Dehydrocholic acid",
    description == "Starch, structure 1 (1,6-{7[1,4-Glc], 4[1,4-Glc]})" ~ "Starch",
    description == "1,2-Diacyl-sn-glycerol (dioctadecanoyl, n-C18:0)" ~ "1,2-Diacyl-sn-glycerol",  description == "Pullulan (n=1200 repeat units, alpha-1,4 and alph-1,6 bounds)" ~ "Pullulan 1200",
    TRUE ~ description
  ))

extra_rows <- tvals_annot %>%
  filter(description %in% c("Propionate (n-C3:0)", "Acetate","Butyrate (n-C4:0)","Formate"))

top30_unique <- bind_rows(top30_unique, extra_rows)

top30_unique <- top30_unique %>%
  arrange(t_value) %>%
  mutate(metabolite = factor(metabolite, levels = unique(metabolite)))
library(RColorBrewer)

top30_unique <- top30_unique %>%
  filter(abs(t_value) <= 2.5)

p <- ggplot(top30_unique, aes(x = t_value, 
                         y = fct_reorder(description, t_value), 
                         fill = class)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed")+
 # scale_fill_brewer(palette = "Set1")  + 
  scale_fill_viridis_d(option = "viridis") +  # options: "viridis", "plasma", "magma", "cividis"
  labs(
    x = "t value (T3cDM vs T1DM)",
    y = "Metabolite",
    fill = "Class"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 18, margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5)
  )

p

#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/metabolite_abundance_barplot_linearreg.svg", plot = p,
#       width = 15, height = 15, units = "in", dpi = 300)

#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/metabolite_abundance_barplot_linearreg.png", plot = p,
#      width = 15, height =15 , units = "in", dpi = 300)



