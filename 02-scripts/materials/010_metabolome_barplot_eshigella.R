library(readr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(forcats)
library(ggplot2)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)

data<- read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/exchange_fluxes_escherichia.csv")

data <- data %>%
  mutate(condition = case_when(
    type == "H"    ~ "Kontrolle",
    type == "T3cDM" ~ "pankreopriver Diabetes",
    type == "T1DM"  ~ "Diabetes mellitus Typ1",
    TRUE ~ NA_character_   # in case other values exist
  ))

data <- data[, -1]
df <-  read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/growthresults_exchanges.tsv")
df <- df[, -1]

ef <- read_csv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/exchange_fluxes.csv")


df_unique <- df %>%
  arrange(sample_id, is.na(sex)) %>%     
  distinct(sample_id, .keep_all = TRUE) %>%
  select(sample_id, sex)


n0 <- nrow(data)
data <- data %>% left_join(df_unique, by = "sample_id")


levels_keep <- c("Diabetes mellitus Typ1", "pankreopriver Diabetes")
dat <- data %>%
  filter(condition %in% levels_keep) 
library(dplyr)

tvals <- dat %>%
  group_by(metabolite) %>%
  filter(n_distinct(condition) == 2) %>%
  group_split() %>%
  map_dfr(function(d) {
    fit <- try(lm(abundance ~ condition + sex , data = d), silent = TRUE)
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

metabolites <- c(
  "ocdcea[e]", "acac[e]", "h2[e]", "pullulan1200[e]", "csn[e]", "12dhchol[e]",
  "dgsn[e]", "akg[e]", "acnam[e]", "man[e]", "malt[e]", "12dgr180[e]",
  "fru[e]", "fum[e]", "meoh[e]", "gthrd[e]", "ptrc[e]", "ncam[e]",
  "glyclt[e]", "xan[e]","ac"
)

tvals_plot <- tvals_annot %>%
  filter(metabolite %in% metabolites) 
library(dplyr)

tvals_plot <- tvals_plot %>% drop_na(description)

tvals_plot <- tvals_plot %>%
  filter(!is.na(metabolite), !is.na(t_value)) %>%     
  mutate(metabolite = fct_reorder(metabolite, t_value, 
                                  .fun = median, .na_rm = TRUE))

top20_metabs <- tvals_annot %>%
  filter(!is.na(t_value)) %>%
  arrange(desc(abs(t_value))) %>%
  slice_head(n = 20) %>%
  pull(metabolite)

keep_set <- union(metabolites, top20_metabs)

data_kept <- tvals_annot %>%
  filter(metabolite %in% keep_set)

extra_rows <- tvals_annot %>%
  filter(description %in% c("Succinate", "Acetate","Formate"))


data_kept <- bind_rows(data_kept, extra_rows)

data_kept <- data_kept %>%
  arrange(t_value) %>%
  mutate(metabolite = factor(metabolite, levels = unique(metabolite)))

top30 <- tvals_plot %>%
  arrange(desc(abs(t_value))) %>%  slice_head(n = 25)


top30_unique <- top30 %>%
  distinct(description, .keep_all = TRUE)

lookup <- c(
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
  "Starch, structure 1 (1,6-{7[1,4-Glc], 4[1,4-Glc]})" = "Carbohydrates",
  "1,2-Diacyl-sn-glycerol (dioctadecanoyl, n-C18:0)" = "Lipids",
  "Pullulan 1200" = "Carbohydrates",
  "12-Dehydrocholic acid" = "Bile acids",
  "2-Oxoglutarate" = "Carboxylic acids",
  "Maltose C12H22O11" = "Carbohydrates",
  "1,2-Diacyl-sn-glycerol" = "Lipids",
  "Methanol" = "Other",
  "Nicotinamide" = "Nucleotides", 
  "Xanthine" = "Nucleotides"  ,
  "Formate" = "SCFA"
)

top30_unique$class <- lookup[top30_unique$description]

data_kept <- data_kept %>%
  filter(!is.na(description))

metab_list <- c( 
"Cellobiose",
"(S)-3-Methyl-2-oxopentanoate", 
"Folate",
"Octadecanoate (n-C18:0)",
"Deoxycytidine",
"Cytidine",
"Deoxyadenosine",
"Uracil",
"Pyridoxal",
"Sulfate",
"Zinc",
"Manganese",
"Magnesium",
"Potassium",
"Copper",
"Co2+",
"Chloride",
"Calcium",
"Thymidine C10H14N2O5",
"Iron (Fe3+)","N-Acetylneuraminate",
"Reduced glutathione",
"Fumarate",
"Deoxyguanosine",
"Acetate", "Formate",
"Maltose C12H22O11",
"2-Oxoglutarate",
"Putrescine","Succinate"
)

lookup <- c(
  "Reduced glutathione"            = "Amino acid related",
  "N-Acetylneuraminate"            = "Carbohydrates",
  "Sulfate"                        = "Inorganic ions",
  "Folate"                         = "Vitamins/Co-factors",
  "Cellobiose"                     = "Carbohydrates",
  "(S)-3-Methyl-2-oxopentanoate"   = "Carboxylic acids",
  "Octadecanoate (n-C18:0)"        = "Lipids",
  "Deoxycytidine"                  = "Nucleotides",
  "Cytidine"                       = "Nucleotides",
  "Deoxyadenosine"                 = "Nucleotides",
  "Uracil"                         = "Nucleotides",
  "Pyridoxal"                      = "Vitamins/Co-factors",
  "Zinc"                           = "Inorganic ions",
  "Manganese"                      = "Inorganic ions",
  "Magnesium"                      = "Inorganic ions",
  "Potassium"                      = "Inorganic ions",
  "Copper"                         = "Inorganic ions",
  "Co2+"                           = "Inorganic ions",
  "Chloride"                       = "Inorganic ions",
  "Calcium"                        = "Inorganic ions",
  "Thymidine C10H14N2O5"          = "Nucleotides",
  "Iron (Fe3+)"                    = "Inorganic ions",
  "Fumarate"                       = "Carboxylic acids",
  "Deoxyguanosine"                 = "Nucleotides",
  "Maltose C12H22O11"              = "Carbohydrates",
  "2-Oxoglutarate"                 = "Carboxylic acids",
  "Putrescine"                     = "Amino acid related",
  "Acetate" = "SCFA",
  "Succinate" = "Carboxylic acids", "Formate" ="SCFA"
)


data_kept$description_ordered <- factor(data_kept$description, levels = rev(metab_list))
data_kept$class <- lookup[data_kept$description]

data_kept <- data_kept %>%
  filter(t_value <= 2.5)
data_kept <- data_kept %>%
  distinct(description, .keep_all = TRUE)



q <- ggplot(data_kept, aes(x = t_value, y = description_ordered, fill = class)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "t value (T3cDM vs T1DM)", y = "Metabolite", fill = "Class") +
  scale_fill_viridis_d(option = "viridis") +  
  theme_minimal(base_size = 18) +  theme(
axis.text.x  = element_text(size = 16),
axis.text.y  = element_text(size = 20),
axis.title.x = element_text(size = 18, margin = margin(t = 10)),
axis.title.y = element_text(size = 18, margin = margin(r = 10)),
legend.title = element_text(size = 16),
legend.text  = element_text(size = 14),
plot.title   = element_text(size = 20, face = "bold", hjust = 0.5)
)
q


#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/metabolite_abundance_barplot_linearreg_eshighella.svg", plot = q,
#       width = 15, height = 15, units = "in", dpi = 300)#

#ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/metabolite_abundance_barplot_linearreg_eshighella.png", plot = q,
#      width = 15, height =15 , units = "in", dpi = 300)



