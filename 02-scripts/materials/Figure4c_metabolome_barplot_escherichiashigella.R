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


# Load data

data<- read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/exchange_fluxes_escherichia.csv")
df <-  read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/growthresults_exchanges.tsv")
ef <- read_csv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/exchange_fluxes.csv")

# Clean data

data <- data %>%
  mutate(condition = case_when(
    type == "H"    ~ "Kontrolle",
    type == "T3cDM" ~ "pankreopriver Diabetes",
    type == "T1DM"  ~ "Diabetes mellitus Typ1",
    TRUE ~ NA_character_   # in case other values exist
  ))

data <- data[, -1]

df <- df[, -1]

df_unique <- df %>%
  arrange(sample_id, is.na(sex)) %>%     
  distinct(sample_id, .keep_all = TRUE) %>%
  select(sample_id, sex)

n0 <- nrow(data)

data <- data %>% left_join(df_unique, by = "sample_id")

levels_keep <- c("Diabetes mellitus Typ1", "pankreopriver Diabetes")

# Log transform abundance 

data <- data %>%
  mutate(abundance_log = log10(abundance + 1e-6))

dat <- data %>%
  filter(condition %in% levels_keep) 

# By sample data

dat_sample <- dat %>%
  group_by(sample_id, metabolite, condition,sex) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE),
    flux = sum(flux, na.rm = TRUE),
    .groups = "drop"
  )


# Log transform abundance 
dat_sample <- dat_sample %>%
  mutate(abundance_log = log10(abundance + 1e-6))
library(dplyr)

## Run linear model covariate sex 
tvals <- dat_sample %>%
  group_by(metabolite) %>%
  filter(n_distinct(condition) == 2) %>%
  group_split() %>%
  map_dfr(function(d) {
    fit <- try(lm(abundance_log ~ condition + sex, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(tibble())
    
    broom::tidy(fit) %>%
      filter(term == "conditionpankreopriver Diabetes") %>%
      transmute(
        metabolite = d$metabolite[1],
        t_value = statistic,
        estimate = estimate,
        p_value = p.value
      )
  }) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

###### Map metabolite names 
get_bigg_metabolite <- function(id) {
  url <- paste0("http://bigg.ucsd.edu/api/v2/universal/metabolites/", id)
  res <- try(GET(url), silent = TRUE)
  if (inherits(res, "try-error") || res$status_code != 200) return(NULL)
  
  cont <- content(res, as = "text", encoding = "UTF-8")
  parsed <- fromJSON(cont, flatten = TRUE)
  
  # KEGG can be NULL, length 1, or length > 1
  kegg <- parsed$database_links$KEGG
  kegg_id <- if (is.null(kegg) || length(kegg) == 0) {
    NA_character_
  } else {
    paste(kegg, collapse = ";")
  }
  
  tibble(
    metabolite  = id,
    description = parsed$name %||% NA_character_,
    formula     = parsed$formula %||% NA_character_,
    kegg_id     = kegg_id,
    category    = parsed$major_bio_category %||% NA_character_
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
ids <- gsub("\\[e\\]$", "", tvals$metabolite)

## Annotate dataframe with metabolite names 
anno <- map_dfr(ids, get_bigg_metabolite)


data_keyed <- dat %>%
  mutate(bigg_id = str_remove(metabolite, "\\[[^\\]]+\\]$"))

tvals_annot <- tvals %>%
  mutate(metabolite_clean = gsub("\\[e\\]$", "", metabolite)) %>%
  left_join(anno, by = c("metabolite_clean" = "metabolite"))

top25 <- tvals_annot %>%
     slice_max(order_by = abs(t_value), n = 25)

