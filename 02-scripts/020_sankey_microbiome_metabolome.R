
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(readr)

library(dplyr)
library(purrr)
library(rstatix)

df <- read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/exchange_fluxes_all_taxa.csv")
df <- df[, -1]



df <- df %>%
  filter(abundance >= 0.4)

taxa_keep <- c("Akkermansia",
               "Blautia",
               "Escherichia",
               "Streptococcus",
               "Faecalibacterium",
               "Subdoligranulum")  # Fusobacterium omitted

df <- df %>%
  filter(taxon %in% taxa_keep)


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


ids <- gsub("\\[e\\]$", "", df$metabolite)


anno <- map_dfr(ids, get_bigg_metabolite)

data_keyed <- df %>%
  mutate(bigg_id = str_remove(metabolite, "\\[[^\\]]+\\]$"))

df <- read_tsv("/data/scratch/kvalem/projects/2024/diabetes_microbe/01-tables/materials/exchange_fluxes_all_taxa.csv")
df <- df[, -1]



df <- df %>%
  filter(abundance >= 0.1)

taxa_keep <- c("Akkermansia",
               "Blautia",
               "Escherichia",
               "Streptococcus",
               "Faecalibacterium",
               "Subdoligranulum")  # Fusobacterium omitted

df <- df %>%
  filter(taxon %in% taxa_keep)

df <- df %>%
  mutate(metabolite_clean = gsub("\\[e\\]$", "", metabolite)) %>%
  left_join(anno, by = c("metabolite_clean" = "metabolite"))


edges <- df %>%
  mutate(abundance = suppressWarnings(as.numeric(abundance))) %>%
  filter(is.finite(abundance)) %>%
  group_by(taxon, description) %>%
  summarise(value = mean(abundance), .groups = "drop") %>%
  group_by(taxon) %>%
  slice_max(value, n = 10, with_ties = FALSE) %>%
  ungroup()

library(dplyr)
library(stringr)





g <- ggplot(edges,
            aes(axis1 = taxon, axis2 = description, y = value)) +
  geom_alluvium(aes(fill = taxon), alpha = 0.6) +
  geom_stratum(width = 0.15) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("", ""), expand = c(.15, .05)) +
  labs(y = "", x = NULL) +
  theme_minimal(base_size = 5) +        # set default font size
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 40),
    plot.title   = element_text(size = 10),
    strip.text   = element_text(size = 10)
  )

g

ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/sankey_microbiome_metaoblite.svg", plot = g,
       width = 10, height = 10, units = "in", dpi = 300)


ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/sankey_microbiome_metaoblite.png", plot = g,
       width =10, height = 10, units = "in", dpi = 300)
