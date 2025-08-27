
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(readr)

library(dplyr)
library(purrr)
library(rstatix)

df <- read_tsv("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/materials/growthresults_exchanges.tsv")
df <- df[, -1]



df <- df %>%
  filter(abundance >= 0.5 )




edges <- df %>%
  mutate(abundance = suppressWarnings(as.numeric(abundance))) %>%
  filter(is.finite(abundance)) %>%
  group_by(taxon, metabolite) %>%
  summarise(value = mean(abundance), .groups = "drop") %>%
  group_by(taxon) %>%
  slice_max(value, n = 10, with_ties = FALSE) %>%
  ungroup()

library(dplyr)
library(stringr)

translate_metabolites <- function(x) {
  # strip compartment like [e]
  id <- str_to_lower(x) %>% str_remove_all("\\[[^\\]]+\\]")
  
  # small dictionary just for your set (edit/extend anytime)
  dict <- c(

    "glc_d"       = "D-glucose",
    "h2"          = "Hydrogen",
    "hdca"        = "Hexadecanoate (palmitate)",
    "12ppd_s"     = "(S)-1,2-propanediol",
    "26dap_m"     = "meso-2,6-diaminopimelate",
    "2obut"       = "2-oxobutanoate (α-ketobutyrate)",
    "ac"          = "Acetate",
    "acald"       = "Acetaldehyde",
    "acgam"       = "N-acetyl-D-glucosamine",
    "4hbz"        = "4-hydroxybenzoate",
    "ade"         = "Adenine",
    "adn"         = "Adenosine",
    "ala_d"       = "D-alanine",
    "ala_l"       = "L-alanine",
    "urea"        = "Urea",
    "amp"         = "AMP (adenosine monophosphate)",
    "for"         = "Formate",
    "cgly"        = "Cysteinylglycine",
    "gal"         = "D-galactose",
    "galmannan"   = "Galactomannan",
    "pullulan1200"= "Pullulan (~1200 Da)",
    "xyluglc"     = "Xyloglucan",
    "3mop"        = "3-methyl-2-oxopentanoate",
    "asn_l"       = "L-asparagine",
    "asp_l"       = "L-aspartate",
    "but"         = "Butyrate",
    "gam"  = "D-Galactosamine",
    "ncam" = "ncam"
  )
  
  out <- dict[id]
  
  # fallback: make something readable if not in dict yet
  missing <- is.na(out)
  if (any(missing)) {
    out[missing] <- str_to_sentence(str_replace_all(id[missing], "_", "-"))
  }
  out
}

# Apply to your edges tibble/data.frame
edges <- edges %>%
  mutate(metabolite_pretty = translate_metabolites(metabolite))


g <- ggplot(edges,
            aes(axis1 = taxon, axis2 = metabolite_pretty, y = value)) +
  geom_alluvium(aes(fill = taxon), alpha = 0.6) +
  geom_stratum(width = 0.15) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 10) +
  scale_x_discrete(limits = c("Taxon", "Metabolite"), expand = c(.15, .05)) +
  labs(y = "Abundance (mean)", x = NULL) +
  theme_minimal(base_size = 10) +        # set default font size
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 40),
    plot.title   = element_text(size = 10),
    strip.text   = element_text(size = 10)
  )

g
ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/sankey_microbiome_metaoblite.svg", plot = g,
       width = 20, height = 20, units = "in", dpi = 300)

ggsave("/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/sankey_microbiome_metaoblite.png", plot = g,
       width = 20, height = 20, units = "in", dpi = 300)
