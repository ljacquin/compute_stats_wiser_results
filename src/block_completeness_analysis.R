# script meant to compute stats for pa and h2 results
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
library(devtools)
# install other requirements from github if necessary
if ("refpop_env" %in% conda_list()$name) {
  use_condaenv("refpop_env")
}
library(MASS)
library(data.table)
library(stringr)
library(tidyr)
library(Matrix)
library(matrixcalc)
library(ggplot2)
library(corrplot)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(stats)
library(fpc)
library(umap)
library(future)
library(future.apply)
library(adegenet)
library(hierfstat)
library(parallel)
computation_mode <- "local"
if (!identical(computation_mode, "cluster")) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}
source("functions.R")

# get phenotypic data for rice, maize and pine datasets which contain blocks
df_rice <- as.data.frame(
  fread("../data/phenotype_data_Rice.csv")
)
length(unique(df_rice$Envir))

df_maize <- as.data.frame(
  fread("../data/phenotype_data_Maize/phenotype_data_Maize.csv")
)
length(unique(df_maize$Envir))

df_pine <- as.data.frame(
  fread("../data/phenotype_data_Pine.csv")
)
length(unique(df_pine$Envir))

df_apple <- as.data.frame(
  fread("../data/phenotype_data_Apple.csv")
)
length(unique(df_apple$Envir))

# rice
# get expected number of genotypes per site
site_geno_counts <- df_rice %>%
  group_by(Site) %>%
  summarise(n_geno_per_site = n_distinct(Genotype), .groups = "drop")

# associate each Envir to its site and count genotypes per Envir
rice_block_completeness <- df_rice %>%
  group_by(Site, Envir, BLOC) %>%
  summarise(
    n_geno_per_block_in_envir = n_distinct(Genotype),
    .groups = "drop"
  ) %>%
  left_join(site_geno_counts, by = "Site") %>%
  mutate(
    is_block_in_envir_complete = n_geno_per_block_in_envir == n_geno_per_site
  )
rice_block_completeness <- as.data.frame(rice_block_completeness)
rice_block_completeness <- rice_block_completeness[
  , c(
    "Site", "Envir", "BLOC", "n_geno_per_site",
    "n_geno_per_block_in_envir",
    "is_block_in_envir_complete"
  )
]
sum(!rice_block_completeness$is_block_in_envir_complete)

# write results for block completeness test
fwrite(rice_block_completeness, "../results/rice_block_completeness_test.csv")

# maize
# get expected number of genotypes per site
site_geno_counts <- df_maize %>%
  group_by(Site) %>%
  summarise(n_geno_per_site = n_distinct(Genotype), .groups = "drop")

# associate each Envir to its site and count genotypes per Envir
maize_block_completeness <- df_maize %>%
  group_by(Site, Envir, block) %>%
  summarise(
    n_geno_per_block_in_envir = n_distinct(Genotype),
    .groups = "drop"
  ) %>%
  left_join(site_geno_counts, by = "Site") %>%
  mutate(
    is_block_in_envir_complete = n_geno_per_block_in_envir == n_geno_per_site
  )
maize_block_completeness <- as.data.frame(maize_block_completeness)
maize_block_completeness <- maize_block_completeness[
  , c(
    "Site", "Envir", "block", "n_geno_per_site",
    "n_geno_per_block_in_envir",
    "is_block_in_envir_complete"
  )
]
sum(!maize_block_completeness$is_block_in_envir_complete)

# write results for block completeness test
fwrite(maize_block_completeness, "../results/maize_block_completeness_test.csv")

# pine
# get expected number of genotypes per site
site_geno_counts <- df_pine %>%
  group_by(Site) %>%
  summarise(n_geno_per_site = n_distinct(Genotype), .groups = "drop")

# associate each Envir to its site and count genotypes per Envir
pine_block_completeness <- df_pine %>%
  group_by(Site, Envir, Block) %>%
  summarise(
    n_geno_per_block_in_envir = n_distinct(Genotype),
    .groups = "drop"
  ) %>%
  left_join(site_geno_counts, by = "Site") %>%
  mutate(
    is_block_in_envir_complete = n_geno_per_block_in_envir == n_geno_per_site
  )
pine_block_completeness <- as.data.frame(pine_block_completeness)
pine_block_completeness <- pine_block_completeness[
  , c(
    "Site", "Envir", "Block", "n_geno_per_site",
    "n_geno_per_block_in_envir",
    "is_block_in_envir_complete"
  )
]
sum(!pine_block_completeness$is_block_in_envir_complete)

# write results for block completeness test
fwrite(pine_block_completeness, "../results/pine_block_completeness_test.csv")

# apple
# get expected number of genotypes per site
site_geno_counts <- df_apple %>%
  group_by(Site) %>%
  summarise(n_geno_per_site = n_distinct(Genotype), .groups = "drop")

# associate each Envir to its site and count genotypes per Envir
apple_manage_completeness <- df_apple %>%
  group_by(Site, Envir, Management) %>%
  summarise(
    n_geno_per_manage_in_envir = n_distinct(Genotype),
    .groups = "drop"
  ) %>%
  left_join(site_geno_counts, by = "Site") %>%
  mutate(
    is_manage_in_envir_complete = n_geno_per_manage_in_envir == n_geno_per_site
  )
apple_manage_completeness <- as.data.frame(apple_manage_completeness)
apple_manage_completeness <- apple_manage_completeness[
  , c(
    "Site", "Envir", "Management", "n_geno_per_site",
    "n_geno_per_manage_in_envir",
    "is_manage_in_envir_complete"
  )
]
sum(!apple_manage_completeness$is_manage_in_envir_complete)

# write results for manage completeness test
fwrite(apple_manage_completeness, "../results/apple_manage_completeness_test.csv")
