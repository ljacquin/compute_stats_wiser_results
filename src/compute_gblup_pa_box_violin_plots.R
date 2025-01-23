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
library(ggplot2)
library(ggpubr)
library(viridisLite)
library(dplyr)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# set color vector for gblup predictive abilities (pa) associated to wiser,
# ls-means and blup phenotype estimation
blue_col_gblup <- c("#90B3E0", "#3D9BC5", "#005AB5", "#00407A", "#002A66")[3]
yellow_col_gblup <- colorRampPalette(c("#FFEA00", "#FF7A00"))(5)[3]
green_col_gblup <- c("#A3E4A7", "#66C266", "#2E8B57", "#006400", "#003200")[3]
pa_colors_ <- c(blue_col_gblup, yellow_col_gblup, green_col_gblup)

# specify number of shuffling and folds for k-folds CV used in genomic
# prediction evaluation scheme
n_shuff_ <- 20
k_folds_ <- 5

# define species
list_species <- c("Rice", "Maize", "Apple", "Pine")

# species <- list_species[4]
for (species in list_species) {
  print(species)
  input_geno_pred_species_data_path <- paste0(
    "../data/genomic_prediction_",
    species, "/"
  )
  output_geno_pred_species_result_path <- paste0(
    "../results/reformatted_result_",
    species, "/"
  )
  # get traits associated to each species
  if (identical(species, "Rice")) {
    list_traits <- c(
      "FL", "PH", "YLD", "ZN"
    )
  } else if (identical(species, "Maize")) {
    list_traits <- c(
      "anthesis", "anthesis.silking.interval", "ear.height",
      "grain.number", "grain.weight", "grain.yield",
      "plant.height", "silking", "tassel.height"
    )
  } else if (identical(species, "Apple")) {
    list_traits <- c(
      "Color_over", "Flowering_begin", "Flowering_end",
      "Flowering_full", "Flowering_intensity", "Fruit_number",
      "Fruit_weight", "Fruit_weight_single", "Harvest_date",
      "Powdery_mildew", "Russet_freq_all", "Scab",
      "Trunk_diameter", "Trunk_increment"
    )
  } else if (identical(species, "Pine")) {
    list_traits <- c(
      "D", "H", "I", "T4", "T5", "T6"
    )
  }

  list_files_ <- list.files(input_geno_pred_species_data_path)
  list_pa_result_traits <- vector("list", length(list_traits))
  list_h2_result_traits <- vector("list", length(list_traits))

  names(list_pa_result_traits) <- list_traits
  names(list_h2_result_traits) <- list_traits

  nb_snp <- as.numeric(str_extract(list_files_[1], "(?<=_)[0-9]+(?=_SNP)"))

  for (trait_ in list_traits) {
    print(trait_)
    # get pa and h2 files for trait
    list_files_pa <- list_files_[str_detect(
      list_files_,
      "genomic_pred"
    )]
    file_pa_trait_ <- list_files_pa[str_detect(
      list_files_pa,
      paste0("SNP_", trait_, "_linear")
    )]

    list_files_h2 <- list_files_[str_detect(
      list_files_,
      "genomic_h2"
    )]
    file_h2_trait_ <- list_files_h2[str_detect(
      list_files_h2,
      paste0("SNP_", trait_, "_linear")
    )]

    # read pa files for trait
    df_pa_trait <- as.data.frame(fread(paste0(
      input_geno_pred_species_data_path,
      file_pa_trait_
    )))

    # remove statistics from distribution data, rename and select columns
    df_pa_trait <- df_pa_trait[-c(nrow(df_pa_trait) - 1, nrow(df_pa_trait)), ]
    colnames(df_pa_trait)[1] <- "Trait"
    df_pa_trait$Trait <- trait_
    df_pa_trait <- df_pa_trait[, c(
      "Trait",
      colnames(df_pa_trait)[
        str_detect(colnames(df_pa_trait), pattern = "GBLUP")
      ]
    )]

    # rename columns associated to pa results
    colnames(df_pa_trait) <- str_replace_all(colnames(df_pa_trait),
      pattern = "_wiser",
      replacement = "_PA_wiser"
    )
    colnames(df_pa_trait) <- str_replace_all(colnames(df_pa_trait),
      pattern = "_ls_means",
      replacement = "_PA_ls_means"
    )
    colnames(df_pa_trait) <- str_replace_all(colnames(df_pa_trait),
      pattern = "_blups",
      replacement = "_PA_blups"
    )

    # read pa files for trait
    df_h2_trait <- as.data.frame(fread(paste0(
      input_geno_pred_species_data_path,
      file_h2_trait_
    )))

    # remove statistics from distribution data, rename and select columns
    df_h2_trait <- df_h2_trait[-c(nrow(df_h2_trait) - 1, nrow(df_h2_trait)), ]
    colnames(df_h2_trait)[1] <- "Trait"
    df_h2_trait$Trait <- trait_
    df_h2_trait <- df_h2_trait[, c(
      "Trait",
      colnames(df_h2_trait)[
        str_detect(colnames(df_h2_trait), pattern = "GBLUP")
      ]
    )]

    # rename columns associated to h2 results
    colnames(df_h2_trait) <- str_replace_all(colnames(df_h2_trait),
      pattern = "_wiser",
      replacement = "_h2_wiser"
    )
    colnames(df_h2_trait) <- str_replace_all(colnames(df_h2_trait),
      pattern = "_ls_means",
      replacement = "_h2_ls_means"
    )
    colnames(df_h2_trait) <- str_replace_all(colnames(df_h2_trait),
      pattern = "_blups",
      replacement = "_h2_blups"
    )

    # store pa results in the list
    list_pa_result_traits[[trait_]] <- df_pa_trait

    # store h2 results in the list
    list_h2_result_traits[[trait_]] <- df_h2_trait
  }
  # end of computation for species

  # rbind pa list
  df_pa_result_traits <- do.call(
    rbind,
    list_pa_result_traits
  )
  colnames(df_pa_result_traits) <- c(
    "Trait",
    "1. GBLUP PA for WISER phenotypes",
    "2. GBLUP PA for LS-means phenotypes",
    "3. GBLUP PA for BLUP phenotypes"
  )

  # rbind h2 list
  df_h2_result_traits <- do.call(
    rbind,
    list_h2_result_traits
  )
  colnames(df_h2_result_traits) <- c(
    "Trait",
    "1. GBLUP h2 for WISER phenotypes",
    "2. GBLUP h2 for LS-means phenotypes",
    "3. GBLUP h2 for BLUP phenotypes"
  )

  # convert pa results to long format
  df_long_pa_result_traits <- df_pa_result_traits %>%
    pivot_longer(
      cols = starts_with(c("1. GBLUP", "2. GBLUP", "3. GBLUP")),
      names_to = "Method",
      values_to = "PA"
    )

  # convert h2 results to long format
  df_long_h2_result_traits <- df_h2_result_traits %>%
    pivot_longer(
      cols = starts_with(c("1. GBLUP", "2. GBLUP", "3. GBLUP")),
      names_to = "Method",
      values_to = "h2"
    )

  # make ggplot for gblup pa across traits and phenotypic estimation methods

  # create the dynamic title
  title_ <- paste0(
    species,
    " GBLUP predictive ability (PA) for traits across WISER, LS-means and BLUP phenotypes, \n",
    "based on ", nb_snp, " SNP across ", n_shuff_,
    " shuffling scenarios for ", k_folds_, "-fold cross-validation"
  )

  ggplot_obj_ <- ggplot(
    data = df_long_pa_result_traits,
    aes(x = Trait, y = PA, fill = Method)
  ) +
    scale_fill_manual(values = pa_colors_) + # Use custom colors
    geom_violin(
      alpha = 0.5, position = position_dodge(width = .75),
      size = 1, color = NA
    ) +
    geom_boxplot(
      notch = TRUE, outlier.size = -1, color = "black", lwd = 1,
      alpha = 0.7, show.legend = F
    ) +
    scale_y_continuous(
      limits = c(min(df_long_pa_result_traits$PA), 1), # set y-axis range
      breaks = seq(min(df_long_pa_result_traits$PA), 1, by = 0.1) # add more ticks with 0.1 intervals
    ) +
    theme_minimal(base_family = "Calibri") +
    ylab("GBLUP PA") +
    xlab(NULL) +
    rremove("legend.title") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.line = element_line(colour = "black", size = 1),
      axis.ticks = element_line(size = 1, color = "black"),
      axis.text = element_text(color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # rotate x-axis labels
      legend.position = c(0.45, 0.9)
    ) +
    font("xylab", size = 15) +
    font("xy", size = 15) +
    font("xy.text", size = 15) +
    font("legend.text", size = 15) +
    guides(fill = guide_legend(override.aes = list(alpha = 1, color = "black"))) +
    ggtitle(title_)

  # save the plot as png
  ggsave(
    filename = paste0(
      output_geno_pred_species_result_path, "/gblup_pa_", species, ".png"
    ),
    plot = ggplot_obj_,
    width = 18,
    height = 8,
    dpi = 300
  )
}
