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
computation_mode <- "cluster"
if (!identical(computation_mode, "cluster")) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}
source("functions.R")

# define number of cores for parallel computations
nb_cores_ <- 12

# specify number of shuffling and folds for k-folds CV used in genomic
# prediction evaluation scheme
n_shuff_ <- 20
k_folds_ <- 5
snp_sample_size_ <- 50e3 # sample same number of markers for Apple as in the study
nb_cluster_ <- 30
nrow_lim <- 5e3

# umap parameters, most sensitive ones
random_state_umap_ <- 15
n_neighbors_umap_ <- 10
min_dist_ <- 0.1
nb_comp_umap_ <- 5

# define species and input path
list_species <- c("Rice", "Maize", "Apple", "Pine")

# initialize empty data frames used for combining results into list
df_pa_median_iqr_all_species_traits <- data.frame()
df_h2_result_all_species_traits <- data.frame()
df_geno_entropy_results_all_species_traits <- data.frame()
df_geno_site_entropy_results_all_species_traits <- data.frame()
df_geno_envir_entropy_results_all_species_traits <- data.frame()
df_kmeans_opt_k_result_all_species_traits <- data.frame()
df_kmeans_opt_ch_result_all_species_traits <- data.frame()
df_kmeans_fst_result_all_species_traits <- data.frame()

# iterate over list of species
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
  list_geno_entropy_result_traits <- vector("list", length(list_traits))
  list_geno_site_entropy_result_traits <- vector("list", length(list_traits))
  list_geno_envir_entropy_result_traits <- vector("list", length(list_traits))
  list_kmeans_opt_k_result_traits <- vector("list", length(list_traits))
  list_kmeans_opt_ch_result_traits <- vector("list", length(list_traits))
  list_kmeans_fst_result_traits <- vector("list", length(list_traits))

  names(list_pa_result_traits) <- list_traits
  names(list_h2_result_traits) <- list_traits
  names(list_geno_entropy_result_traits) <- list_traits
  names(list_geno_site_entropy_result_traits) <- list_traits
  names(list_geno_envir_entropy_result_traits) <- list_traits
  names(list_kmeans_opt_k_result_traits) <- list_traits
  names(list_kmeans_opt_ch_result_traits) <- list_traits
  names(list_kmeans_fst_result_traits) <- list_traits

  # extract number of snps for plots with plotly
  nb_snp <- as.numeric(str_extract(list_files_[1], "(?<=_)[0-9]+(?=_SNP)"))

  # get phenotype data associated to species
  df_pheno_species <- as.data.frame(fread(
    paste0(
      "../data/phenotype_data_",
      species, ".csv"
    )
  ))

  # get genomic data associated to species
  df_genomic_species <- as.data.frame(fread(
    paste0(
      "../data/genomic_data_",
      species, ".csv"
    )
  ))
  colnames(df_genomic_species)[1] <- "Genotype"

  # sample markers for apple dataset to make analysis comparable to study
  if (species == "Apple") {
    set.seed(123)
    idx_snp_sample_size_ <- sample(
      2:ncol(df_genomic_species),
      size = snp_sample_size_, replace = F
    )
    df_genomic_species <- df_genomic_species[
      , c(1, idx_snp_sample_size_)
    ]
  }

  # compute umap for genomic data
  geno_umap_ <- data.frame(
    umap(apply(df_genomic_species[, -1], 2, as.numeric),
      n_components = nb_comp_umap_, random_state = random_state_umap_,
      n_neighbors = n_neighbors_umap_, min_dist = min_dist_
    )[["layout"]]
  )
  geno_umap_$Genotype <- df_genomic_species$Genotype

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

    # read pa and h2 files for trait, and format results

    #- pa
    df_pa_trait <- as.data.frame(fread(paste0(
      input_geno_pred_species_data_path,
      file_pa_trait_
    )))
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
    df_pa_trait_stats <- as.data.frame(tail(df_pa_trait, 2)[, -1])

    # format values to two decimal places before combining
    df_pa_trait_stats[1, ] <- sprintf(
      "%.2f",
      as.numeric(df_pa_trait_stats[1, ])
    )
    df_pa_trait_stats[2, ] <- sprintf(
      "%.2f",
      as.numeric(df_pa_trait_stats[2, ])
    )

    # combine rows into the desired format
    df_pa_trait_stats[1, ] <- paste0(
      df_pa_trait_stats[1, ], " (",
      df_pa_trait_stats[2, ], ")"
    )
    df_pa_trait_stats <- df_pa_trait_stats[-2, ]
    rownames(df_pa_trait_stats) <- NULL

    # add trait column
    df_pa_trait_stats$Trait <- trait_

    # store result in the list
    list_pa_result_traits[[trait_]] <- df_pa_trait_stats

    #- h2
    df_h2_trait <- as.data.frame(fread(paste0(
      input_geno_pred_species_data_path,
      file_h2_trait_
    )))
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
    df_h2_trait_stats <- as.data.frame(tail(df_h2_trait, 2)[, -1])

    # format values to two decimal places
    df_h2_trait_stats[1, ] <- sprintf(
      "%.2f",
      as.numeric(df_h2_trait_stats[1, ])
    )
    df_h2_trait_stats[2, ] <- sprintf(
      "%.2f",
      as.numeric(df_h2_trait_stats[2, ])
    )

    # combine rows into the desired format
    df_h2_trait_stats[1, ] <- paste0(
      df_h2_trait_stats[1, ], " (",
      df_h2_trait_stats[2, ], ")"
    )
    df_h2_trait_stats <- df_h2_trait_stats[-2, ]
    rownames(df_h2_trait_stats) <- NULL

    # add trait column
    df_h2_trait_stats$Trait <- trait_

    # store result in the list
    list_h2_result_traits[[trait_]] <- df_h2_trait_stats

    # get genotype, site and envir data for trait
    df_trait_genotype_site_envir <- df_pheno_species[
      , c(trait_, "Genotype", "Site", "Envir")
    ]
    df_trait_genotype_site_envir$Genotype <- as.factor(
      df_trait_genotype_site_envir$Genotype
    )

    # keep genotype, site and envir values for which individual phenotypic
    # trait value is non missing
    df_trait_genotype_site_envir <- df_trait_genotype_site_envir[
      !is.na(suppressWarnings(
        as.numeric(df_trait_genotype_site_envir[, trait_])
      )),
    ]

    # get genomic data for genotypes for which individual phenotypic
    # trait value is non missing
    df_trait_genomic <- df_genomic_species[
      df_genomic_species$Genotype %in% df_trait_genotype_site_envir$Genotype,
    ]
    df_trait_geno_umap <- geno_umap_[
      geno_umap_$Genotype %in% df_trait_genotype_site_envir$Genotype,
    ]

    # ascertain that df_trait_genotype_site_envir has same genotypes as df_trait_genomic
    df_trait_genotype_site_envir <- df_trait_genotype_site_envir[
      df_trait_genotype_site_envir$Genotype %in% df_trait_genomic$Genotype,
    ]
    df_trait_genotype_site_envir <- droplevels(df_trait_genotype_site_envir)

    tryCatch(
      {
        # expand df_trait_genomic and geno_umap_ according to repetitions
        # in the experimental design
        repetition_counts <- table(df_trait_genotype_site_envir$Genotype)

        df_trait_genomic_expanded <- df_trait_genomic[
          rep(1:nrow(df_trait_genomic),
            times = repetition_counts[as.character(df_trait_genomic$Genotype)]
          ),
        ]

        df_trait_geno_umap_expanded <- df_trait_geno_umap[
          rep(1:nrow(df_trait_geno_umap),
            times = repetition_counts[as.character(df_trait_geno_umap$Genotype)]
          ),
        ]

        # if expanded datasets are too large, apply a stratified sampling to preserve
        # data structure between genotypes and reduce computation burden linked
        # to memory and ressources (limit number of rows to nrow_lim approx.)
        if (nrow(df_trait_genomic_expanded) > nrow_lim) {
          df_trait_genomic_stratified_sampled <- df_trait_genomic_expanded %>%
            group_by(across(1)) %>% # group by genotype column
            sample_frac(size = nrow_lim / nrow(df_trait_genomic_expanded))
          df_trait_genomic_stratified_sampled <- as.data.frame(
            df_trait_genomic_stratified_sampled
          )

          df_trait_geno_umap_stratified_sampled <- df_trait_geno_umap_expanded %>%
            group_by(across(1)) %>% # group by genotype column
            sample_frac(size = nrow_lim / nrow(df_trait_geno_umap_expanded))
          df_trait_geno_umap_stratified_sampled <- as.data.frame(
            df_trait_geno_umap_stratified_sampled
          )
        } else {
          df_trait_genomic_stratified_sampled <- df_trait_genomic_expanded
          df_trait_geno_umap_stratified_sampled <- df_trait_geno_umap_expanded
        }

        # estimate "optimal" number of cluster (k), for trait-species umap data,
        # using k-means
        umap_kmeans_opt_ <- compute_ch_opt_k_parallel(
          X = apply(df_trait_geno_umap_stratified_sampled[
            ,
            -ncol(df_trait_geno_umap_stratified_sampled)
          ], 2, as.numeric),
          k_range = 2:nb_cluster_
        )

        # replace fist column with cluster indices
        df_trait_genomic_stratified_sampled[, 1] <- umap_kmeans_opt_$opt_clust_vect

        # rename columns for fst computation
        colnames(df_trait_genomic_stratified_sampled) <- c(
          "pop",
          paste0("M", 1:(ncol(df_trait_genomic_stratified_sampled) - 1))
        )

        # split genomic data matrix into K sub-matrices for fst calculation
        # NB. K = (ncol(df_trait_genomic_stratified_sampled)-1)/chunk_size
        list_sub_matrices <- split_matrix(
          df_trait_genomic_stratified_sampled,
          chunk_size = 1000
        )

        if (!identical(computation_mode, "cluster")) {
          cl <- makeCluster(nb_cores_)
          clusterEvalQ(cl, {
            library(adegenet)
            library(hierfstat)
          })
          clusterExport(cl,
            varlist = c("calculate_fst", "list_sub_matrices"),
            envir = environment()
          )
          fst_values <- parLapply(cl, list_sub_matrices, calculate_fst)
          stopCluster(cl)
        } else {
          fst_values <- mclapply(list_sub_matrices, calculate_fst,
            mc.cores = nb_cores_
          )
        }
        fst_ <- mean(unlist(fst_values))

        # assign fst value to associated list of results
        list_kmeans_fst_result_traits[[
          trait_
        ]] <- fst_
      },
      error = function(e) {
        cat(
          "Error with : ", conditionMessage(e), "\n"
        )
      }
    )

    # assign "optimal" number of cluster (k) to associated list of results
    list_kmeans_opt_k_result_traits[[
      trait_
    ]] <- umap_kmeans_opt_$opt_k

    # assign "optimal" ch value to associated list of results
    list_kmeans_opt_ch_result_traits[[
      trait_
    ]] <- umap_kmeans_opt_$opt_ch

    # compute entropy of genotype frequencies
    list_geno_entropy_result_traits[[
      trait_
    ]] <- compute_entropy_factor(df_trait_genotype_site_envir$Genotype)

    # compute joint entropy between genotype and site
    list_geno_site_entropy_result_traits[[trait_]] <- compute_joint_entropy_factors(
      df_trait_genotype_site_envir$Genotype,
      df_trait_genotype_site_envir$Site
    )

    # compute joint entropy between genotype and environment
    list_geno_envir_entropy_result_traits[[trait_]] <- compute_joint_entropy_factors(
      df_trait_genotype_site_envir$Genotype,
      df_trait_genotype_site_envir$Envir
    )
  }

  # rbind pa list across traits
  df_pa_result_traits <- do.call(
    rbind,
    list_pa_result_traits
  )

  # add species and make last two columns become two first ones
  df_pa_result_traits$Species <- species
  df_pa_result_traits <- df_pa_result_traits[, c(
    ncol(df_pa_result_traits), ncol(df_pa_result_traits) - 1,
    1:(ncol(df_pa_result_traits) - 2)
  )]

  # get gblup pa results for species and traits across phenotype estimation methods
  df_gblup_pa_species_traits <- df_pa_result_traits[, str_detect(
    colnames(df_pa_result_traits), "Species|Trait|GBLUP"
  )]

  # rbind h2 list (note: only h2 for gblup have been computed and reported)
  df_gblup_h2_species_traits <- do.call(
    rbind,
    list_h2_result_traits
  )
  df_gblup_h2_species_traits$Species <- species

  # add species and make last two columns become two first ones
  df_gblup_h2_species_traits <- df_gblup_h2_species_traits[, c(
    ncol(df_gblup_h2_species_traits), ncol(df_gblup_h2_species_traits) - 1,
    1:(ncol(df_gblup_h2_species_traits) - 2)
  )]

  # convert genotype entropy list into data frame
  df_geno_entropy_results_traits <- data.frame(
    "Species" = species,
    "Trait" = names(list_geno_entropy_result_traits),
    "Genotype_entropy" = unlist(list_geno_entropy_result_traits)
  )

  # convert genotype-site joint entropy list into data frame
  df_geno_site_entropy_results_traits <- data.frame(
    "Species" = species,
    "Trait" = names(list_geno_site_entropy_result_traits),
    "Genotype_site_entropy" = unlist(list_geno_site_entropy_result_traits)
  )

  # convert genotype-envir joint entropy list into data frame
  df_geno_envir_entropy_results_traits <- data.frame(
    "Species" = species,
    "Trait" = names(list_geno_envir_entropy_result_traits),
    "Genotype_envir_entropy" = unlist(list_geno_envir_entropy_result_traits)
  )

  # convert detected number of clusters list into data frame
  df_kmeans_opt_k_results_traits <- data.frame(
    "Species" = species,
    "Trait" = names(list_kmeans_opt_k_result_traits),
    "Number_of_clusters" = unlist(list_kmeans_opt_k_result_traits)
  )

  # convert computed CH value list into data frame
  df_kmeans_opt_ch_results_traits <- data.frame(
    "Species" = species,
    "Trait" = names(list_kmeans_opt_ch_result_traits),
    "CH_value" = unlist(list_kmeans_opt_ch_result_traits)
  )

  # convert computed fst value list into data frame
  df_kmeans_fst_results_traits <- data.frame(
    "Species" = species,
    "Trait" = names(list_kmeans_fst_result_traits),
    "Fst_value" = unlist(list_kmeans_fst_result_traits)
  )

  # write pa and h2 results for all predictive models,
  # and for gblup specifically
  fwrite(df_pa_result_traits, file = paste0(
    output_geno_pred_species_result_path,
    "all_predictive_models_median_iqr_pa_", species,
    ".csv"
  ))
  fwrite(df_gblup_pa_species_traits, file = paste0(
    output_geno_pred_species_result_path,
    "gblup_median_iqr_pa_", species,
    ".csv"
  ))
  fwrite(df_gblup_h2_species_traits, file = paste0(
    output_geno_pred_species_result_path,
    "gblup_median_iqr_h2_", species,
    ".csv"
  ))

  # for each species, extract gblup pa medians across traits and methods
  df_median_pa <- df_gblup_pa_species_traits %>%
    mutate(
      GBLUP_wiser = as.numeric(sub(" \\(.*", "", GBLUP_PA_wiser)),
      GBLUP_ls_means = as.numeric(sub(" \\(.*", "", GBLUP_PA_ls_means)),
      GBLUP_blups = as.numeric(sub(" \\(.*", "", GBLUP_PA_blups))
    ) %>%
    select(Trait, GBLUP_wiser, GBLUP_ls_means, GBLUP_blups) %>%
    tidyr::pivot_longer(-Trait, names_to = "Method", values_to = "Median")

  # change labels
  df_median_pa$Method[
    df_median_pa$Method == "GBLUP_blups"
  ] <- "GBLUP median PA for BLUP phenotypes"
  df_median_pa$Method[
    df_median_pa$Method == "GBLUP_ls_means"
  ] <- "GBLUP median PA for LS-means phenotypes"
  df_median_pa$Method[
    df_median_pa$Method == "GBLUP_wiser"
  ] <- "GBLUP median PA for WISER phenotypes"

  median_pa_plot_ <- plot_ly(df_median_pa,
    x = ~Trait, y = ~Median, color = ~Method,
    colors = c("#2E8B57", "#FFB200", "#005AB5"),
    type = "scatter", mode = "lines+markers"
  ) %>%
    layout(
      title = paste0(
        species,
        " GBLUP median predictive ability (PA) for traits across BLUP, LS-means, and WISER phenotypes,
      based on ", nb_snp, " SNP across ", n_shuff_,
        " shuffling scenarios for ", k_folds_, "-fold cross-validation"
      ),
      xaxis = list(title = "trait"),
      yaxis = list(title = "GBLUP median PA"),
      legend = list(
        x = 1, # always on the right
        y = 0.98, # slightly lower to avoid overlapping with the title
        xanchor = "right",
        yanchor = "top",
        bgcolor = "rgba(255, 255, 255, 0.7)" # semi-transparent background for better readability
      )
    )
  saveWidget(median_pa_plot_, file = paste0(
    output_geno_pred_species_result_path, "/gblup_median_iqr_pa_", species, ".html"
  ))

  # get gblup median and iqr pa for all species and traits
  df_pa_median_iqr_all_species_traits <- rbind(
    df_pa_median_iqr_all_species_traits, df_gblup_pa_species_traits
  )
  df_h2_result_all_species_traits <- rbind(
    df_h2_result_all_species_traits, df_gblup_h2_species_traits
  )

  # get detected number of clusters for all species and traits
  df_kmeans_opt_k_result_all_species_traits <- rbind(
    df_kmeans_opt_k_result_all_species_traits,
    df_kmeans_opt_k_results_traits
  )

  # get computed CH value for all species and traits
  df_kmeans_opt_ch_result_all_species_traits <- rbind(
    df_kmeans_opt_ch_result_all_species_traits,
    df_kmeans_opt_ch_results_traits
  )

  # get computed fst value for all species and traits
  df_kmeans_fst_result_all_species_traits <- rbind(
    df_kmeans_fst_result_all_species_traits,
    df_kmeans_fst_results_traits
  )

  # get genotype entropy for all species and traits
  df_geno_entropy_results_all_species_traits <- rbind(
    df_geno_entropy_results_all_species_traits,
    df_geno_entropy_results_traits
  )

  # get genotype-site joint entropy for all species and traits
  df_geno_site_entropy_results_all_species_traits <- rbind(
    df_geno_site_entropy_results_all_species_traits,
    df_geno_site_entropy_results_traits
  )

  # get genotype-envir joint entropy for all species and traits
  df_geno_envir_entropy_results_all_species_traits <- rbind(
    df_geno_envir_entropy_results_all_species_traits,
    df_geno_envir_entropy_results_traits
  )
}
# end of computation for species

# data treatment section
colnames(df_pa_median_iqr_all_species_traits) <- c(
  "Species",
  "Trait",
  "GBLUP_median_PA_wiser",
  "GBLUP_median_PA_ls_means",
  "GBLUP_median_PA_blups"
)
colnames(df_h2_result_all_species_traits) <- c(
  "Species",
  "Trait",
  "GBLUP_median_h2_wiser",
  "GBLUP_median_h2_ls_means",
  "GBLUP_median_h2_blups"
)

# merge data frames to garantee data analyses integrity
merged_data <- merge(
  df_pa_median_iqr_all_species_traits,
  df_h2_result_all_species_traits,
  by = c("Species", "Trait")
)
merged_data <- merge(merged_data,
  df_geno_entropy_results_all_species_traits,
  by = c("Species", "Trait")
)
merged_data <- merge(merged_data,
  df_geno_site_entropy_results_all_species_traits,
  by = c("Species", "Trait")
)
merged_data <- merge(merged_data,
  df_geno_envir_entropy_results_all_species_traits,
  by = c("Species", "Trait")
)
merged_data <- merge(merged_data,
  df_kmeans_opt_k_result_all_species_traits,
  by = c("Species", "Trait")
)
merged_data <- merge(merged_data,
  df_kmeans_opt_ch_result_all_species_traits,
  by = c("Species", "Trait")
)
merged_data <- merge(merged_data,
  df_kmeans_fst_result_all_species_traits,
  by = c("Species", "Trait")
)

# extract gblup pa medians and compute differences between these across all
# species and traits
df_pa_h2_diff_and_stats <- merged_data %>%
  mutate(
    GBLUP_median_PA_wiser = as.numeric(sub(" \\(.*", "", GBLUP_median_PA_wiser)),
    GBLUP_median_PA_ls_means = as.numeric(sub(" \\(.*", "", GBLUP_median_PA_ls_means)),
    GBLUP_median_PA_gain = GBLUP_median_PA_wiser - GBLUP_median_PA_ls_means,
    GBLUP_median_h2_wiser = as.numeric(sub(" \\(.*", "", GBLUP_median_h2_wiser)),
    GBLUP_median_h2_ls_means = as.numeric(sub(" \\(.*", "", GBLUP_median_h2_ls_means)),
    GBLUP_median_h2_difference = GBLUP_median_h2_wiser - GBLUP_median_h2_ls_means
  ) %>%
  select(
    Species,
    Trait,
    GBLUP_median_PA_gain,
    GBLUP_median_h2_difference,
    GBLUP_median_h2_wiser,
    Genotype_entropy,
    Genotype_site_entropy,
    Genotype_envir_entropy,
    CH_value,
    Fst_value
  )
fwrite(
  df_pa_h2_diff_and_stats,
  "../results/gblup_median_pa_gain_h2_deviations_and_computed_stats.csv"
)

# df_pa_h2_diff_and_stats <- as.data.frame(fread(
#   "../results/gblup_median_pa_gain_h2_deviations_and_computed_stats.csv"
# ))

# compute average of median PA gain across all species and traits
mean(df_pa_h2_diff_and_stats$GBLUP_median_PA_gain)

# compute median and iqr
median_iqr_df <- as.data.frame(df_pa_h2_diff_and_stats %>%
  group_by(Species) %>%
  summarise(
    Genotype_entropy_median = round(median(Genotype_entropy), 3),
    Genotype_entropy_iqr = round(IQR(Genotype_entropy), 3),
    Fst_value_median = round(median(Fst_value), 3),
    Fst_value_iqr = round(IQR(Fst_value), 3)
  ))
median_iqr_df

# define variables to keep for correlation analyses
vars_to_keep <- c(
  "GBLUP_median_PA_gain",
  "GBLUP_median_h2_difference",
  "Genotype_entropy",
  "Genotype_site_entropy",
  "Genotype_envir_entropy",
  "Fst_value"
)

# plot a correlation matrix for df_pa_h2_diff_and_stats, with values inside the circles
# and only the lower triangular part without the diagonal
df_pa_h2_diff_and_stats_all_traits <- df_pa_h2_diff_and_stats[, vars_to_keep]
corr_matrix <- cor(df_pa_h2_diff_and_stats_all_traits)

# start the pdf device to save the plot
pdf("../results/diff_median_pa_corr_plot.pdf",
  width = 10, height = 10
) # increased dimensions for better space

# reduce the margins (bottom, left, top, right)
par(mar = c(4, 6, 4, 4)) # adjusted margins to ensure labels fit
custom_corrplot_with_scatter(corr_matrix, df_pa_h2_diff_and_stats_all_traits,
  tl.col = "black",
  tl.srt = 45
)
# add a title to the plot
title(
  main = "Correlation plots between GBLUP median PA gain, median h2 deviations,
  Shannon entropies, and Fst values across all species and traits.",
  line = 1, cex.main = 1.2
)

# close the device to save the plot
dev.off()


# plot correlation matrices for df_pa_h2_diff_and_stats associated with specific
# species (spec_), here "Apple"
spec_ <- "Apple"

# compute the correlation matrix for the specific species
df_pa_h2_diff_and_stats_trait_ <- df_pa_h2_diff_and_stats[
  df_pa_h2_diff_and_stats$Species == spec_,
  vars_to_keep
]
corr_matrix_spec_ <- cor(df_pa_h2_diff_and_stats_trait_)

# create a pdf file with adjusted dimensions
pdf(paste0("../results/resdiff_median_pa_corr_plot_", spec_, ".pdf"),
  width = 10, height = 10
) # dimensions in inches

# plot the correlation matrix
custom_corrplot_with_scatter(corr_matrix_spec_,
  df_pa_h2_diff_and_stats_trait_,
  tl.col = "black",
  tl.srt = 45
)
# define title
title(
  main = paste0("Correlation plots between GBLUP median PA gain, median h2 deviations,
  Shannon entropies, and Fst values for ", tolower(spec_), " traits."),
  line = 1, cex.main = 1.2
)

# close the pdf device
dev.off()

# save reformatted pa and h2 results
fwrite(
  df_pa_median_iqr_all_species_traits,
  "../results/gblup_pa_median_iqr_all_species_traits_and_methods.csv"
)
fwrite(
  df_h2_result_all_species_traits,
  "../results/gblup_h2_median_iqr_all_species_traits_and_methods.csv"
)

# df_pa_median_iqr_all_species_traits <- as.data.frame(fread(
#   "../results/gblup_pa_median_iqr_all_species_traits_and_methods.csv"
# ))
# df_h2_result_all_species_traits <- as.data.frame(fread(
#   "../results/gblup_h2_median_iqr_all_species_traits_and_methods.csv"
# ))

# function to extract iqr values and compute their means
calculate_mean_parentheses <- function(column) {
  values <- gsub(".*\\((.*?)\\).*", "\\1", column)
  numeric_values <- as.numeric(values)
  mean(numeric_values, na.rm = TRUE)
}

# apply function to compute mean of iqr values
columns_to_calculate <- c(
  "GBLUP_median_PA_wiser",
  "GBLUP_median_PA_ls_means",
  "GBLUP_median_PA_blups"
)
mean_pa_iqr_results <- sapply(
  df_pa_median_iqr_all_species_traits[columns_to_calculate],
  calculate_mean_parentheses
)
print(mean_pa_iqr_results)

columns_to_calculate <- c(
  "GBLUP_median_h2_wiser",
  "GBLUP_median_h2_ls_means",
  "GBLUP_median_h2_blups"
)
mean_h2_iqr_results <- sapply(
  df_h2_result_all_species_traits[columns_to_calculate],
  calculate_mean_parentheses
)
print(mean_h2_iqr_results)
