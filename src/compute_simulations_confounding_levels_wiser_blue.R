# load libraries
library(mvtnorm)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(MASS)
library(data.table)
library(stringr)
library(whitening)
library(lme4)
library(tidyr)
# library(rstudioapi)
# setwd(dirname(getActiveDocumentContext()$path))
setwd("compute_stats_wiser_results/src/")
source("functions.R")

# initialize parameters
set.seed(123)
n_sim <- 100
species_ <- "Maize" # Rice, Maize, Apple, Pine
sig2_total <- 1
h2 <- 0.5

# GxE interaction variance as a factor of genetic variance
ge_sig2u_factor <- 0.5

# whitening method and default regularization parameter
whiten_method_ <- "ZCA-cor"
alpha_ <- 0.1

# optimize computation time
snp_sample_size_ <- 50e3
nrow_lim <- 1000

# tolerance for target Cramer's V
cramer_tolerance_ <- 0.02

# maximum number of attempts to obtain target Cramer's V
max_cramer_iter_ <- 5000

# result path
output_result_path <- "../results/wiser_blue_ea_gain/"

# load datasets
omic_df <- fread(
  paste0("../data/genomic_data_", species_, ".csv")
) |>
  as.data.frame()

colnames(omic_df)[1] <- "Genotype"
omic_df$Genotype <- as.character(omic_df$Genotype)

pheno_df <- fread(
  paste0("../data/phenotype_data_", species_, ".csv")
) |>
  as.data.frame()

pheno_df$Genotype <- as.character(pheno_df$Genotype)
pheno_df$Envir <- as.character(pheno_df$Envir)

# sample datasets if necessary to manage computation resources
if (nrow(pheno_df) > nrow_lim) {
  pheno_df <- pheno_df[
    sample(
      seq_len(nrow(pheno_df)),
      size = nrow_lim,
      replace = FALSE
    ),
  ]
}
if (species_ == "Apple" && ncol(omic_df) - 1 > snp_sample_size_) {
  idx_snp_sample_size_ <- sample(
    2:ncol(omic_df),
    size = snp_sample_size_,
    replace = FALSE
  )
  omic_df <- omic_df[
    , c(1, idx_snp_sample_size_)
  ]
}

# harmonize genomic and phenotypic datasets
sel_geno <- intersect(
  omic_df$Genotype,
  pheno_df$Genotype
)
omic_df <- omic_df[
  omic_df$Genotype %in% sel_geno,
]
pheno_df <- pheno_df[
  pheno_df$Genotype %in% sel_geno,
]

# define specified Cramer's V
vect_cramer_v <- seq(
  from = 0.1,
  to = 1,
  length.out = n_sim
)

# initialize vectors for results
vect_rho_u_v_hat <- rep(NA_real_, n_sim)
vect_rho_u_geno_hat_blue <- rep(NA_real_, n_sim)

vect_mse_u_v_hat <- rep(NA_real_, n_sim)
vect_mse_u_geno_hat_blue <- rep(NA_real_, n_sim)

vect_observed_cramer_v <- rep(NA_real_, n_sim)

# simulations
for (sim_ in seq_len(n_sim)) {
  cat(
    "\nSimulation",
    sim_,
    "/",
    n_sim,
    "- target Cramer's V =",
    round(vect_cramer_v[sim_], 3),
    "\n"
  )

  # resample phenotype data according to target Cramer's V
  resampled_ <- resample_to_cramer_v(
    pheno_df_ = pheno_df,
    target_v_ = vect_cramer_v[sim_],
    tolerance_ = cramer_tolerance_,
    max_iter_ = max_cramer_iter_
  )

  pheno_df_sampled <- resampled_$pheno_df
  vect_observed_cramer_v[sim_] <- resampled_$cramer_v
  cat(
    "Observed Cramer's V =",
    round(
      vect_observed_cramer_v[sim_],
      3
    ),
    "\n"
  )

  # harmonize omic and phenotype datasets
  sel_geno <- intersect(
    omic_df$Genotype,
    pheno_df_sampled$Genotype
  )
  omic_df_sampled <- omic_df[
    omic_df$Genotype %in% sel_geno,
  ]
  pheno_df_sampled <- pheno_df_sampled[
    pheno_df_sampled$Genotype %in% sel_geno,
  ]

  # define genotype order
  geno_list <- sort(
    unique(
      pheno_df_sampled$Genotype
    )
  )

  # reorder genomic data accordingly
  omic_df_sampled <- omic_df_sampled[
    match(
      geno_list,
      omic_df_sampled$Genotype
    ),
  ]

  # safety check
  stopifnot(
    identical(
      omic_df_sampled$Genotype,
      geno_list
    )
  )

  # compute kinship matrix
  snp_mat <- apply(
    omic_df_sampled[, -1, drop = FALSE], 2, as.numeric
  )

  if (is.null(dim(snp_mat))) {
    snp_mat <- matrix(
      snp_mat,
      ncol = 1
    )
  }

  K <- snp_mat %*% t(snp_mat)
  K <- K / mean(diag(K))
  rownames(K) <- geno_list
  colnames(K) <- geno_list

  # compute incidence matrices
  env_list <- sort(
    unique(
      pheno_df_sampled$Envir
    )
  )
  n_geno_ <- length(geno_list)
  n_env_ <- length(env_list)
  len_y <- nrow(pheno_df_sampled)
  Z <- matrix(0,
    nrow = len_y,
    ncol = n_geno_
  )
  X_env <- matrix(0,
    nrow = len_y,
    ncol = n_env_
  )
  colnames(Z) <- geno_list
  colnames(X_env) <- env_list

  for (i in seq_len(len_y)) {
    Z[
      i,
      pheno_df_sampled$Genotype[i]
    ] <- 1

    X_env[
      i,
      pheno_df_sampled$Envir[i]
    ] <- 1
  }

  # remove one environment column because intercept +
  # all environment dummy variables are linearly dependent
  if (n_env_ > 1) {
    X <- cbind(
      Intercept = 1,
      X_env[, -1, drop = FALSE]
    )
  } else {
    X <- matrix(1,
      nrow = len_y,
      ncol = 1
    )
    colnames(X) <- "Intercept"
  }

  # compute variance components
  sig2u <- as.numeric(
    h2 * sig2_total
  )
  sig2ge <- as.numeric(
    ge_sig2u_factor * sig2u
  )
  sig2e <- as.numeric(
    sig2_total - sig2u - sig2ge
  )

  if (sig2e < 0) {
    message(
      "Residual variance is negative. ",
      "Check h2 and ge_sig2u_factor."
    )
    next
  }

  # compute G x E interaction incidence matrix
  Z_ge <- matrix(0,
    nrow = len_y,
    ncol = n_geno_ * n_env_
  )
  colnames(Z_ge) <- paste0(
    "GxE_",
    rep(
      geno_list,
      each = n_env_
    ), "_",
    rep(
      env_list,
      times = n_geno_
    )
  )

  for (i in seq_len(len_y)) {
    g <- pheno_df_sampled$Genotype[i]
    e <- pheno_df_sampled$Envir[i]
    Z_ge[
      i,
      paste0(
        "GxE_", g, "_", e
      )
    ] <- 1
  }

  # compute Whitening matrix
  Sigma_u <- as.numeric(sig2u) * Z %*% K %*% t(Z)

  if (!is.positive.definite(Sigma_u)) {
    Sigma_u <- frobenius_norm_regularization(
      Sigma_u,
      alpha_ = alpha_
    )
  }

  W <- whiteningMatrix(
    Sigma_u,
    method = whiten_method_
  )
  X_tilde <- W %*% X

  # simulate fixed effects
  mu_ <- runif(1, 0, 1)

  # one reference environment is absorbed into intercept
  if (ncol(X) > 1) {
    env_eff_ <- rnorm(
      ncol(X) - 1, 0, 3
    )
  } else {
    env_eff_ <- numeric(0)
  }
  beta_ <- c(
    mu_,
    env_eff_
  )

  # simulate residuals, genetic and GxE effects
  epsilon <- rnorm(
    len_y,
    0,
    sqrt(sig2e)
  )

  u <- drop(
    mvrnorm(1,
      mu = rep(0, n_geno_),
      Sigma = sig2u * K
    )
  )

  ge_ <- rnorm(
    n_geno_ * n_env_, 0,
    sqrt(sig2ge)
  )

  Y <- drop(
    X %*% beta_ +
      Z %*% u +
      Z_ge %*% ge_ +
      epsilon
  )

  pheno_df_sampled$trait <- Y

  geno_u_df <- data.frame(
    Genotype = geno_list,
    u = u - mean(u)
  )

  # wiser
  beta_hat <-
    ginv(
      t(X_tilde) %*% X_tilde
    ) %*%
    t(X_tilde) %*%
    Y


  v_hat <- ginv(t(Z) %*% Z) %*%
    t(Z) %*% (
      Y - X_tilde %*% beta_hat
    )
  df_v_hat <- data.frame(
    Genotype = geno_list,
    v_hat = drop(v_hat)
  )

  # wiser without whitening classical linear model
  pheno_df_sampled$Genotype <- factor(
    pheno_df_sampled$Genotype,
    levels = geno_list
  )
  pheno_df_sampled$Envir <- factor(
    pheno_df_sampled$Envir
  )
  lm_ <- lm(
    trait ~ Genotype + Envir,
    data = pheno_df_sampled
  )

  # obtain genotype blues including reference genotype
  beta_lm_ <- coef(lm_)
  intercept_ <- beta_lm_["(Intercept)"]
  geno_hat_blue <- rep(
    intercept_,
    n_geno_
  )
  names(geno_hat_blue) <- geno_list

  if (n_geno_ > 1) {
    for (g_ in geno_list[-1]) {
      coef_name_ <- paste0(
        "Genotype",
        g_
      )
      if (coef_name_ %in% names(beta_lm_) &&
        !is.na(beta_lm_[coef_name_])
      ) {
        geno_hat_blue[g_] <-
          intercept_ +
          beta_lm_[coef_name_]
      }
    }
  }
  geno_hat_df <- data.frame(
    Genotype = geno_list,
    geno_hat_blue = as.numeric(
      geno_hat_blue
    )
  )

  # center genotype estimates because genetic effects u
  # were simulated with expectation zero
  geno_hat_df$geno_hat_blue <-
    geno_hat_df$geno_hat_blue -
    mean(
      geno_hat_df$geno_hat_blue,
      na.rm = TRUE
    )

  # center wiser estimates for the same reason
  df_v_hat$v_hat <-
    df_v_hat$v_hat -
    mean(
      df_v_hat$v_hat,
      na.rm = TRUE
    )

  # merge data frames
  df_merge_wiser <- merge(
    geno_u_df,
    df_v_hat,
    by = "Genotype"
  )

  df_merge_blue <- merge(
    geno_u_df,
    geno_hat_df,
    by = "Genotype"
  )

  # get results
  vect_rho_u_v_hat[sim_] <- cor(
    df_merge_wiser$u,
    df_merge_wiser$v_hat,
    use = "complete.obs"
  )

  vect_mse_u_v_hat[sim_] <- mean(
    (df_merge_wiser$u - df_merge_wiser$v_hat)^2,
    na.rm = TRUE
  )

  vect_rho_u_geno_hat_blue[sim_] <- cor(
    df_merge_blue$u,
    df_merge_blue$geno_hat_blue,
    use = "complete.obs"
  )

  vect_mse_u_geno_hat_blue[sim_] <- mean(
    (df_merge_blue$u -
      df_merge_blue$geno_hat_blue
    )^2,
    na.rm = TRUE
  )
}

# create results data frame
results_df <- data.frame(
  Simulation = seq_len(n_sim),
  Target_Cramer_V = vect_cramer_v,
  Observed_Cramer_V = vect_observed_cramer_v,
  WISER_Correlation = vect_rho_u_v_hat,
  BLUE_Correlation = vect_rho_u_geno_hat_blue,
  WISER_MSE = vect_mse_u_v_hat,
  BLUE_MSE = vect_mse_u_geno_hat_blue,
  EA_gain = vect_rho_u_v_hat - vect_rho_u_geno_hat_blue
)

fwrite(
  results_df,
  file = paste0(
    output_result_path,
    species_,
    "_sim_results_cramer_v_h2_",
    h2,
    ".csv"
  )
)

# reshape correlation results
df_rho <- data.frame(
  Cramer_V = rep(vect_observed_cramer_v, 2),
  Method = rep(
    c(
      "WISER",
      "BLUE"
    ),
    each = n_sim
  ),
  Correlation = c(
    vect_rho_u_v_hat,
    vect_rho_u_geno_hat_blue
  )
)
df_rho <- df_rho |>
  drop_na()

# plot estimation accuracy versus cramer v
title_rho_ <- paste0(
  "Estimation accuracy (EA) between simulated and estimated genetic values ",
  "for v_hat (WISER) and u_hat (BLUE)\n",
  "for a simulated trait based on ",
  tolower(species_),
  " genomic and experimental design data (h2 = ",
  h2,
  "), 
  across different levels of genotype-environment association (i.e., confounding)"
)

ggplot_rho <- ggplot(
  df_rho,
  aes(
    x = Cramer_V,
    y = Correlation,
    color = Method
  )
) +
  geom_point(
    alpha = 0.6,
    size = 2
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_family = "Calibri") +
  labs(
    title = title_rho_,
    x = "Cramer's V",
    y = "Estimation accuracy (EA)",
    color = NULL
  ) +
  coord_cartesian(ylim = c(
    min(
      df_rho$Correlation,
      0,
      na.rm = TRUE
    ),
    1
  )) +
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 1
    ),
    axis.ticks = element_line(
      linewidth = 1,
      color = "black"
    ),
    axis.text = element_text(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  font("xylab", size = 15) +
  font("xy.text", size = 15) +
  font("legend.text", size = 15)


ggsave(
  filename = paste0(
    output_result_path, "sim_rho_cramer_v_",
    species_, "_h2_", h2, ".png"
  ),
  plot = ggplot_rho, width = 12, height = 8, dpi = 300
)
ggplot_rho

# reshape mse results
df_mse <- data.frame(
  Cramer_V = rep(
    vect_observed_cramer_v,
    2
  ),
  Method = rep(
    c(
      "WISER",
      "BLUE"
    ),
    each = n_sim
  ),
  MSE = c(
    vect_mse_u_v_hat,
    vect_mse_u_geno_hat_blue
  )
)
df_mse <- df_mse |>
  drop_na()

# plot mse versus cramer v
title_mse_ <- paste0(
  "Mean squared error (MSE) between simulated and estimated genetic values ",
  "for v_hat (WISER) and u_hat (BLUE)\n",
  "for a simulated trait based on ",
  tolower(species_),
  " genomic and experimental design data (h2 = ",
  h2,
  "), 
  across different levels of genotype-environment association (i.e., confounding)"
)

ggplot_mse <- ggplot(
  df_mse,
  aes(
    x = Cramer_V,
    y = MSE,
    color = Method
  )
) +
  geom_point(
    alpha = 0.6,
    size = 2
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_family = "Calibri") +
  labs(
    title = title_mse_, x = "Cramer's V",
    y = "Mean squared error (MSE)",
    color = NULL
  ) +
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 1
    ), axis.ticks = element_line(
      linewidth = 1,
      color = "black"
    ),
    axis.text = element_text(color = "black"),
    axis.ticks.length = unit(
      0.2,
      "cm"
    )
  ) +
  font("xylab", size = 15) +
  font("xy.text", size = 15) +
  font("legend.text",
    size = 15
  )

ggsave(
  filename = paste0(
    output_result_path,
    "sim_mse_cramer_v_",
    species_,
    "_h2_",
    h2,
    ".png"
  ),
  plot = ggplot_mse, width = 12, height = 8, dpi = 300
)

# summary statistics
medians_df <- data.frame(
  species = species_,
  statistic = c(
    "Correlation",
    "MSE"
  ),
  WISER = c(
    format_median_iqr(
      df_rho$Correlation[
        df_rho$Method == "WISER"
      ]
    ),
    format_median_iqr(
      df_mse$MSE[
        df_mse$Method == "WISER"
      ]
    )
  ),
  BLUE = c(
    format_median_iqr(
      df_rho$Correlation[
        df_rho$Method == "BLUE"
      ]
    ),
    format_median_iqr(
      df_mse$MSE[
        df_mse$Method == "BLUE"
      ]
    )
  )
)

# create EA gain data frame
df_ea_gain <- data.frame(
  Cramer_V = vect_observed_cramer_v,
  EA_gain = vect_rho_u_v_hat - vect_rho_u_geno_hat_blue
)

df_ea_gain <- df_ea_gain |>
  drop_na()

# plot EA gain versus Cramer's V
title_ea_gain_ <- paste0(
  "Estimation accuracy gain of WISER over BLUE ",
  "(EA gain = EA[WISER] - EA[BLUE])\n",
  "for a simulated trait based on ",
  tolower(species_),
  " genomic and experimental design data (h2 = ",
  h2,
  "), 
  across different levels of genotype-environment association (i.e., confounding)"
)

ggplot_ea_gain <- ggplot(
  df_ea_gain,
  aes(
    x = Cramer_V,
    y = EA_gain
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_point(
    alpha = 0.6,
    size = 2
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  theme_minimal(base_family = "Calibri") +
  labs(
    title = title_ea_gain_,
    x = "Cramer's V",
    y = "EA gain (WISER - BLUE)"
  ) +
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 1
    ),
    axis.ticks = element_line(
      linewidth = 1,
      color = "black"
    ),
    axis.text = element_text(
      color = "black"
    ),
    axis.ticks.length = unit(
      0.2,
      "cm"
    )
  ) +
  font("xylab", size = 15) +
  font("xy.text", size = 15)

ggsave(
  filename = paste0(
    output_result_path,
    "sim_ea_gain_cramer_v_",
    species_,
    "_h2_",
    h2,
    ".png"
  ),
  plot = ggplot_ea_gain,
  width = 12,
  height = 8,
  dpi = 300
)

ggplot_ea_gain

# write summary statistics
fwrite(
  x = medians_df,
  file = paste0(
    output_result_path,
    species_,
    "_sim_stats_medians_h2_",
    h2,
    ".csv"
  )
)
medians_df
cor(results_df$Observed_Cramer_V, results_df$EA_gain)