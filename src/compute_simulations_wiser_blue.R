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
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
source("../src/functions.R")

# initialize parameters
set.seed(123)
n_sim <- 100
species_ <- "Apple"
sig2_total <- 1        # set total phenotypic variance to 1
h2 <- 0.5              # set heritability
ge_sig2u_factor <- 0.5 # set gxe interaction variance as a factor of genetic variance
whiten_method_ <- "ZCA-cor"
alpha_ <- 0.1
snp_sample_size_ <- 50e3
nrow_lim <- 1000       # optimize computation time

# define result path
output_result_path <- "../results/"

# load datasets
omic_df <- fread(paste0("../data/genomic_data_", species_, ".csv")) |>
  as.data.frame()
colnames(omic_df)[1] <- "Genotype"
omic_df$Genotype <- as.character(omic_df$Genotype)

pheno_df <- fread(paste0("../data/phenotype_data_", species_, ".csv")) |>
  as.data.frame()
pheno_df$Genotype <- as.character(pheno_df$Genotype)

# sample datasets to manage computation ressources
if (nrow(pheno_df) > nrow_lim) {
  pheno_df <- pheno_df[
    sample(1:nrow(pheno_df), size = nrow_lim, replace = F),
  ]
}
if (species_ == "Apple") {
  idx_snp_sample_size_ <- sample(2:ncol(omic_df),
    size = snp_sample_size_, replace = F
  )
  omic_df <- omic_df[, c(1, idx_snp_sample_size_)]
}

# harmonize omic and phenotype datasets
sel_geno <- intersect(omic_df$Genotype, pheno_df$Genotype)
omic_df <- omic_df[omic_df$Genotype %in% sel_geno, ]
pheno_df <- pheno_df[pheno_df$Genotype %in% sel_geno, ]

# compute kinship matrix
snp_mat <- apply(omic_df[, -1], 2, as.numeric)
K <- snp_mat %*% t(snp_mat)
K <- K / (sum(diag(K)) / nrow(K))

# compute incidence matrices
geno_list <- sort(unique(pheno_df$Genotype))
env_list <- sort(unique(pheno_df$Envir))
n_geno_ <- length(geno_list)
n_env_ <- length(env_list)
len_y <- nrow(pheno_df)

Z <- matrix(0, len_y, n_geno_)
X <- matrix(0, len_y, n_env_)
colnames(Z) <- geno_list
colnames(X) <- env_list

for (i in 1:len_y) {
  Z[i, pheno_df$Genotype[i]] <- 1
  X[i, pheno_df$Envir[i]] <- 1
}
X <- cbind(Intercept = 1, X)

# compute variance components
sig2u <- as.numeric(h2 * sig2_total)
sig2ge <- as.numeric(ge_sig2u_factor * sig2u)
sig2e <- as.numeric(sig2_total - sig2u - sig2ge)

# compute g x e interaction incidence matrix
Z_ge <- matrix(0, len_y, n_geno_ * n_env_)
colnames(Z_ge) <- paste0(
  "GxE_", rep(geno_list, each = n_env_),
  "_", rep(env_list, times = n_geno_)
)
for (i in 1:len_y) {
  g <- pheno_df$Genotype[i]
  e <- pheno_df$Envir[i]
  Z_ge[i, paste0("GxE_", g, "_", e)] <- 1
}

# compute Whitening matrix
Sigma_u <- as.numeric(sig2u) * Z %*% K %*% t(Z)
if (!is.positive.definite(Sigma_u)) {
  Sigma_u <- frobenius_norm_regularization(
    Sigma_u,
    alpha_ = alpha_
  )
}
W <- whiteningMatrix(Sigma_u, method = whiten_method_)
X_tilde <- W %*% X

# initialize vectors for results
vect_rho_u_v_hat <- vect_rho_u_geno_hat_blue <- rep(0, n_sim)
vect_mse_u_v_hat <- vect_mse_u_geno_hat_blue <- rep(0, n_sim)

for (sim_ in 1:n_sim) {
  cat("Simulation ", sim_, "\n")

  # simulate fixed effects parameters
  mu_ <- runif(1, 0, 1)
  env_eff_ <- rnorm(n_env_, 0, 3)
  beta_ <- c(mu_, env_eff_)

  # simulate residuals, genetic and g x e effects
  epsilon <- rnorm(len_y, 0, sqrt(sig2e))
  u <- drop(mvrnorm(1, mu = rep(0, n_geno_), Sigma = sig2u * K))
  ge_ <- rnorm(n_geno_ * n_env_, 0, sqrt(sig2ge))

  Y <- X %*% beta_ + Z %*% u + Z_ge %*% ge_ + epsilon
  pheno_df$trait <- Y
  geno_u_df <- data.frame(Genotype = geno_list, u = u)

  # WISER
  beta_hat <- ginv(t(X_tilde) %*% X_tilde) %*% t(X_tilde) %*% Y
  v_hat <- ginv(t(Z) %*% Z) %*% t(Z) %*% (Y - X_tilde %*% beta_hat)
  df_v_hat <- data.frame(Genotype = geno_list, v_hat = v_hat)

  # WISER without whitening, i.e. classical linear model with Envir and
  # Genotype as fixed effects
  lm_ <- lm(trait ~ Genotype + Envir, data = pheno_df)
  geno_hat <- coef(lm_)[grep("Genotype", names(coef(lm_)))]
  geno_hat_df <- data.frame(
    Genotype = gsub("Genotype", "", names(geno_hat)),
    geno_hat_blue = as.numeric(geno_hat)
  )

  # merge data frames for integrity of analyzes
  df_merge_wiser <- merge(geno_u_df, df_v_hat, by = "Genotype")
  df_merge_blue <- merge(geno_u_df, geno_hat_df, by = "Genotype")

  # get results
  vect_rho_u_v_hat[sim_] <- cor(df_merge_wiser$u, df_merge_wiser$v_hat)
  vect_mse_u_v_hat[sim_] <- mean((df_merge_wiser$u - df_merge_wiser$v_hat)^2)

  vect_rho_u_geno_hat_blue[sim_] <- cor(df_merge_blue$u, df_merge_blue$geno_hat_blue)
  vect_mse_u_geno_hat_blue[sim_] <- mean((df_merge_blue$u - df_merge_blue$geno_hat_blue)^2)
}

# create data frames for graphical results
df_rho <- data.frame(
  Method = rep(c("WISER", "BLUE"), each = n_sim),
  Correlation = c(vect_rho_u_v_hat, vect_rho_u_geno_hat_blue)
)
df_mse <- data.frame(
  Method = rep(c("WISER", "BLUE"), each = n_sim),
  MSE = c(vect_mse_u_v_hat, vect_mse_u_geno_hat_blue)
)

# remove NA values
df_rho <- df_rho %>% drop_na()
title_ <- paste0(
  "Predictive ability (PA) distributions between simulated and estimated genetic values,
  for v_hat (WISER) and u_hat (BLUE), for simulated trait based on ", tolower(species_),
  " genomic and experimental design data (h2 = ",
  h2, ")"
)

# make ggplot
ggplot_rho <- ggplot(df_rho, aes(x = Method, y = Correlation, fill = Method)) +
  geom_violin(alpha = 0.2, trim = FALSE) + # Violin plot
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) + # Boxplot
  theme_minimal() +
  labs(
    title =
    ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_family = "Calibri") +
  ylab("Predictive ability (PA)") +
  xlab("Genetic value estimation method") +
  rremove("legend.title") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
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
  coord_cartesian(ylim = c(min(df_rho$Correlation, 0), 1)) +
  ggtitle(title_)

# save the plot as png
ggsave(
  filename = paste0(
    output_result_path, "sim_rho_", species_, "_h2_",
    h2, ".png"
  ),
  plot = ggplot_rho,
  width = 18,
  height = 8,
  dpi = 300
)
ggplot_rho

# remove NA values
df_mse <- df_mse %>% drop_na()

title_ <- paste0(
  "Mean squared error (MSE) between simulated and estimated genetic values,
  for v_hat (WISER) and u_hat (BLUE), for simulated trait based on ", tolower(species_),
  " genomic and experimental design data (h2 = ",
  h2, ")"
)

# make ggplot
ggplot_mse <- ggplot(df_mse, aes(x = Method, y = MSE, fill = Method)) +
  geom_violin(alpha = 0.2, trim = FALSE) + # Violin plot
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) + # Boxplot
  theme_minimal() +
  labs(
    title =
    ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_family = "Calibri") +
  ylab("Mean squared error (MSE)") +
  xlab("Genetic value estimation method") +
  rremove("legend.title") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
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
    output_result_path, "sim_mse_", species_, "_h2_",
    h2, ".png"
  ),
  plot = ggplot_mse,
  width = 18,
  height = 8,
  dpi = 300
)
ggplot_rho

# compute medians for MSE and correlations for each genetic value estimation method
medians_df <- data.frame(
  species = species_,
  statistic = c("Correlation", "MSE"),
  BLUE = c(
    format_median_iqr(df_rho$Correlation[df_rho$Method == "BLUE"]),
    format_median_iqr(df_mse$MSE[df_mse$Method == "BLUE"])
  ),
  WISER = c(
    format_median_iqr(df_rho$Correlation[df_rho$Method == "WISER"]),
    format_median_iqr(df_mse$MSE[df_mse$Method == "WISER"])
  )
)

# write data frame of statistics medians
fwrite(x = medians_df, file = paste0(
  output_result_path, species_,
  "_sim_stats_medians_h2_", h2,
  ".csv"
))
medians_df
