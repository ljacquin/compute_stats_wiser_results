# function to compute entropy for a factor
compute_entropy_factor <- function(v) {
  freq_tab <- table(v)
  prob <- freq_tab / sum(freq_tab)
  return(-sum(prob * log2(prob)))
}

# function to compute joint entropy between two factors x and y
compute_joint_entropy_factors <- function(x, y) {
  # create the contingency table
  cont_table <- table(x, y)
  # compute relative frequencies
  p_ij <- as.vector(cont_table / sum(cont_table))
  # return joint entropy
  return(-sum(p_ij * log2(p_ij), na.rm = TRUE))
}

# function to determine the optimal k based on Calinsky-Habaraz index maximization
compute_ch_opt_k <- function(X, k_range = 2:30, n_init = 500) {
  # X : data frame of n x p data
  # k_range : vector of possible values for K (e.g., 2 to 10)
  # n_init : number of random starts for K-means (default is 25)

  # store the calinski-harabasz scores for each K
  ch_scores <- numeric(length(k_range))

  # loop through the different k values in the range
  for (k in k_range) {
    # perform k-means with k initial centers
    kmeans_result <- kmeans(X, centers = k, nstart = n_init)

    # compute the calinski-harabasz index for each k
    ch_scores[k - min(k_range) + 1] <- cluster.stats(
      dist(X), kmeans_result$cluster
    )$ch
  }

  # select the k with the highest ch score
  opt_k <- k_range[which.max(ch_scores)]

  # return the best k and the associated calinski-harabasz score
  return(list(opt_k = opt_k, ch_score = max(ch_scores)))
}

# function to determine the optimal k, in a paralllel fashion, based
# on Calinsky-Habaraz index maximization
compute_ch_opt_k_parallel <- function(X, k_range = 2:30, n_init = 500) {
  # X: data frame of n x p data
  # k_range: vector of possible values for k (e.g., 2 to 10)
  # n_init: number of random starts for k-means (default is 25)

  # set up parallel execution (use multisession or another backend)
  plan(multisession)

  # calculate calinski-harabasz scores for each k in parallel
  kmeans_list <- future_lapply(k_range, function(k) {
    # perform k-means with k initial centers
    kmeans_result <- kmeans(X, centers = k, nstart = n_init)

    # create output list
    out_list_ <- list(
      "cluster_idx" = kmeans_result$cluster,
      "ch" = cluster.stats(dist(X), kmeans_result$cluster)$ch
    )
    out_list_
  }, future.seed = TRUE) # ensure parallel-safe random numbers

  # get ch vector of values
  ch_vect <- sapply(kmeans_list, function(x) x$ch)

  # select the k with the highest calinski-harabasz score
  opt_k <- k_range[which.max(ch_vect)]

  # get the associated indices for the clustering
  opt_clust_vect <- kmeans_list[[which.max(ch_vect)]]$cluster_idx

  # reset to sequential execution
  plan(sequential)

  # return the best k and the associated calinski-harabasz score
  return(list(
    opt_ch = max(ch_vect),
    opt_k = opt_k,
    opt_clust_vect = opt_clust_vect,
    ch_vect = ch_vect
  ))
}

# function which removes monomorphic markers
remove_monomorphic_markers <- function(geno_df) {
  if ("Genotype" %in% colnames(geno_df)) {
    # remove genotype column
    geno_df_ <- geno_df[, -match("Genotype", colnames(geno_df))]
  } else {
    geno_df_ <- geno_df
  }
  # identify monomorphic markers
  monomorphic_markers <- apply(
    geno_df_,
    2, function(col) length(unique(col)) == 1
  )
  # get the names of the monomorphic markers
  monomorphic_marker_names <- colnames(geno_df_)[
    monomorphic_markers
  ]

  if (length(monomorphic_markers) > 0) {
    # filter the monomorphic markers
    geno_df_filtered <- geno_df_[, !monomorphic_markers]

    if ("Genotype" %in% colnames(geno_df)) {
      # add genotype column
      geno_df_filtered <- cbind(geno_df$Genotype, geno_df_filtered)
      colnames(geno_df_filtered)[1] <- "Genotype"
    }

    # return the filtered data frame and the list of monomorphic markers
    return(list(
      "filtered_df" = geno_df_filtered,
      "monomorphic_markers" = monomorphic_marker_names
    ))
  } else {
    # return the filtered data frame and the list of monomorphic markers
    return(list(
      "filtered_df" = geno_df,
      "monomorphic_markers" = NULL
    ))
  }
}


# custom corrplot function with ellipses and only lower triangle
custom_corrplot_with_scatter <- function(corr_matrix, data, ...) {
  # base corrplot with ellipses and only the lower triangle
  corrplot(corr_matrix,
    method = "ellipse", type = "lower",
    diag = FALSE, addCoef.col = "black", ...
  )

  # get the number of variables
  n <- ncol(corr_matrix)

  # overlay scatter plots in the lower triangle
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) { # only lower triangle
        # get the data for the current cell
        x_data <- data[[j]]
        y_data <- data[[i]]

        # normalize data to fit within the cell
        x_norm <- (x_data - min(x_data)) / (max(x_data) - min(x_data))
        y_norm <- (y_data - min(y_data)) / (max(y_data) - min(y_data))

        # scale normalized data to the cell size
        x_scaled <- 0.8 * x_norm + j - 0.4
        y_scaled <- 0.8 * y_norm + (n - i + 1) - 0.4

        # add scatter points
        points(x_scaled, y_scaled, pch = 16, col = "palegreen3", cex = 0.5)
      }
    }
  }
}

# functions which performs imputation using column means
impute_mean <- function(x) {
  # treat NaN as NA
  x[is.nan(x)] <- NA
  mean_value <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- mean_value
  return(x)
}

# function to select the 1000 snps with the highest variability
high_variance_snp <- function(df_genomic_species,
                              nb_snp_ = 1000) {
  # exclude the first column if it contains identifiers
  genomic_data <- df_genomic_species[, -1]

  # calculate the variance for each snp
  snp_var_ <- apply(genomic_data, 2, var)

  # identify the indices of the 1000 snps with the highest variance
  highest_var_snps_idx <- order(snp_var_, decreasing = TRUE)[1:nb_snp_]

  # select the columns corresponding to the top snps
  high_variance_snp_geno_df <- df_genomic_species[
    , c(1, highest_var_snps_idx + 1)
  ]

  return(high_variance_snp_geno_df)
}

# function to compute entropy based on SVD and select enough singular values
# to explain alpha% of the variance
compute_eigen_entropy <- function(matrix, alpha = 0.90) {
  # perform Singular Value Decomposition
  svd_result <- fast.svd(matrix)

  # extract the singular values (the diagonal elements of the sigma matrix)
  singular_values <- svd_result$d

  # compute the total variance (sum of squared singular values)
  total_variance <- sum(singular_values^2)

  # compute the cumulative sum of squared singular values
  cumulative_variance <- cumsum(singular_values^2) / total_variance

  # find the number of singular values needed to explain at least 'alpha' of the variance
  num_singular <- which(cumulative_variance >= alpha)[1]

  # select the top 'num_singular' singular values
  singular_values_top <- singular_values[1:num_singular]

  # normalize the selected singular values so that their sum equals 1
  singular_values_norm <- singular_values_top / sum(singular_values_top)

  # compute Shannon entropy using the normalized singular values
  entropy <- -sum(singular_values_norm * log2(singular_values_norm), na.rm = TRUE)

  return(entropy)
}

# function to filter snps by maf for a given threshold
filter_snp_by_maf <- function(df_, maf_threshold = 0.01) {
  # remove genotype id
  geno_df <- df_[, -1]

  # compute maf for each snp
  maf <- apply(geno_df, 2, function(snp) {
    allele_counts <- table(factor(snp, levels = 0:2)) # count alleles
    total_alleles <- sum(allele_counts) * 2 # total number of alleles

    # allele frequencies for 0 and 1
    freq_allele1 <- (allele_counts[1] * 2 + allele_counts[2]) / total_alleles
    freq_allele2 <- 1 - freq_allele1

    # return maf
    min(freq_allele1, freq_allele2)
  })

  # filter snps with maf > maf_threshold
  snps_to_keep <- which(maf > maf_threshold)

  # keep snps with maf > maf_threshold
  filtered_data <- df_[, c(1, snps_to_keep + 1)] # add id column

  return(filtered_data)
}

# function to filter snps, by correlation with others, for a given threshold
filter_snps_by_correlation <- function(geno_df, threshold = 0.8) {
  # compute the correlation matrix between snps
  # note that monomorphic markers should be removed for correlation calculations
  geno_df <- remove_monomorphic_markers(geno_df)$filtered_df[, -1]
  fbm_obj <- as_FBM(apply(geno_df, 2, as.numeric))
  cor_matrix <- big_cor(fbm_obj)
  cor_matrix <- cor_matrix[1:nrow(cor_matrix), 1:ncol(cor_matrix)]
  rownames(cor_matrix) <- colnames(geno_df)
  snps_to_keep <- c()
  # identify which snp has a correlation lower than threshold with other snps
  for (i in 1:nrow(cor_matrix)) {
    # extract correlation of snp i with others (excluding diagonal)
    correlations <- cor_matrix[i, -i]

    # verify if the correlation with all other snp is less than threshold
    if (all(abs(correlations) < threshold)) {
      # add snp to list if it is the case
      snps_to_keep <- c(snps_to_keep, rownames(cor_matrix)[i])
    }
  }
  return(snps_to_keep)
}

# function to compute Fst
calculate_fst <- function(sub_matrix) {
  # create a genind object
  genind_obj <- df2genind(
    apply(sub_matrix[, -1], 2, as.numeric),
    pop = as.factor(sub_matrix$pop),
    ncode = 1
  )
  # compute basic stats
  stats_obj_ <- basic.stats(genind_obj)
  # get Fst value
  fst_ <- as.numeric(stats_obj_$overall[match(
    "Fst",
    names(stats_obj_$overall)
  )])
  return(fst_)
}

# function to split a big SNP matrix into a list of smaller SNP sub-matrices
# of n x chunk_size
split_matrix <- function(df_, chunk_size = 1000) {
  # nb of columns
  p <- ncol(df_) - 1 # exclude pop column
  # compute number K of sub-matrices
  K <- p %/% chunk_size
  # create a list of sub-matrices
  sub_matrices <- lapply(1:K, function(i) {
    if (i < K) {
      df_subset <- df_[, c(1, ((i - 1) * chunk_size + 2):(i * chunk_size + 1))]
    } else {
      df_subset <- df_[, c(1, ((i - 1) * chunk_size + 2):ncol(df_))]
    }
    return(df_subset)
  })
  return(sub_matrices)
}

# function to split a big matrix into sampled sub-matrices without replacement
split_matrix_by_sampling <- function(df_, chunk_size = 1000) {
  # initialize the available indices (excluding the first column, which is fixed)
  all_indices <- 2:ncol(df_) # exclude the "pop" column (column 1)

  # list to store the sub-matrices
  sub_matrices <- list()

  # repeat until all indices are exhausted
  while (length(all_indices) > 0) {
    # if the number of remaining indices is less than the chunk_size, take all remaining indices
    if (length(all_indices) <= chunk_size) {
      sampled_indices <- all_indices
    } else {
      # otherwise, randomly sample chunk_size indices without replacement
      sampled_indices <- sample(all_indices, chunk_size, replace = FALSE)
    }

    # create the sub-matrix including the first column and the sampled indices
    df_subset <- df_[, c(1, sampled_indices)]

    # add the sub-matrix to the list
    sub_matrices <- append(sub_matrices, list(df_subset))

    # remove the sampled indices from the list of available indices
    all_indices <- setdiff(all_indices, sampled_indices)
  }

  # return the list of sub-matrices
  return(sub_matrices)
}

