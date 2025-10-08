Sys.setenv(OPENBLAS_NUM_THREADS = "1")
library(data.table)
library(dplyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(patchwork)
library(stringr)
#library(ggtext)
library(tidyr)
rm(list = ls())

trait_name1 <- "T2D"
trait_name2 = "Height"

data_dir <- "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/"
annot_file <- paste0(data_dir, "cell_type_annot_human_total162.filtered.with_rest.annot")
annot_df <- as.data.frame(fread(annot_file, header = T))
sum(colSums(annot_df[,-1]))
nrow(annot_df)
annotationName <- colnames(annot_df)[-1]
res_dir <- "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/analysis/real_data/cell_type_human_total162_wRest/MTSBayesCC_3K_corrR_tuned_v3_estSigma/"
nCat <- ncol(annot_df) - 1
print(paste("Number of annotations:", nCat))

annot_weight = 1 / rowSums(annot_df[,-1])
X <- as.matrix(annot_df[, -1])   # 0/1 matrix
# weighted sum per annotation = sum_i (X[i,j] * w[i])
weight_sum <- as.numeric(crossprod(X, annot_weight))   # same as t(X) %*% w
names(weight_sum) <- annotationName

total_seeds_T2D_FG = c(631,910,2374,11243,11578,15679,25727,31409,31642,35563,52643,57449,57515,57623,59688,
62081,66074,66609,67309,67533,69516,71532,78866,80640,85209,88270,89157,93198,94188,94238)
seed_groups <- list(
  Group1 = total_seeds_T2D_FG[1:10],
  Group2 = total_seeds_T2D_FG[11:20],
  Group3 = total_seeds_T2D_FG[21:30]
)
seed_groups <- list(
  Group1 = total_seeds_T2D_FG[1:15],
  Group2 = total_seeds_T2D_FG[16:30]
)
total_seeds_T2D_FG = c(30:39)

total_seeds <- as.numeric(fread(sprintf("%s/%s_%s/seeds_total162_wRest.txt", data_dir, trait_name1, trait_name2), header = FALSE)[1,])
seed_groups <- list(
  Group1 = total_seeds[1:10]
)

# ---- Function: Load and extract off-diagonal covariances from MCMC sample ----
read_covariance_samples <- function(mcmc_file, n_annot, iter_step = NULL) {
    lines <- readLines(mcmc_file)
    matrix_values <- do.call(rbind, strsplit(lines, ","))
    matrix_values <- apply(matrix_values, 2, as.numeric)

    n_matrices <- nrow(matrix_values) / 2
    stopifnot(n_matrices %% n_annot == 0)

    n_samples <- n_matrices / n_annot
    result <- matrix(NA, nrow = n_samples, ncol = n_annot + 1)

    for (i in 1:n_samples) {
        for (j in 1:n_annot) {
            mat_idx <- (i - 1) * n_annot + j
            row_start <- (mat_idx - 1) * 2 + 1
            mat <- matrix(matrix_values[row_start:(row_start + 1), ], nrow = 2, byrow = TRUE)
            result[i, j + 1] <- mat[1, 2]
        }
        result[i, 1] <- i * iter_step
    }

    df <- as.data.frame(result)
    colnames(df) <- c("Iteration", paste0("Cov_Annot", 1:n_annot))
    return(df)
}

read_h2_samples <- function(mcmc_file, n_annot, iter_step = NULL, trait = NULL) {
    lines <- readLines(mcmc_file)
    matrix_values <- do.call(rbind, strsplit(lines, ","))
    matrix_values <- apply(matrix_values, 2, as.numeric)

    n_matrices <- nrow(matrix_values) / 2
    stopifnot(n_matrices %% n_annot == 0)

    n_samples <- n_matrices / n_annot
    result <- matrix(NA, nrow = n_samples, ncol = n_annot + 1)

    for (i in 1:n_samples) {
        for (j in 1:n_annot) {
            mat_idx <- (i - 1) * n_annot + j
            row_start <- (mat_idx - 1) * 2 + 1
            mat <- matrix(matrix_values[row_start:(row_start + 1), ], nrow = 2, byrow = TRUE)
            result[i, j + 1] <- mat[trait, trait]
        }
        result[i, 1] <- i * iter_step
    }

    df <- as.data.frame(result)
    colnames(df) <- c("Iteration", paste0("Cov_Annot", 1:n_annot))
    return(df)
}

read_correlation_samples <- function(mcmc_file, n_annot, iter_step = NULL) {
  lines <- readLines(mcmc_file)
  matrix_values <- do.call(rbind, strsplit(lines, ","))
  matrix_values <- apply(matrix_values, 2, as.numeric)

  n_matrices <- nrow(matrix_values) / 2
  stopifnot(n_matrices %% n_annot == 0)

  n_samples <- n_matrices / n_annot
  result <- matrix(NA, nrow = n_samples, ncol = n_annot + 1)

  for (i in 1:n_samples) {
    for (j in 1:n_annot) {
      mat_idx <- (i - 1) * n_annot + j
      row_start <- (mat_idx - 1) * 2 + 1
      mat <- matrix(matrix_values[row_start:(row_start + 1), ], nrow = 2, byrow = TRUE)

      # extract gcov and diagonal variances
      gcov <- mat[1, 2]
      varg1 <- mat[1, 1]
      varg2 <- mat[2, 2]

      # compute gcor (with safety for division by zero)
      if (varg1 > 0 && varg2 > 0) {
        gcor <- gcov / sqrt(varg1 * varg2)
      } else {
        gcor <- NA # or 0, depending on how you want to handle this
      }

      result[i, j + 1] <- gcor
    }
    result[i, 1] <- i * iter_step
  }

  df <- as.data.frame(result)
  colnames(df) <- c("Iteration", paste0("Gcor_Annot", 1:n_annot))
  return(df)
}

read_total_gcov_gcor_h2_samples <- function(mcmc_file, iter_step = 50) {
  lines <- readLines(mcmc_file)
  matrix_values <- do.call(rbind, strsplit(lines, ","))
  matrix_values <- apply(matrix_values, 2, as.numeric)

  stopifnot(nrow(matrix_values) %% 2 == 0)

  n_samples <- nrow(matrix_values) / 2
  result <- matrix(NA, nrow = n_samples, ncol = 5) # Now includes iteration

  for (i in 1:n_samples) {
    row_start <- (i - 1) * 2 + 1
    mat <- matrix(matrix_values[row_start:(row_start + 1), ], nrow = 2, byrow = TRUE)

    varg1 <- mat[1, 1]
    varg2 <- mat[2, 2]
    gcov <- mat[1, 2]
    gcor <- if (varg1 > 0 && varg2 > 0) gcov / sqrt(varg1 * varg2) else NA

    result[i, ] <- c(i * iter_step, gcov, varg1, varg2, gcor)
  }

  df <- as.data.frame(result)
  colnames(df) <- c("Iteration", "coh2", "h21", "h22", "gcor")
  return(df)
}

run_enrichment_by_chain_length <- function(df_samples, df_sample_total, max_iter, min_iter, annotationName, annot_df, out_dir, para, weight_sum, absGtotal = FALSE) {
    mcmc_Gcov_c <- df_samples %>% filter(Iteration <= max_iter) %>% 
        filter(Iteration > min_iter) %>%
        select(-Iteration) # Remove Iteration column
    mcmc_Gcov_total <- df_sample_total %>% filter(Iteration <= max_iter) %>%
        filter(Iteration > min_iter) %>%
        pull(para) # Extract the relevant column as a vector
    if (absGtotal) {
        mcmc_Gcov_total <- abs(mcmc_Gcov_total)
    }
    calculateConditionalEnrichments2(mcmc_Gcov_total, mcmc_Gcov_c, annotationName, annot_df, out_dir, para, max_iter, weight_sum)
}

calculateGcovEnrichment <- function(annoIndex, mcmcGcovEnr_c, NtotalEffs, annoSize, annotationName) {
    annoti <- annotationName[annoIndex]
    message("Calculating enrichment for ", annoti)
    mcmcGcov_enrichment <- as.numeric(mcmcGcovEnr_c[[annoIndex]]) * (NtotalEffs / annoSize[annoIndex])
    return(list(mean = mean(mcmcGcov_enrichment), sd = sd(mcmcGcov_enrichment)))
}

# Helper function to calculate proportion of posterior probabilities
calculateGcovEnrichmentPP <- function(annoIndex, mcmcGcovEnr_c, NtotalEffs, annoSize, annotationName) {
    annoti <- annotationName[annoIndex]
    message("Calculating enrichment for ", annoti)
    mcmcGcov_enrichment <- as.numeric(mcmcGcovEnr_c[[annoIndex]]) * (NtotalEffs / annoSize[annoIndex])
    return(sum(mcmcGcov_enrichment > 1 | mcmcGcov_enrichment < -1 ) / length(mcmcGcov_enrichment))
}

calculateGcovDepletionPP <- function(annoIndex, mcmcGcovEnr_c, NtotalEffs, annoSize, annotationName) {
  annoti <- annotationName[annoIndex]
  message("Calculating depletion PP for ", annoti)
  # Scale the MCMC samples
  mcmcGcov_enrichment <- as.numeric(mcmcGcovEnr_c[[annoIndex]]) * (NtotalEffs / annoSize[annoIndex])
  return(sum(mcmcGcov_enrichment > -1 & mcmcGcov_enrichment < 1) / length(mcmcGcov_enrichment))

}

calculateConditionalEnrichments2 <- function(mcmcGcov_total, mcmcGcov_c, annotationName, annot_df, res_dir, para, max_iter, weight_sum) {
    # Check if mcmcGcov_c has correct column names
    colnames(mcmcGcov_c) <- annotationName

    # Conditional genetic covariance
    conditional_gcov <- colMeans(mcmcGcov_c)
    conditional_gcov_df <- data.frame(annotationName, conditional_gcov)
    conditional_gcov_df <- conditional_gcov_df[order(conditional_gcov_df$conditional_gcov, decreasing = TRUE), ]
    write.csv(conditional_gcov_df, file = file.path(res_dir, sprintf("conditional_%s_%d.csv", para, max_iter)), row.names = FALSE, quote = FALSE)

    # Calculate enrichment
    mcmcGcovEnr_c <- as.data.table(mcmcGcov_c) / mcmcGcov_total

    # Annotation sizes and names
    annoSize <- weight_sum
    NtotalEffs = nrow(annot_df)

    # Initialize results data frame
    Genr_results <- data.frame(
        Annotation = annotationName,
        Mean = numeric(length(annotationName)),
        SD = numeric(length(annotationName))
    )

    # Calculate enrichment for each annotation and store results
    for (i in seq_along(annotationName)) {
        enrichment <- calculateGcovEnrichment(i, mcmcGcovEnr_c, NtotalEffs, annoSize, annotationName)
        Genr_results$Mean[i] <- enrichment$mean
        Genr_results$SD[i] <- enrichment$sd
    }

    # Add Enrichment Posterior Probability (EnrPP)
    Genr_results$EnrPP <- sapply(1:length(annotationName), function(i) {
        calculateGcovEnrichmentPP(i, mcmcGcovEnr_c, NtotalEffs, annoSize, annotationName)
    })
    Genr_results$DeplPP <- sapply(1:length(annotationName), function(i) {
      calculateGcovDepletionPP(i, mcmcGcovEnr_c, NtotalEffs, annoSize, annotationName)
    })

    # Confidence intervals
    Genr_results$lb <- Genr_results$Mean - Genr_results$SD
    Genr_results$ub <- Genr_results$Mean + Genr_results$SD

    # Classify results based on Mean and confidence intervals
    Genr_results <- Genr_results %>%
        mutate(Significance = case_when(
            Mean > 0 & EnrPP > 0.90 ~ "Enrichment",
            Mean > 0 & DeplPP > 0.90 ~ "Depletion",
            Mean < 0 & EnrPP > 0.90 ~ "Enrichment",
            Mean < 0 & DeplPP > 0.90 ~ "Depletion",
            TRUE ~ "Non-significant"
        ))

    # Save the results
    write.csv(Genr_results, file = file.path(res_dir, sprintf("conditional_%s_%d_enrichments_res.csv", para, max_iter)), row.names = FALSE, quote = FALSE)

    # Reorder annotations based on Mean
    Genr_results <- Genr_results[order(Genr_results$Mean, decreasing = TRUE), ]
    Genr_results$Annotation <- factor(
        Genr_results$Annotation,
        levels = Genr_results$Annotation
    )

    # Return results
    return(Genr_results)
}

process_seed_chain <- function(seed, res_dir, trait_name1, trait_name2, annotationName, annot_df, iter_step = NULL, chain_lengths = NULL, min_iter = NULL, single_trait = FALSE, weight_sum = NULL) {
    dir_seed <- file.path(res_dir, sprintf("%s_%s_seed%s", trait_name1, trait_name2, seed))
    print(sprintf("Processing seed %d in directory %s", seed, dir_seed))
    nCat <- length(annotationName)
    print(paste("Number of annotations:", nCat))
    if (single_trait) {
        mcmc_file_h21 <- file.path(dir_seed, "Trait1/MCMC_samples_genetic_effects_variance.txt")
        mcmc_file_h22 <- file.path(dir_seed, "Trait2/MCMC_samples_genetic_effects_variance.txt")
        h21_df <- read_h2_singletrait_samples(mcmc_file_h21, iter_step = iter_step)
        h22_df <- read_h2_singletrait_samples(mcmc_file_h22, iter_step = iter_step)
        cov_df <- NULL
    } else {
        mcmc_file <- file.path(dir_seed, "MCMC_samples_genetic_effects_variance.txt")
        h21_df <- read_h2_samples(mcmc_file, n_annot = nCat, iter_step = iter_step, trait = 1)
        h22_df <- read_h2_samples(mcmc_file, n_annot = nCat, iter_step = iter_step, trait = 2)
        cov_df <- read_covariance_samples(mcmc_file, n_annot = nCat, iter_step = iter_step)
        mcmc_file_total <- file.path(dir_seed, "MCMC_samples_total_genetic_effects_variance.txt")
        total_df <- read_total_gcov_gcor_h2_samples(mcmc_file_total, iter_step = iter_step)
    }

    # Run enrichment
    enrich_results <- function(df, df_total, param_name) {
        result_list <- lapply(chain_lengths, function(len) {
            run_enrichment_by_chain_length(df, df_total, len, min_iter, annotationName, annot_df, dir_seed, param_name, weight_sum)
        })
        names(result_list) <- paste0(param_name, "_", chain_lengths / 1000, "K")
        return(result_list)
    }

    res_list <- list(
        seed = seed,
        h21 = enrich_results(h21_df, total_df, "h21"),
        h22 = enrich_results(h22_df, total_df, "h22")
    )
    if (!is.null(cov_df)) {
      res_list$coh2 <- enrich_results(cov_df, total_df, "coh2")
    }
    return(res_list)
}

extract_param_df <- function(all_results, param = "coh2", chains = NULL, seed_filter = NULL) {
    dfs <- list()

    # Create a filtered version without modifying all_results
    filtered_results <- if (!is.null(seed_filter)) {
        Filter(function(res) res$seed == seed_filter, all_results)
    } else {
        all_results
    }

    if (is.null(seed_filter)) {
        for (res in filtered_results) {
            for (chain in chains) {
                df <- res[[param]][[paste0(param, "_", chain)]] %>%
                    select(Annotation, Mean) %>%
                    mutate(Annotation = gsub(".bed", "", Annotation)) %>%
                rename(!!paste0(param, "_enr_", chain) := Mean) %>%
                rename_with(~ paste0(.x, "_seed", res$seed), -Annotation)
                dfs <- append(dfs, list(df))
              }
            }
          }else {
            for (res in filtered_results) {
            for (chain in chains) {
                df <- res[[param]][[paste0(param, "_", chain)]] %>%
                    select(Annotation, Mean) %>%
                    mutate(Annotation = gsub(".bed", "", Annotation)) %>%
                rename(!!paste0(param, "_enr_", chain) := Mean) 
                dfs <- append(dfs, list(df))
              }
            }
          }

    # Full join across all (Annotation must exist in all)
    merged_df <- reduce(dfs, full_join, by = "Annotation")
    return(merged_df)
}
total_seeds = c(32)
all_results <- lapply(total_seeds, function(s, chain_lengths = c(1000, 2000, 3000, 4000, 5000, 6000), min_iter = 0, iter_step = 50) {
    process_seed_chain(s, res_dir, trait_name1, trait_name2, annotationName, annot_df, iter_step = iter_step, chain_lengths = chain_lengths, min_iter = min_iter, single_trait = FALSE, weight_sum = weight_sum)
})
merged_coh2_df <- extract_param_df(all_results, "coh2", chains = c("1K", "2K", "3K", "4K","5K", "6K"), seed_filter = 32)
head(merged_coh2_df)
# plot correlations between enrichment estimates across different chain lengths (1K–6K)
coh2_mat <- merged_coh2_df %>%
  select(starts_with("coh2_enr_")) %>%
  # remove the first row
  slice(-1)
cor_mat <- cor(coh2_mat, use = "pairwise.complete.obs")
cor_df <- melt(cor_mat)  # reshape into long format
p_cor_within_chain = ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1),
                       name = "Correlation") +
  theme_minimal(base_size = 12) +
  labs(x = "Chain length", y = "Chain length",
       title = "Correlation  of coheritability enrichment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(res_dir, sprintf("corBtwChains_%s_%s/cor_heatmap_within_chain_seed32.jpg", trait_name1, trait_name2)), plot = p_cor_within_chain, width = 6, height = 5)


merged_coh2_df <- extract_param_df(all_results, "coh2", chains = c("2K"))
head(merged_coh2_df)
merged_coh2_df_75annot = merged_coh2_df
merged_coh2_df_156annot = merged_coh2_df
merged_coh2_df_75annot$n_annot = "Annot75"
merged_coh2_df_156annot$n_annot = "Annot156"
merged_coh2_df = rbind(merged_coh2_df_75annot, merged_coh2_df_156annot)
head(merged_coh2_df)


df_plot <- merged_coh2_df %>% 

  pivot_longer(
    cols = starts_with("coh2_enr"),
    names_to = "Seed",
    values_to = "Value"
  ) 

p1 = ggplot(df_plot, aes(x = Annotation, y = Value)) + 
  geom_jitter(
    size = 1, alpha = 0.7
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Scatterplot of coh2 enrichment",
    x = "Annotation",
    y = "Enrichment"#,
    #color = "Number of Annotations"
  ) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")

ggsave(file.path(res_dir, "corBtwChains_SCZ_EA/enr_scatterplot_Group1.jpg"), plot = p1, width = 10, height = 6)




# Get unique annotation names
annotations <- unique(merged_coh2_df$Annotation)
# Loop over each annotation
walk(annotations, function(ann) {
  print(paste("Processing annotation:", ann))
  # Extract values for this annotation
  values <- merged_coh2_df %>%
    filter(Annotation == ann) %>%
    select(-Annotation) %>%
    as.numeric()
  
  # Split into chains (here assuming 10 per chain, adjust if different)
  df <- data.frame(
    chain1 = values[1:10]#,
    #chain2 = values[16:30],
    #chain3 = values[21:30]
  )
  
  # Reshape to long format
  df_long <- pivot_longer(df, cols = starts_with("chain"),
                          names_to = "Chain", values_to = "Value")
  
  # Make scatter plot
  p <- ggplot(df_long, aes(x = Chain, y = Value, color = Chain)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
    labs(title = paste("Scatter plots for", ann),
         x = "Chain", y = "Value") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")
  
  # Save plot
  ggsave(file.path(res_dir, "corBtwChains_T2D_FG/enr_scatterplot_15chains/",
                   paste0("enr_scatterplot_2chains_", ann, ".jpg")),
         plot = p, width = 6, height = 5)
})


#######################################################
# Construct Enrichment Matrix from Averaged MCMC Samples
#######################################################
build_grouped_enrichment_list_from_mcmc <- function(seed_groups, para, res_dir, saved_dir, trait_name1, trait_name2, annotationName, annot_df,
                                                  iter_step = 100, chain_length = 10000, min_iter = 1000, weight_sum, absGtotal) {
  max_iter <- chain_length
  nCat <- length(annotationName)
  print(paste("Number of annotations:", nCat))
  group_enrichments <- list()

  for (i in seq_along(seed_groups)) {
    seeds <- seed_groups[[i]]
    message("Processing Group ", i, " with seeds: ", paste(seeds, collapse = ", "))

    #Read MCMC samples
    mcmc_list <- lapply(seeds, function(seed) {
      mcmc_file <- file.path(res_dir, sprintf("%s_%s_seed%d/MCMC_samples_genetic_effects_variance.txt", trait_name1, trait_name2, seed))
      if (para == "coh2") {
        read_covariance_samples(mcmc_file, n_annot = nCat, iter_step = iter_step)
      } else {
        trait <- ifelse(para == "h21", 1, 2)
        read_h2_samples(mcmc_file, n_annot = nCat, iter_step = iter_step, trait = trait)
      }
    })

    # Average MCMC samples across chains
    sample_mat_list <- lapply(mcmc_list, function(df) {
      df_filtered <- df[df$Iteration %% iter_step == 0 & df$Iteration <= max_iter & df$Iteration > min_iter, ]
      df_filtered[, -1]  # remove Iteration column
    })

    avg_sample_mat <- Reduce("+", sample_mat_list) / length(sample_mat_list)
    n_samples <- (chain_length - min_iter) / iter_step
    iterations <- seq_len(n_samples) * iter_step + min_iter
    avg_df <- data.frame(Iteration = iterations, avg_sample_mat)
    tmp_dir <- sprintf("%s/Group%d_%dchains_enr",saved_dir, i, length(seeds))
    dir.create(tmp_dir, showWarnings = FALSE)
    write.csv(avg_df, file.path(tmp_dir, sprintf("Group%d_%dchains_avg_%s_MCMC.csv", i, length(seeds), para)), row.names = FALSE)
    
    # Read MCMC total gcor samples
    mcmc_list_total <- lapply(seeds, function(seed) {
      mcmc_file <- file.path(res_dir, sprintf("%s_%s_seed%d/MCMC_samples_total_genetic_effects_variance.txt", trait_name1, trait_name2, seed))
      read_total_gcov_gcor_h2_samples(mcmc_file, iter_step = iter_step)
    })

    # Extract matrices, remove Iteration column
    sample_mat_list_total <- lapply(mcmc_list_total, function(df) {
      df_filtered <- df[df$Iteration %% iter_step == 0 & df$Iteration <= max_iter & df$Iteration > min_iter, ]
      as.matrix(df_filtered[, -1])
    })

    # Combine into [iteration, chain]
    n_chain <- length(sample_mat_list_total)
    n_iter <- nrow(sample_mat_list_total[[1]])
    n_metrics <- ncol(sample_mat_list_total[[1]])

    sample_array_total <- array(NA, dim = c(n_iter, n_metrics, n_chain))
    for (k in seq_len(n_chain)) {
      sample_array_total[, , k] <- sample_mat_list_total[[k]]
    }

    # Mean across chains for each iteration
    avg_total_df <- apply(sample_array_total, c(1, 2), mean)
    colnames(avg_total_df) <- colnames(mcmc_list_total[[1]])[-1]
    avg_total_df <- data.frame(Iteration = iterations, avg_total_df)

    # Compute enrichment
    enrich_df <- run_enrichment_by_chain_length(avg_df, avg_total_df, max_iter, min_iter, annotationName, annot_df, out_dir = tmp_dir, para, weight_sum, absGtotal = absGtotal)
    enrich_df$Annotation <- gsub("\\.bed$", "", enrich_df$Annotation)
    group_enrichments[[paste0("Group", i)]] <- enrich_df
  }
  return(group_enrichments)
}


build_grouped_genetic_correlation_list_from_mcmc <- function(seed_groups, res_dir, saved_dir, trait_name1, trait_name2, annotationName, annot_df,
                                                  iter_step = 100, chain_length = 10000, min_iter = 1000, isgcor = TRUE) {
  max_iter <- chain_length
  nCat <- length(annotationName)
  print(paste("Number of annotations:", nCat))
  group_correlation <- list()

  for (i in seq_along(seed_groups)) {
    seeds <- seed_groups[[i]]
    message("Processing Group ", i, " with seeds: ", paste(seeds, collapse = ", "))

    #Read MCMC samples
    mcmc_list <- lapply(seeds, function(seed) {
      mcmc_file <- file.path(res_dir, sprintf("%s_%s_seed%d/MCMC_samples_genetic_effects_variance.txt", trait_name1, trait_name2, seed))
      if (isgcor){
        read_correlation_samples(mcmc_file, n_annot = nCat, iter_step = iter_step)
      } else {
        read_covariance_samples(mcmc_file, n_annot = nCat, iter_step = iter_step)
      }
    })

    # Average MCMC samples across chains
    sample_mat_list <- lapply(mcmc_list, function(df) {
      df_filtered <- df[df$Iteration %% iter_step == 0 & df$Iteration <= max_iter & df$Iteration > min_iter, ]
      df_filtered[, -1]  # remove Iteration column
    })

    avg_sample_mat <- Reduce("+", sample_mat_list) / length(sample_mat_list)
    n_samples <- (chain_length - min_iter) / iter_step
    iterations <- seq_len(n_samples) * iter_step + min_iter
    avg_df <- data.frame(Iteration = iterations, avg_sample_mat)
    tmp_dir <- sprintf("%s/Group%d_%dchains_enr",saved_dir, i, length(seeds))
    dir.create(tmp_dir, showWarnings = FALSE)
    if (isgcor){
      write.csv(avg_df, file.path(tmp_dir, sprintf("Group%d_%dchains_avg_%s_MCMC.csv", i, length(seeds), "gcor")), row.names = FALSE)
    } else {
      write.csv(avg_df, file.path(tmp_dir, sprintf("Group%d_%dchains_avg_%s_MCMC.csv", i, length(seeds), "gcov")), row.names = FALSE)
    }
    group_correlation[[paste0("Group", i)]] <- avg_df
  }
  return(group_correlation)
}


#######################################################
# Plot enrichment barplots for each group using average enrichments across chains
#######################################################
# get enrichment for all T2D trait-pairs 
trait_name1 = "T2D"
trait_name2s = c("AF", "BMI", "CHOL", "DEP", "EA", "Height", "IBD", "INS", "PC", "PD", "RA", "RBC", "SBP", "SCZ")
all_data <- list()

for (trait_name2 in trait_name2s) {
    saved_dir <- sprintf("%s/corBtwChains_%s_%s_absGtotal/", res_dir, trait_name1, trait_name2)
    dir.create(saved_dir, showWarnings = FALSE)

    enrich_list <- build_grouped_enrichment_list_from_mcmc(
      #seed_groups_FG,
      seed_groups,
      "coh2",
      res_dir = res_dir,
      saved_dir = saved_dir,
      trait_name1 = trait_name1,
      trait_name2 = trait_name2,
      annotationName,
      annot_df,
      iter_step = 50,
      chain_length = 2000,
      weight_sum = weight_sum,
      absGtotal = TRUE,
      min_iter = 0
    )

    enriched_df = enrich_list[[1]] # Assuming we only have one group for now
    enriched_df$Trait <- trait_name2
    all_data[[trait_name2]] <- enriched_df
}


combined_df <- bind_rows(all_data)
head(combined_df)
dim(combined_df)
write.csv(combined_df, file.path(res_dir, "concatedGenr_absGtotal.csv"), row.names = FALSE)

saved_dir <- sprintf("%s/corBtwChains_%s_%s/", res_dir, trait_name1, trait_name2)
dir.create(saved_dir, showWarnings = FALSE)
para = "coh2"
# seed_groups = list(Group1 = c(30:39))
enrich_list <- build_grouped_enrichment_list_from_mcmc(
#enrich_list <- build_grouped_enrichment_list_from_mcmc_totalSamples(
  seed_groups,
  para,
  res_dir = res_dir,
  saved_dir = saved_dir,
  trait_name1 = trait_name1,
  trait_name2 = trait_name2,
  annotationName,
  annot_df,
  iter_step = 50,
  chain_length = 2000,
  weight_sum = weight_sum,
  absGtotal = FALSE, 
  min_iter = 0
)
enrich_6K = enrich_list$Group1
enrich_3K = enrich_list$Group1
cor(enrich_6K$Mean, enrich_3K$Mean)

selected_anno = enrich_3K[enrich_3K$EnrPP > 0.9,"Annotation"]

# enrich_list$Group1
write.csv(enrich_list$Group1, sprintf("%s/Group1_10chains_enr/%s_enrichment_results_Group1.csv", saved_dir, para), row.names = FALSE, quote = FALSE)
write.csv(enrich_list$Group2, sprintf("%s/Group2_10chains_enr/%s_enrichment_results_Group2.csv", saved_dir, para), row.names = FALSE, quote = FALSE)
write.csv(enrich_list$Group3, sprintf("%s/Group3_10chains_enr/%s_enrichment_results_Group3.csv", saved_dir, para), row.names = FALSE, quote = FALSE)


saved_dir <- sprintf("%s/corBtwChains_%s_%s/", res_dir, trait_name1, trait_name2)
dir.create(saved_dir, showWarnings = FALSE)
gcor_list <- build_grouped_genetic_correlation_list_from_mcmc(
  seed_groups,
  res_dir = res_dir,
  saved_dir = saved_dir,
  trait_name1 = trait_name1,
  trait_name2 = trait_name2,
  annotationName,
  annot_df,
  iter_step = 50,
  chain_length = 2000,
  min_iter = 0,
  isgcor = TRUE
)
colnames(gcor_list$Group1)[-1] <- annotationName
gcor_df <- gcor_list$Group1
# column means of gcor_df excluding Iteration column
gcor_means <- colMeans(gcor_df[, -1])
gcor_sds <- apply(gcor_df[, -1], 2, sd)
gcor_summary <- data.frame(
  Annotation = names(gcor_means),
  gcor = gcor_means,
  gcor_sd = gcor_sds
)

gcov_list <- build_grouped_genetic_correlation_list_from_mcmc(
  seed_groups,
  res_dir = res_dir,
  saved_dir = saved_dir,
  trait_name1 = trait_name1,
  trait_name2 = trait_name2,
  annotationName,
  annot_df,
  iter_step = 50,
  chain_length = 2000,
  min_iter = 0,
  isgcor = FALSE
)
colnames(gcov_list$Group1)[-1] <- annotationName
gcov_df <- gcov_list$Group1
# column means of gcov_df excluding Iteration column
gcov_means <- colMeans(gcov_df[, -1])
gcov_sds <- apply(gcov_df[, -1], 2, sd)
gcov_summary <- data.frame(
  Annotation = names(gcov_means),
  gcov = gcov_means,
  gcov_sd = gcov_sds
)

# combine gcor_summary and gcov_summary
gcor_gcov_summary <- merge(gcor_summary, gcov_summary, by = "Annotation")
gcor_gcov_summary <- gcor_gcov_summary[order(gcor_gcov_summary$gcor, decreasing = TRUE), ]
write.csv(gcor_gcov_summary, file = file.path(saved_dir, sprintf("Group1_10chains_enr/gcor_gcov_results_Group1_%s_%s.csv", trait_name1, trait_name2)), row.names = FALSE, quote = FALSE)
# sort by gcor
gcor_gcov_summary
#######################################################
# Group specific SNP effect covaraince
#######################################################
build_grouped_SNPcov_list_from_mcmc <- function(seed_groups, res_dir, trait_name1, trait_name2,
                                              annotationName, n_annot,
                                              iter_step = 50, chain_length = 3000, min_iter = 1000, compute_cor = FALSE) {
  max_iter <- chain_length
  group_SNPcov_list <- list()

  for (i in seq_along(seed_groups)) {
    seeds <- seed_groups[[i]]
    message("Processing Group ", i, " with seeds: ", paste(seeds, collapse = ", "))

    # Read SNPcov samples for each seed
    mcmc_list <- lapply(seeds, function(seed) {
      mcmc_file <- file.path(res_dir, sprintf("%s_%s_seed%d/MCMC_samples_marker_effects_variance.txt", trait_name1, trait_name2, seed))
      if (compute_cor) {
        read_correlation_samples(mcmc_file, n_annot = n_annot, iter_step = iter_step)
      } else {
        # Read covariance samples
        read_covariance_samples(mcmc_file, n_annot = n_annot, iter_step = iter_step)
      }
    })

    # Filter by iteration and remove Iteration column
    chain_matrices <- lapply(mcmc_list, function(df) {
      df_filtered <- df[df$Iteration %% iter_step == 0 & df$Iteration <= max_iter & df$Iteration > min_iter, ]
      as.matrix(df_filtered[, -1]) # drop Iteration column
    })

    #Check dimensions
    n_chain <- length(chain_matrices)
    n_iter <- nrow(chain_matrices[[1]])
    n_annot <- ncol(chain_matrices[[1]])

    # Combine manually into 3D array: [iteration, annotation, chain]
    sample_array <- array(NA, dim = c(n_iter, n_annot, n_chain))
    for (k in seq_len(n_chain)) {
      sample_array[, , k] <- chain_matrices[[k]]
    }

    # mean across chains for each [iteration, annotation]
    per_iter_chain_avg <- apply(sample_array, c(1, 2), mean, na.rm = TRUE)

    # mean and sd across iterations
    SNPcov_mean <- colMeans(per_iter_chain_avg)
    SNPcov_sd <- apply(per_iter_chain_avg, 2, sd)

    # Final result for this trait and group
    SNPcov_df <- data.frame(
      Trait = sprintf("%s_%s", trait_name1, trait_name2),
      Annotation = gsub("\\.bed$", "", annotationName),
      Mean = SNPcov_mean,
      SD = SNPcov_sd
    )

    # Use trait name as part of the list name
    group_label <- sprintf("Trait_%s_%s_Group_%d", trait_name1, trait_name2, i)
    group_SNPcov_list[[group_label]] <- SNPcov_df
  }

  return(group_SNPcov_list)
}

# SNPcov_df = build_grouped_SNPcov_list_from_mcmc(seed_groups,
#   res_dir, trait_name, annotationName,
#   iter_step = 50, chain_length = 3000)[[1]]
seed_groups <- list(
  Group1 = total_seeds_T2D_FG[1:10]
)
seed_groups <- list(
  Group1 = total_seeds_T2D_FG[1:15]
)
SNPcor_df = build_grouped_SNPcov_list_from_mcmc(seed_groups,
  res_dir, trait_name1, trait_name2, annotationName,n_annot = nCat,
  iter_step = 50, chain_length = 2000, min_iter = 0, compute_cor = TRUE
)[[1]]
write.csv(SNPcor_df, file = file.path(saved_dir, sprintf("Group1_10chains_enr/SNPcor_results_Group1.csv")), row.names = FALSE, quote = FALSE)


#######################################################
# Group specific SNP Pi
#######################################################
summarize_pi <- function(pi_array, annotationName, label, min_iter, max_iter, iter_step) {
  # Determine valid iteration indices
  n_samples <- dim(pi_array)[1]
  iter_sample <- (1:n_samples) * iter_step
  iter_indices <- which(iter_sample > min_iter & iter_sample <= max_iter)

  # Subset the array
  pi_array_sub <- pi_array[iter_indices, , , drop = FALSE]

  # Compute average across chains first → mean over chains per iteration
  pi_mean_mat <- apply(pi_array_sub, c(1, 2), mean, na.rm = TRUE) # [iter, annot]

  # Now get summary over iterations
  mean_vec <- colMeans(pi_mean_mat, na.rm = TRUE) # length = n_annot
  sd_vec <- apply(pi_mean_mat, 2, sd, na.rm = TRUE)

  data.frame(
    Annotation = annotationName,
    Mean = mean_vec,
    SD = sd_vec,
    PiType = label
  )
}

build_grouped_SNPpi_df <- function(pi_files, n_annot, annotationName, iter_step = 50, chain_length = 3000, min_iter = 1000) {
  n_iter <- chain_length / iter_step
  n_chain <- length(pi_files)
  max_iter <- chain_length

  # 3D arrays: [iter, annot, chain]
  pi00_array <- array(NA, dim = c(n_iter, n_annot, n_chain))
  pi11_array <- array(NA, dim = c(n_iter, n_annot, n_chain))
  pi10_array <- array(NA, dim = c(n_iter, n_annot, n_chain))
  pi01_array <- array(NA, dim = c(n_iter, n_annot, n_chain))
  pi11_standardized_array <- array(NA, dim = c(n_iter, n_annot, n_chain))
  pi10_standardized_array <- array(NA, dim = c(n_iter, n_annot, n_chain))
  pi01_standardized_array <- array(NA, dim = c(n_iter, n_annot, n_chain))

  for (c in seq_along(pi_files)) {
    lines <- readLines(pi_files[c])
    probs <- as.numeric(sub(".*,(.*)", "\\1", lines)) # extract probability

    if (length(probs) != n_iter * n_annot * 4) {
      stop("Number of lines in file does not match expected [n_iter * n_annot * 4]")
    }

    # Reshape to: [n_iter, n_annot, 4]
    prob_matrix <- array(probs, dim = c(4, n_annot, n_iter))
    prob_matrix <- aperm(prob_matrix, c(3, 2, 1)) # now [n_iter, n_annot, 4]

    pi00_array[, , c] <- prob_matrix[, , 1]
    pi11_array[, , c] <- prob_matrix[, , 2]
    pi10_array[, , c] <- prob_matrix[, , 3]
    pi01_array[, , c] <- prob_matrix[, , 4]
    pi11_standardized_array[, , c] <- pi11_array[, , c] / (1 - pi00_array[, , c])
    pi10_standardized_array[, , c] <- pi10_array[, , c] / (1 - pi00_array[, , c])
    pi01_standardized_array[, , c] <- pi01_array[, , c] / (1 - pi00_array[, , c])
  }

  # Summarize all 4
  df_pi00 <- summarize_pi(pi00_array, annotationName, "Pi00", min_iter, max_iter, iter_step)
  df_pi11 <- summarize_pi(pi11_array, annotationName, "Pi11", min_iter, max_iter, iter_step)
  df_pi10 <- summarize_pi(pi10_array, annotationName, "Pi10", min_iter, max_iter, iter_step)
  df_pi01 <- summarize_pi(pi01_array, annotationName, "Pi01", min_iter, max_iter, iter_step)

  df_pi11_standardized <- summarize_pi(pi11_standardized_array, annotationName, "Pi11_standardized", min_iter, max_iter, iter_step)
  df_pi10_standardized <- summarize_pi(pi10_standardized_array, annotationName, "Pi10_standardized", min_iter, max_iter, iter_step)
  df_pi01_standardized <- summarize_pi(pi01_standardized_array, annotationName, "Pi01_standardized", min_iter, max_iter, iter_step)

  # Combine
  final_df <- rbind(
    df_pi00, df_pi11, df_pi10, df_pi01,
    df_pi11_standardized, df_pi10_standardized, df_pi01_standardized
  )
  return(final_df)
}

pi_files <- sprintf("%s/%s_%s_seed%s/MCMC_samples_pi.txt", res_dir, trait_name1, trait_name2, total_seeds)
pi_summary_df <- build_grouped_SNPpi_df(
  pi_files,
  n_annot = nCat,
  annotationName = annotationName,
  iter_step = 50,
  chain_length = 2000,
  min_iter = 0
)
pi_summary_df$Annotation <- gsub("\\.bed$", "", pi_summary_df$Annotation)

# ###################################################################
# #### Get the pi_summary_df for all traits
# ###################################################################
# Placeholder to collect all data
trait_name1 = "T2D"
trait_name2s = c("AF","BMI", "CHOL", "DEP","EA", "Height", "IBD", "INS", "PC", "PD", "RA", "RBC", "SBP", "SCZ")

trait_name_pairs = paste0(trait_name1, "_", trait_name2s)
combined_pi_summary_df <- list()

for (trait_pair in trait_name_pairs) {
  print(paste("Processing trait pair:", trait_pair))

  pi_files <- paste0(
    res_dir,  trait_pair, "_seed", total_seeds,
    "/MCMC_samples_pi.txt"
  )

  pi_summary_df <- build_grouped_SNPpi_df(
    pi_files,
    n_annot = nCat,
    annotationName = annotationName,
    iter_step = 50,
    chain_length = 2000,
    min_iter = 0
  )

  pi_summary_df$Annotation <- gsub("\\.bed$", "", pi_summary_df$Annotation)
  pi_summary_df$trait_pair <- trait_pair

  combined_pi_summary_df[[trait_pair]] <- pi_summary_df
}

# Combine into one data frame
final_pi_df <- bind_rows(combined_pi_summary_df)
final_pi_wide_df <- final_pi_df %>%
  pivot_wider(
    names_from = PiType,
    values_from = c(Mean, SD),
    names_glue = "{PiType}_{.value}"
  )
nrow(final_pi_wide_df)
write.csv(final_pi_wide_df, file.path(res_dir, "concatedSNPpi_results.csv"), row.names = FALSE)

##########################################################
# Plot mirror plot for Pi11 and SNPcov
##########################################################
# Pivot pi_summary_df to wide format
pi_wide_df <- pi_summary_df %>%
  pivot_wider(
    names_from = PiType,
    values_from = c(Mean, SD),
    names_glue = "{PiType}_{.value}"
  )
write.csv(pi_wide_df, file = file.path(saved_dir, sprintf("Group1_10chains_enr/SNPpi_results_Group1.csv")), row.names = FALSE, quote = FALSE)


#######################################################
# Get total genetic correlations dataframe
#######################################################
read_total_gcov_gcor_h2_samples <- function(mcmc_file, iter_step = 50) {
  lines <- readLines(mcmc_file)
  matrix_values <- do.call(rbind, strsplit(lines, ","))
  matrix_values <- apply(matrix_values, 2, as.numeric)

  stopifnot(nrow(matrix_values) %% 2 == 0)

  n_samples <- nrow(matrix_values) / 2
  result <- matrix(NA, nrow = n_samples, ncol = 5) # Now includes iteration

  for (i in 1:n_samples) {
    row_start <- (i - 1) * 2 + 1
    mat <- matrix(matrix_values[row_start:(row_start + 1), ], nrow = 2, byrow = TRUE)

    varg1 <- mat[1, 1]
    varg2 <- mat[2, 2]
    gcov <- mat[1, 2]
    gcor <- if (varg1 > 0 && varg2 > 0) gcov / sqrt(varg1 * varg2) else NA

    result[i, ] <- c(i * iter_step, gcov, varg1, varg2, gcor)
  }

  df <- as.data.frame(result)
  colnames(df) <- c("Iteration", "Gcov", "H2_1", "H2_2", "Gcor")
  return(df)
}

build_total_summary_list_from_mcmc <- function(seed_groups, res_dir, trait_pair,
                                            iter_step = 50, chain_length = 3000, min_iter = 1000) {
  max_iter <- chain_length
  group_summary_list <- list()

  for (i in seq_along(seed_groups)) {
    seeds <- seed_groups[[i]]
    message("Processing Group ", i, " for Trait Pair ", trait_pair)

    # Read MCMC total gcor samples
    mcmc_list <- lapply(seeds, function(seed) {
      mcmc_file <- file.path(res_dir, sprintf("%s_seed%d/MCMC_samples_total_genetic_effects_variance.txt", trait_pair, seed))
      read_total_gcov_gcor_h2_samples(mcmc_file, iter_step = iter_step)
    })

    # Extract matrices, remove Iteration column
    chain_matrices <- lapply(mcmc_list, function(df) {
      df_filtered <- df[df$Iteration %% iter_step == 0 & df$Iteration <= max_iter & df$Iteration > min_iter, ]
      as.matrix(df_filtered[, -1])
    })

    # Combine into [iteration, chain]
    n_chain <- length(chain_matrices)
    n_iter <- nrow(chain_matrices[[1]])
    n_metrics <- ncol(chain_matrices[[1]])

    sample_array <- array(NA, dim = c(n_iter, n_metrics, n_chain))
    for (k in seq_len(n_chain)) {
      sample_array[, , k] <- chain_matrices[[k]]
    }

    # Mean across chains for each iteration
    per_iter_chain_avg <- apply(sample_array, c(1, 2), mean)

    # Mean aross iterations for each chain
    per_chain_iter_avg <- apply(sample_array, c(2, 3), mean)
    per_chain_iter_avg_gcor_samples <- per_chain_iter_avg[4, ]
    # scatter plot of per_chain_iter_avg_gcor_samples using ggplot
    

    # Final mean and sd
    metric_mean <- colMeans(per_iter_chain_avg)
    metric_sd <- apply(per_iter_chain_avg, 2, sd)

    # Create result df
    summary_df <- data.frame(
      Trait = trait_pair,
      Metric = c("Gcov", "H2_1", "H2_2", "Gcor"),
      Mean = metric_mean,
      SD = metric_sd
    )
    group_label <- paste0("Trait_", trait_pair, "_Group", i)
    group_summary_list[[group_label]] <- summary_df
  }
  return(group_summary_list)
}

# Placeholder to collect all data
trait_name1 = "T2D"
trait_name2s = c("AF", "BMI", "CHOL", "DEP", "EA", "Height", "IBD", "INS", "PC", "PD", "RA", "RBC", "SBP", "SCZ")
trait_pairs = c(paste0(trait_name1, "_", trait_name2s))
#trait_pairs = c(sprintf("%s_%s", trait_name1, trait_name2))

all_total_summary_lists <- lapply(trait_pairs, function(tr) {
  build_total_summary_list_from_mcmc(seed_groups, res_dir, tr, iter_step = 50, chain_length = 2000, min_iter = 0)
})
#all_total_summary_lists[[length(all_total_summary_lists)+1]] <- build_total_summary_list_from_mcmc(seed_groups_FG, res_dir, "T2D_FG", iter_step = 50, chain_length = 2000, min_iter = 0)
combined_total_df <- do.call(rbind, unlist(all_total_summary_lists, recursive = FALSE))
dim(combined_total_df)

combined_total_df[combined_total_df$Metric == "Gcor",]
write.csv(combined_total_df, file = file.path(saved_dir, sprintf("Group1_10chains_enr/grouped_total_summary_results_%s.csv", trait_pairs)), row.names = FALSE)
write.csv(combined_total_df, file = file.path(res_dir, sprintf("grouped_total_summary_results_T2D.csv")), row.names = FALSE)

MCMC_path1 = "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/analysis/real_data/MTwoAnnot/MTSBayesCC_3K_corrR_tuned_v3_estSigma/T2D_FG_seed123/MCMC_samples_genetic_effects_variance.txt"
samples1 = read_total_gcov_gcor_h2_samples(MCMC_path1)
samples1_af_burnin = samples1[samples1$Iteration > 1000,] 
# # column mean after burn-in
apply(samples1_af_burnin[2:ncol(samples1_af_burnin)],2,mean)
apply(samples1_af_burnin[2:ncol(samples1_af_burnin)],2,sd)






plot_corr_heatmap <- function(df, para, saved_dir, which_chains = NULL, seeds = NULL, exclude_top_annotation = FALSE) {
  # Define seeds and desired column order
  ordered_colnames <- unlist(lapply(seeds, function(s) {
      paste0(para, "_enr_", which_chains, "_seed", s)
  }))

  # Keep only columns that exist
  selected_cols <- ordered_colnames[ordered_colnames %in% colnames(df)]

  # Subset and convert to matrix
  df_numeric <- df %>%
      select(all_of(selected_cols)) %>%
      mutate(across(everything(), as.numeric)) 

  # Optionally exclude the most enriched annotation
  if (exclude_top_annotation) {
      max_row <- which.max(apply(df_numeric, 1, max, na.rm = TRUE))
      df_numeric <- df_numeric[-max_row, , drop = FALSE]
  }

  # Convert to matrix and compute correlation
  enr_mat <- as.matrix(df_numeric)
  cor_mat <- cor(enr_mat, use = "pairwise.complete.obs")
  cor_df <- reshape2::melt(cor_mat)
  cor_df$label <- sprintf("%.2f", cor_df$value)

  # Set factor levels to control plot order
  cor_df$Var1 <- factor(cor_df$Var1, levels = selected_cols)
  cor_df$Var2 <- factor(cor_df$Var2, levels = selected_cols)

  # Plot
  p <- ggplot(cor_df, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = label), size = 2.5) +
      scale_fill_gradient2(
          low = "blue", high = "red", mid = "white",
          midpoint = 0, limit = c(-1, 1),
          name = "Correlation"
      ) +
      theme_minimal() +
      theme(
          axis.text.x = element_text(size = 10, angle = 30, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_blank(),
          legend.position = "right",
          plot.margin = unit(c(1, 1, 2, 1), "cm")
      )

  # Save filename
  extra_suffix <- if (exclude_top_annotation) "_rmTop" else ""
  if (length(seeds) > 1) {
      filename <- sprintf("%s/cor_%s_enr_%s%s.png", saved_dir, para, tolower(paste(which_chains, collapse = "_")), extra_suffix)
  } else {
      filename <- sprintf("%s/cor_%s_enr_%s_seed%d%s.png", saved_dir, para, tolower(paste(which_chains, collapse = "_")), seeds[1], extra_suffix)
  }
  ggsave(filename, plot = p, width = 8, height = 6)

  return(p)
}






