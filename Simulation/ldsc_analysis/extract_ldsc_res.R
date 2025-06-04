library(dplyr)

# Define the root directory
root_dir <- "./"

# Get all subdirectories containing `ldsc/totalrg/ldsc_res.log`
log_files <- list.files(root_dir, pattern = "ldsc_res.log", recursive = TRUE, full.names = TRUE)

# Initialize a data frame to store results
results <- data.frame(Subfolder = character(), h2_trait1 = numeric(), h2_trait2 = numeric(), genetic_corr = numeric(), genetic_cov = numeric(), stringsAsFactors = FALSE)

# Loop through each log file
for (file in log_files) {
    print(file)

    content <- readLines(file)

    # Extract values using regex
    line_index_trait1 <- grep("Heritability of phenotype 1", content)
    h2_trait1_line <- content[line_index_trait1 + 2]
    h2_trait1 <- as.numeric(sub(".*Total Observed scale h2: ([0-9.]+).*", "\\1", h2_trait1_line))

    line_index_trait2 <- grep("Heritability of phenotype 2", content)
    h2_trait2_line <- content[line_index_trait2 + 2]
    h2_trait2 <- as.numeric(sub(".*Total Observed scale h2: ([0-9.]+).*", "\\1", h2_trait2_line))

    line_index_corr <- grep("Genetic Correlation:", content)
    genetic_corr_line <- content[line_index_corr]
    genetic_corr <- as.numeric(sub(".*Genetic Correlation: ([0-9.-]+).*", "\\1", genetic_corr_line))

    line_index_gcov <- grep("Total Observed scale gencov:", content)
    genetic_cov_line <- content[line_index_gcov]
    genetic_cov <- as.numeric(sub(".*Total Observed scale gencov: ([0-9.-]+).*", "\\1", genetic_cov_line))

    # Append results
    results <- results %>%
        add_row(Subfolder = dirname(file), h2_trait1 = h2_trait1, h2_trait2 = h2_trait2, genetic_corr = genetic_corr, genetic_cov = genetic_cov)
}

results$Subfolder <- sub("^.//", "", results$Subfolder) # Remove leading './'
results$Subfolder <- sub("/ldsc/totalrg$", "", results$Subfolder) # Remove trailing '/ldsc/totalrg'


# Write the summarized results to a file
write.csv(results, file = "summary_ldsc_results.csv", row.names = FALSE)
