# Load necessary library
library(dplyr)

# Define the parent folder
parent_folder <- "/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/"

# Get a list of all sub-folders
sub_folders <- list.dirs(parent_folder, recursive = TRUE, full.names = TRUE)

# Initialize an empty data frame to store results
results <- data.frame(SubFolder = character(), Rho = numeric(), stringsAsFactors = FALSE)

# Loop through each sub-folder
for (folder in sub_folders) {
  file_path <- file.path(folder, "gnova.txt")
  if (file.exists(file_path)) {
    # Read the file and extract rho values
    data <- read.table(file_path, header = TRUE)
    category <- data$annot_name
    rho_values <- data$rho
    folder_name <- basename(folder) # Extract only the folder name
    results <- rbind(results, data.frame(SubFolder = folder_name, Rho = rho_values, Category = category))
  }
}

# Save the concatenated results to a file
write.csv(results, "gnova_rho_values.csv", row.names = FALSE)


