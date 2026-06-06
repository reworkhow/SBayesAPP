# Load necessary libraries
library(dplyr)

base_dir = "/common/zhao/jyqqu/MTSBayesCC/data/real_data/"
# Read the sumstats file (b1_formatted.munge.sumstats)

# Read the annotation file (cell_type_annot_human_total_unoverlap.txt)
cell_types <- read.table(sprintf("%s/cell_type_annot_human_total_oneAnnotperSNP_BetaCellPriori.txt", base_dir), header = TRUE, stringsAsFactors = FALSE)

# Read 1000G bim file
bim <- read.table(sprintf("%s/../CELLECT/data/ldsc/1000G_EUR_Phase3_plink_1M/1000G.EUR.QC.merged.filtered.bim", base_dir), stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

# Merge based on the SNP column
merged_data <- bim %>% 
select(SNP) %>% 
left_join(cell_types, by = c("SNP" = "SNP"))

# Remove the SNP column
final_data <- merged_data %>%
  select(-SNP)

# Save the final reordered data without the SNP column
write.table(final_data, sprintf("%s/cell_type_annot_human_total_oneAnnotperSNP_BetaCellPriori.nosnpinfo.reordered4gnova.txt", base_dir), sep = "\t", row.names = FALSE, quote = FALSE)
