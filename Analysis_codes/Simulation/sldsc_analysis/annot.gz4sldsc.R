options(scipen = 999)
library(MASS)
library(data.table)
library(optparse)

##############################################################
#####         Step 1. Parameter setting option         #######
##############################################################
option_list <- list(
    ##  type
    make_option(c("--seed"),
        action = "store", default = 123, type = "numeric",
        help = ""
    ),
    make_option(c("--pleio_percent"),
        action = "store", default = 0.1, type = "numeric",
        help = ""
    ), 
    make_option(c("--sample_size"),
        action = "store", default = 10, type = "numeric",
        help = ""
    ),
    # make_option(c("--annotation_size"),
    #     action = "store", default = 0.1, type = "numeric",
    #     help = ""
    # ),
    make_option(c("--h21"),
        action = "store", default = 0.01, type = "numeric",
        help = ""
    ),
    make_option(c("--h22"),
        action = "store", default = 0.01, type = "numeric",
        help = ""
    )
)
opt <- parse_args(OptionParser(option_list = option_list))

print("option is: ")
opt

sample_size <- opt$sample_size
annotation_size <- opt$annotation_size
pleio_percent <- opt$pleio_percent
seed <- opt$seed
h21 <- opt$h21 # h2 for trait 1
h22 <- opt$h22 # h2 for trait 2

plink19 <- "/home/zhao/jyqqu/plink-1.07-x86_64/plink"
data_folder_path <- "/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/"
data_folder_name <- paste0("h2_trait1.", h21, ".h2_trait2.", h22, "_pleioPercent", pleio_percent, "_sampleSize", sample_size, "_seed", seed)
data_folder = paste0(data_folder_path, data_folder_name, "/") 

setwd(data_folder)
dir.create("sldsc_annot")
# Concatenate files SNPsc1.txt and SNPsc2.txt
system("cat ../SNPc1.txt ../SNPc2.txt > snp_list.txt")

# Read the SNP list into R
snp_list <- readLines("snp_list.txt")

# Read SNPs from each annotation file
snps_annot1 <- readLines("../SNPc1.txt")
snps_annot2 <- readLines("../SNPc2.txt")

bfiles_folder = "/common/zhao/jyqqu/MTSBayesCC/data/CELLECT/data/ldsc/1000G_EUR_Phase3_plink/"
saved_bfiles = "/common/zhao/jyqqu/MTSBayesCC/data/bfiles/sldsc/"
if (!file.exists(saved_bfiles)) {
    dir.create(saved_bfiles)
}

chr = 1
# Define file paths
bim_file <- paste0(bfiles_folder, "1000G.EUR.QC.", chr, ".bim")

# Load .bim file
bim_data <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bim_data) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

# Find SNPs in bim file that are also in snp_list.txt
intersect_snps <- bim_data$SNP[bim_data$SNP %in% snp_list]

# Filter bim_data to only include intersecting SNPs
filtered_bim <- bim_data[bim_data$SNP %in% intersect_snps, ]


# Update bed file with filtered SNPs using PLINK
new_bfile_prefix <- paste0(saved_bfiles, "1000G.EUR.QC.merged.", chr)

# check if the new_bfile_prefix exists
if (!file.exists(paste0(new_bfile_prefix, ".bim"))) {
    system(paste0(
        plink19, " --noweb --bfile ", sub("\\.bim$", "", bim_file),
        " --extract snp_list.txt",
        " --make-bed --out ", new_bfile_prefix
    ))   
}

# Create annotation data frame for this chromosome
annot_data <- data.frame(
    SNP = filtered_bim$SNP,
    Annot1 = as.integer(filtered_bim$SNP %in% snps_annot1),
    Annot2 = as.integer(filtered_bim$SNP %in% snps_annot2)
)

annot_data_expanded <- data.frame(
    SNP = bim_data$SNP,
    Annot1 = as.integer(bim_data$SNP %in% snps_annot1),
    Annot2 = as.integer(bim_data$SNP %in% snps_annot2)
)

# Save annotation data frame as a temporary uncompressed file
annot_file <- paste0("sldsc_annot/SimData.", chr, ".annot")
write.table(annot_data[, 2:3], annot_file,
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

annot_file_expanded <- paste0("sldsc_annot/SimData.expanded.", chr, ".annot")
write.table(annot_data_expanded[, 2:3], annot_file_expanded,
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

# Compress the file using gzip
system(paste("gzip", annot_file))
system(paste("gzip", annot_file_expanded))


