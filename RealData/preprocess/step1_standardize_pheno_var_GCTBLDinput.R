library(data.table)

# Parsing arguments for data path
args <- commandArgs(trailingOnly = TRUE)
gwas_path <- args[1] # /group/qtlchenggrp/jiayi/MTSBayesC/data/UKBioBank/RealData/T2D_Comorbidities/GWAS_summary_statistics/ImputedSumstat/
                     # /common/zhao/tianjing/GWAS_summary_statistics/ImputedSumstat/
LD_info_path <- args[2] # /group/qtlchenggrp/jiayi/MTSBayesC/data/UKBioBank/GCTB_ldm/ukbEUR_HM3/
                        # /common/zhao/jyqqu/MTSBayesCC/data/eigen_data_whole_genome/
preprocess_file_path <- args[3] # /group/qtlchenggrp/jiayi/MTSBayesC/data/UKBioBank/RealData/T2D_Comorbidities/T2D_SBP/
                                # /common/zhao/jyqqu/MTSBayesCC/data/real_data/T2D_FG/
trait1_file <- args[4] # T2D.imputed.ma
trait2_file <- args[5] # FG.imputed.convertsign.ma
LDinfo_file <- args[6] # snp.info

# read gwas data files
Trait1 <- fread(paste0(gwas_path, trait1_file))
Trait2 <- fread(paste0(gwas_path, trait2_file))
LDinfo <- fread(paste0(LD_info_path, LDinfo_file), header = TRUE)

# After GCTB Imputation, the SNPs in Traiti and LDinfo matches with nSNP = 1154522
cat("nSNPs in LD matrix:", nrow(LDinfo), "\n")
cat("nSNPs having bhat:", nrow(Trait1), "\n")

# Preparing folder to save SNPs per block 
dir.create(paste0(LD_info_path, "SNPsPerBlock/"), recursive = TRUE, showWarnings = FALSE)

nBlocks <- 591

# Function to adjust SNP effects - standardize the phenotypic variance
# The order of bhat and BlockInfo is consistent
adjust_effects <- function(bhat, BlockInfo) {
    bhat$freq <- ifelse(bhat$A1 == BlockInfo$A2, 1 - bhat$freq, bhat$freq)
    bhat$b <- ifelse(bhat$A1 == BlockInfo$A2, -bhat$b, bhat$b)
    bhat[, c("A1", "A2") := .(BlockInfo$A1, BlockInfo$A2)]
    bhat[, sj := sqrt(1 / (N * se^2 + b^2))]
    bhat[, bAdj := b * sj]
    bhat[, seAdj := se * sj]
    bhat[, D := 2 * freq * (1 - freq) * N]
    bhat[, varps := N * seAdj^2 + bAdj^2]
    return(bhat)
}

for (i in 1:nBlocks) {
    print(sprintf("block %d", i))
    # read the LD info for block i
    BlockInfoi <- LDinfo[LDinfo$Block == i, ]
    SNPsLDOrder <- BlockInfoi$ID

    bhat1i <- Trait1[Trait1$SNP %in% BlockInfoi$ID, ]
    bhat2i <- Trait2[Trait2$SNP %in% BlockInfoi$ID, ]
    # make the order of BlockInfoi consistent with bhat1i and bhat2i
    BlockInfoi <- BlockInfoi[match(bhat1i$SNP, BlockInfoi$ID), ]

    if (all(BlockInfoi$ID == bhat1i$SNP) != TRUE) {
        stop("Error: SNPs in BlockInfoi and bhat1 are not the same")
    }
    if (all(BlockInfoi$ID == bhat2i$SNP) != TRUE) {
        stop("Error: SNPs in BlockInfoi and bhat2 are not the same")
    }

    cat(length(BlockInfoi$ID), "SNPs in common with LD information \n")

    # Adjusting SNP effects
    bhat1i <- adjust_effects(bhat1i, BlockInfoi)
    bhat2i <- adjust_effects(bhat2i, BlockInfoi)

    # # Convert to data.table if not already
    setDT(bhat1i)
    setDT(bhat2i)

    # Ensure both data.tables have 'SNP' as a key for joining
    setkey(Trait1, SNP)
    setkey(bhat1i, SNP)

    setkey(Trait2, SNP)
    setkey(bhat2i, SNP)

    # Perform an update join.
    # This updates Trait1/Trait2 directly with values from bhat1i/bhat2i for matching SNPs.
    Trait1[bhat1i, `:=`(
        bAdj = i.bAdj,
        varps = i.varps,
        A1 = i.A1,
        A2 = i.A2,
        b = i.b,
        freq = i.freq
    ), on = .(SNP)]

    Trait2[bhat2i, `:=`(
        bAdj = i.bAdj,
        varps = i.varps,
        A1 = i.A1,
        A2 = i.A2,
        b = i.b,
        freq = i.freq
    ), on = .(SNP)]

    write.table(SNPsLDOrder, paste0(LD_info_path, sprintf("SNPsPerBlock/SNPs_block%d.csv", i)), col.names = FALSE, row.names = FALSE)
}

dir.create(paste0(preprocess_file_path, "standardPhenoVar/"), recursive = TRUE)
write.csv(Trait1, file = paste0(preprocess_file_path, "standardPhenoVar/b1_complete.txt"), row.names = FALSE) 
write.csv(Trait2, file = paste0(preprocess_file_path, "standardPhenoVar/b2_complete.txt"), row.names = FALSE)

message("mean(b1$varps) = ", mean(Trait1$varps))
message("mean(b2$varps) = ", mean(Trait2$varps))