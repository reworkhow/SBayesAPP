#Simulation with annotation
# 3.4.1. R: generate bivarte-annotation
# whole-genome version

options(scipen = 999)
library(MASS)
library(data.table)
library(optparse)

##############################################################
#####         Step 1. Parameter setting option         #######
##############################################################
option_list = list(
  ##  type
  make_option(c("--seed"), action="store", default=123, type='numeric',
              help=""),
  make_option(c("--pleio_percent"), action="store", default=0.1, type='numeric',
              help=""),
  make_option(c("--sample_size"), action="store", default=10, type='numeric',
              help=""),
  make_option(c("--h21"), action="store", default=0.01, type='numeric',
              help=""),
  make_option(c("--h22"), action="store", default=0.01, type='numeric',
              help="")
)
opt = parse_args(OptionParser(option_list=option_list))

print("option is: ")
opt

###############################################
### Below are parameters to be changed
###############################################

gcta64 <- "/home/uqjzeng1/wd/proj/TJ/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta64"
plink19 <- "/home/uqjzeng1/wd/proj/TJ/software/plink_linux_x86_64_20231211/plink"

#############################################
sample_size <- opt$sample_size
pleio_percent <- opt$pleio_percent
seed <- opt$seed
h21 <- opt$h21 # h2 for trait 1
h22 <- opt$h22 # h2 for trait 2

# parameters
# seed <- 123 # set seed for reproducibility
# sample_size <- 500
# annotation_size <- 0.1 # 10% of SNPs in category 1
# pleio_percent <- 1.0
# h21 <- 0.1 # h2 for trait 1
# h22 <- 0.07 # h2 for trait 2
nCategory <- 2 # number of annotation groups
pis <- rep(0.01, nCategory) # proportion of non-zero effects in each annotation group
QTL_cor <- c(0.8, 0.2)
genetic_cor <- QTL_cor * pleio_percent / (pleio_percent + (1 - pleio_percent) / 2) # Genetic correlation bewteen two traits for the two category -> half nonpleiotropic QTL affect trait 1, the other half affect trait 2
proportion_h2_explained <- c(0.5, 0.5) # proportion of h2 explained by each annotation group
Rescor <- 0.1
nTrait <- 2

repNum <- 1 # number of replicates in gcta


rawDataPath <- paste0("/home/uqjzeng1/wd/proj/TJ/chr1_real_genotype_ukb/") # Folder containing the UKBB bfiles
outPath <- paste0("/home/uqjzeng1/wd/proj/TJ/sim_chr1_output_v2/") # Folder to save simulation results
rawBimFileSuffix <- paste0("ukbV3_eur_unrel_info0.6_maf0.01_hwe1e10_hm3_chr1")
outPath <- paste0(outPath, "h2_trait1.", h21, ".h2_trait2.", h22, "_pleioPercent", pleio_percent, "_sampleSize", sample_size, "_seed", seed, "/")
# Read in the provided SNPs for category 1 and category 2
SNPsc1 <- read.table("/home/uqjzeng1/wd/proj/TJ/sim_chr1_output_v2/SNPc1.txt", header = F)$V1
SNPsc2 <- read.table("/home/uqjzeng1/wd/proj/TJ/sim_chr1_output_v2/SNPc2.txt", header = F)$V1

###############################################
### Above are parameters to be changed
###############################################
# output parameters
cat("seed: ", seed, "\n")
cat("sample size: ", sample_size, "\n")
cat("pleiotropic percent: ", pleio_percent, "\n")
cat("h2 for trait 1: ", h21, "\n")
cat("h2 for trait 2: ", h22, "\n")
cat("number of annotation groups: ", nCategory, "\n")
cat("proportion of non-zero effects in each annotation group: ", paste(pis, collapse = ", "), "\n")
cat("Genetic correlation in each annotation category (for pleiotropic QTLs): ", paste(genetic_cor, collapse = ", "), "\n")
cat("Residual correlation: ", Rescor, "\n")
cat("folder to save simulation results: ", outPath, "\n")
dir.create(outPath)


rawBimFile <- paste0(rawDataPath, rawBimFileSuffix, ".bim")
rawFamFile <- paste0(rawDataPath, rawBimFileSuffix, ".fam")
bim <- fread(rawBimFile)
names(bim) <- c("CHR", "SNP", "POS", "CM", "A1", "A2")
fam <- fread(rawFamFile)
names(fam) <- c("FID", "IID", "PID", "MID","SEX", "PHENOTYPE")

set.seed(seed)
##############################################
# 1. Sample Size (100K or 300K) ⇒ #individuals
# select individuals based on sample size parameters
##############################################
nInd <- nrow(fam) 
nMarker <- nrow(bim) 
cat("number of Inds in raw file: ", nInd, "\n")
cat("number of SNPs in raw file: ", nMarker, "\n")

# Generate random indices for selected individuals without replacement
random_ind_index <- sort(sample(nInd, sample_size, replace = FALSE))
# selected individuals
selected_ind <- fam[random_ind_index, c("FID", "IID")]

selected_ind_file <- paste0(outPath, rawBimFileSuffix, "-selectedInd.txt")
fwrite(data.frame(selected_ind), file = selected_ind_file, col.names = FALSE, sep = "\t")
selected_ind_out <- sprintf("%s%s.%dInd", outPath, rawBimFileSuffix, sample_size)

# plinkInd <- paste0(
#     plink19, " --bfile ", rawDataPath, rawBimFileSuffix,
#     " --keep ", selected_ind_file,
#     " --make-bed ",
#     "--threads 100 ",
#     " --out ", selected_ind_out
# )
# system(plinkInd)

##  Calculate allele frequency
plinkFreq <- paste0(
    plink19, " --bfile ", rawDataPath, rawBimFileSuffix,
    " --keep ", selected_ind_file,
    " --freq ",
    " --out ", selected_ind_out
)
system(plinkFreq)

##############################################
# 2. Annotation size (10%,50%) ⇒ #SNPs in each annotation group
# select SNPs into annotation group based on provided SNP list
##############################################
# # Read in the provided SNPs for category 1 and category 2
nSNPc1 <- length(SNPsc1)
nSNPc2 <- length(SNPsc2)
cat("number of SNPs selected in category 1: ", nSNPc1, "\n")
cat("number of SNPs selected in category 2: ", nSNPc2, "\n")

##############################################
# 3. Pi ⇒ number of QTL in each annotation group
# select QTL into annotation group based on Pi parameters
##############################################
# number of QTLs in each category
nSNPc <- c(nSNPc1, nSNPc2)
nQTL <- floor(nSNPc * pis)
cat("number of QTL in each annotation group: ", paste(nQTL, collapse = ", "), "\n")

##############################################
# 4. %pleiotropic QTL for 2 traits
# select pleiotropic QTL based on %pleiotropic parameters
##############################################
npleioQTL <- floor(nQTL * pleio_percent) # number of pleiotropic QTL in each annotation group
cat("number of pleiotropic QTL in each annotation group: ", paste(npleioQTL, collapse = ", "), "\n")
if (pleio_percent < 1.0) {
    nnonpleioQTL <- nQTL - npleioQTL # number of non-pleiotropic QTL in each annotation group
    cat("number of non-pleiotropic QTL in each annotation group: ", paste(nnonpleioQTL, collapse = ", "), "\n")
}else{
   nnonpleioQTL = c(0,0)
}

##############################################
# 5. QTL effect correlation
##############################################

# only half affect trait 1
nnonpleioQTLT1 <- floor(nnonpleioQTL / 2)
# only half affect trait 2
nnonpleioQTLT2 <- nnonpleioQTL - nnonpleioQTLT1

QTLmean <- c(0, 0) # mean of QTL effect for trait 1 and trait 2

# Calculate QTL covariance matrix for category 1 and category 2 based on h2 and genetic correlation
QTLcovmat <- list()
for (c in 1:nCategory) {
    genetic_var_11 = h21 * proportion_h2_explained[c]
    genetic_var_22 = h22 * proportion_h2_explained[c]
    genetic_cov_12 = sqrt(genetic_var_11 * genetic_var_22) * genetic_cor[c]
    nQTL11 = npleioQTL[c] + nnonpleioQTLT1[c]
    nQTL22 = npleioQTL[c] + nnonpleioQTLT2[c]
    nQTL12 = nQTL21 = npleioQTL[c]
    QTLcovmat[[c]] <- matrix(c(genetic_var_11 / nQTL11, genetic_cov_12 / nQTL12, genetic_cov_12 / nQTL21, genetic_var_22 / nQTL22), nrow = 2, ncol = 2)
    write.table(QTLcovmat[[c]], sprintf("%s/QTLcovmatc%d.txt", outPath, c), row.names = F, col.names = F, quote = F)
    # write.table([genetic_var_11 genetic_cov_12; genetic_cov_12 genetic_var_22], sprintf("%s/genetic_covmatc%d.txt", outPath, c), row.names = F, col.names = F, quote = F)
    cat("QTL covariance matrix for category ", c, ": \n")
    print(QTLcovmat[[c]])
    cat("genetic covariance matrix for category ", c, ": \n")
    print(matrix(c(genetic_var_11, genetic_cov_12, genetic_cov_12, genetic_var_22), nrow = 2, ncol = 2))
}

QTLpos <- list(
    QTLposc1 = sample(SNPsc1, nQTL[1], replace = FALSE), # QTL positions for category 1
    QTLposc2 = sample(SNPsc2, nQTL[2], replace = FALSE) # QTL positions for category 2
) # QTL positions for category 1 and category 2

# pleiotropic QTL positions for category 1 and category 2
if (pleio_percent < 1.0) {
    pleioQTLpos <- list(
        QTLposc1 = sample(QTLpos[["QTLposc1"]], npleioQTL[1], replace = FALSE),
        QTLposc2 = sample(QTLpos[["QTLposc2"]], npleioQTL[2], replace = FALSE)
    ) # pleiotropic QTL positions for category 1 and category 2
} else if (pleio_percent == 1.0) {
    pleioQTLpos <- QTLpos
}
# save pleiotropic QTL positions
if (pleio_percent != 0) {
    for (i in 1:nCategory) {
        write.table(pleioQTLpos[[i]], sprintf("%s/pleioQTLposc%d.txt", outPath, i), row.names = F, col.names = F, quote = F)
    }
}
# pleiotropic QTL effects for category 1 and category 2
if (pleio_percent != 0){
    pleioQTLeff <- list(
        QTLeffc1 = mvrnorm(npleioQTL[1], QTLmean, QTLcovmat[[1]]),
        QTLeffc2 = mvrnorm(npleioQTL[2], QTLmean, QTLcovmat[[2]])
    )
}

cat("Check Point 1 \n")
# show correlation between pleiotropic QTL effects
if (pleio_percent != 0) {
    cat("cor(QTLeffc1,QTLeffc2) of pleiotropic QTL for category 1: ", cor(pleioQTLeff$QTLeffc1[, 1], pleioQTLeff$QTLeffc1[, 2]), "\n")
    cat("cor(QTLeffc1,QTLeffc2) of pleiotropic QTL for category 2: ", cor(pleioQTLeff$QTLeffc2[, 1], pleioQTLeff$QTLeffc2[, 2]), "\n")
}

if (pleio_percent < 1.0) {
    # non-pleiotropic QTL positions for category 1 and category 2
    nonpleioQTLposc1 <- setdiff(QTLpos[["QTLposc1"]], pleioQTLpos[["QTLposc1"]])
    nonpleioQTLposc2 <- setdiff(QTLpos[["QTLposc2"]], pleioQTLpos[["QTLposc2"]])
    nonpleioQTLpos <- list(
        QTLposc1 = list(T1 = sample(nonpleioQTLposc1, nnonpleioQTLT1[1], replace = FALSE)),
        QTLposc2 = list(T1 = sample(nonpleioQTLposc2, nnonpleioQTLT1[2], replace = FALSE))
    )
    # the other half affect trait 2
    nonpleioQTLpos$QTLposc1$T2 <- setdiff(nonpleioQTLposc1, nonpleioQTLpos$QTLposc1$T1)
    nonpleioQTLpos$QTLposc2$T2 <- setdiff(nonpleioQTLposc2, nonpleioQTLpos$QTLposc2$T1)

    # save non-pleiotropic QTL positions
    for (i in 1:nCategory) {
        for (j in 1:nTrait) {
            write.table(nonpleioQTLpos[[i]][j], sprintf("%s/nonpleioQTLposc%dt%d.txt", outPath, i, j), row.names = F, col.names = F, quote = F)
        }
    }
    # non-pleiotropic QTL effects for category 1 and category 2
    nonpleioQTLeff <- list(
        QTLeffc1 = list(
            T1 = rnorm(nnonpleioQTLT1[1], mean=0, sd=sqrt(QTLcovmat[[1]][1, 1])),
            T2 = rnorm(nnonpleioQTLT2[1], mean=0, sd=sqrt(QTLcovmat[[1]][2, 2]))
        ),
        QTLeffc2 = list(
            T1 = rnorm(nnonpleioQTLT1[2], mean=0, sd=sqrt(QTLcovmat[[2]][1, 1])),
            T2 = rnorm(nnonpleioQTLT2[2], mean=0, sd=sqrt(QTLcovmat[[2]][2, 2]))
        )
    )
}

cat("Check Point 2 \n")

##############################################
# 6. Simulate Breeding Values 
##############################################
totalQTL <- unique(c(QTLpos$QTLposc1, QTLpos$QTLposc2))
totalQTLFile <- paste0(sprintf("%s/totalQTL.txt", outPath))
write.table(totalQTL, totalQTLFile, row.names = F, col.names = F, quote = F)

## select bfile subset based on caucal snplist
plinkQTL <- paste0(
    plink19, " --bfile ", rawDataPath, rawBimFileSuffix,
    " --keep ", selected_ind_file,
    " --extract ", totalQTLFile,
    " --make-bed ",
    "--threads 100 ",
    " --out ", paste0(selected_ind_out, ".QTL")
)
system(plinkQTL)

# Then 4 columns corresponding to the QTLtotal, effects: T1c1, T1c2, T2c1, T2c2
# (a little tedious for non-overlapping category, but it's a general framework for overlapping category)
QTLeffects = data.frame(matrix(0, nrow = length(totalQTL), ncol = 4))
rownames(QTLeffects) = totalQTL
colnames(QTLeffects) = c("T1c1", "T1c2", "T2c1", "T2c2")

# T1c1
if (pleio_percent != 0) {
    QTLeffects[pleioQTLpos$QTLposc1, "T1c1"] = pleioQTLeff$QTLeffc1[, 1]
    # T1c2
    QTLeffects[pleioQTLpos$QTLposc2, "T1c2"] = pleioQTLeff$QTLeffc2[, 1]
    # T2c1
    QTLeffects[pleioQTLpos$QTLposc1, "T2c1"] = pleioQTLeff$QTLeffc1[, 2]
    # T2c2
    QTLeffects[pleioQTLpos$QTLposc2, "T2c2"] = pleioQTLeff$QTLeffc2[, 2]
} 
if (pleio_percent < 1.0) {
    # T1c1
    QTLeffects[nonpleioQTLpos$QTLposc1$T1, "T1c1"] = nonpleioQTLeff$QTLeffc1$T1
    # T1c2
    QTLeffects[nonpleioQTLpos$QTLposc2$T1, "T1c2"] = nonpleioQTLeff$QTLeffc2$T1
    # T2c1
    QTLeffects[nonpleioQTLpos$QTLposc1$T2, "T2c1"] = nonpleioQTLeff$QTLeffc1$T2
    # T2c2
    QTLeffects[nonpleioQTLpos$QTLposc2$T2, "T2c2"] = nonpleioQTLeff$QTLeffc2$T2
}
QTLeffects$SNP = rownames(QTLeffects)

# Use gcta to calculate the breeding values
catTraitcombination = c("T1c1", "T1c2", "T2c1", "T2c2")
for (c in catTraitcombination) {
    print(c)
    # generate the eff file for gcta
    effFile = paste0(sprintf("%s/effFile_%s.txt", outPath, c))
    fwrite(QTLeffects[, c("SNP", c)], file = effFile, col.names = FALSE, sep = "\t")
}

### calcualte breeding values (geno*beta) for "T1c1", "T1c2", "T2c1", "T2c2"
for (c in catTraitcombination){
    print(c)
    gctaMlm <- paste0(
        gcta64, " --bfile ", paste0(selected_ind_out, ".QTL"),
        " --simu-qt ",
        " --simu-causal-loci ", sprintf("%s/effFile_%s.txt", outPath, c),
        " --simu-hsq ", 1,
        " --simu-rep ", repNum,
        " --thread-num 100",
        " --out ", sprintf("%s/bv4%s.gcta",outPath, c)
    )
    system(gctaMlm)
}
cat("Check Point 3 \n")

##############################################
#Delete bfiles once it is not needed
##############################################
# Create a pattern to match the file types
pattern <- "\\.bed$|\\.bim$|\\.fam$"

# List all files in the directory that match the pattern
files_to_delete <- list.files(path = outPath, pattern = pattern, full.names = TRUE)

cat("delete large bfiles: \n")
print(files_to_delete)

# Delete the files
file.remove(files_to_delete)
##############################################
# Delete bfiles once it is not needed
##############################################

# Standardize breeding values
bvs = data.frame(matrix(NA, nrow = sample_size, ncol = 4))
colnames(bvs) = catTraitcombination
for (c in catTraitcombination) {
    print(c)
    # read the breeding value file
    bv_df = fread(paste0(sprintf("%s/bv4%s.gcta", outPath, c), ".phen"))
    bvs[, c] = bv_df$V3
}
bvs$FID = bv_df$V1
bvs$IID = bv_df$V2

cat(sprintf("cor(bv1,bv2) for category 1: %f\n", cor(bvs$T1c1, bvs$T2c1)))
cat(sprintf("cor(bv1,bv2) for category 2: %f\n", cor(bvs$T1c2, bvs$T2c2)))

gcov_c1 = matrix(c(var(bvs$T1c1), cov(bvs$T1c1, bvs$T2c1), cov(bvs$T1c1, bvs$T2c1), var(bvs$T2c1)), nrow = 2, ncol = 2) 
write.table(gcov_c1, sprintf("%s/bvc1_cov_var.txt", outPath), row.names = F, col.names = F, quote = F)

gcov_c2 = matrix(c(var(bvs$T1c2), cov(bvs$T1c2, bvs$T2c2), cov(bvs$T1c2, bvs$T2c2), var(bvs$T2c2)), nrow = 2, ncol = 2) 
write.table(gcov_c2, sprintf("%s/bvc2_cov_var.txt", outPath), row.names = F, col.names = F, quote = F)

# total breeding values for trait 1
bv1 = bvs$T1c1 + bvs$T1c2
# total breeding values for trait 2
bv2 = bvs$T2c1 + bvs$T2c2
# sdbv <- c(sd(bv1), sd(bv2))
# Standardize breeding value
# bv1 <- bv1 / sd(bv1) * sqrt(h21)
# bv2 <- bv2 / sd(bv2) * sqrt(h22)
cat("Without standardization:\n")
cat("var(g1) = ", var(bv1), "\n")
cat("var(g2) = ", var(bv2), "\n")
cat("cov(bv1,bv2) = ", cov(bv1, bv2), "\n")
cat("cor(bv1,bv2) = ", cor(bv1, bv2), "\n")

# save breeding value
#write.table(bv1, sprintf("%s/bv1.txt", outPath), row.names = F, col.names = F, quote = F)
#write.table(bv2, sprintf("%s/bv2.txt", outPath), row.names = F, col.names = F, quote = F)
gcov=matrix(c(var(bv1), cov(bv1, bv2), cov(bv1, bv2), var(bv2)), nrow = 2, ncol = 2)
write.table(gcov, sprintf("%s/bv_cov_var.txt", outPath), row.names = F, col.names = F, quote = F)


##############################################
# 7. Residual correlation
##############################################
res_mean <- c(0, 0) # mean of residual effect for trait 1 and trait 2
res_11 = (1 - h21)
res_22 = (1 - h22)
res_12 = sqrt(res_11 * res_22) * Rescor
res_covmat <- matrix(c(res_11, res_12, res_12, res_22), nrow = 2, ncol = 2, byrow = TRUE) # residual covariance matrix
eAll <- mvrnorm(sample_size, res_mean, res_covmat) # residual effect
e1 <- eAll[, 1] # residual effect for trait 1
e2 <- eAll[, 2] # residual effect for trait 2
# e1 <- e1 / sd(e1) * sqrt(1 - h21)
# e2 <- e2 / sd(e2) * sqrt(1 - h22)
cat("Without standardization:\n")
cat("var(e1) = ", var(e1), "\n")
cat("var(e2) = ", var(e2), "\n")
cat("cov(e1,e2) = ", cov(e1, e2), "\n")
cat("cor(e1,e2) = ", cor(e1, e2), "\n")

# Standardize phenotypes
y1 <- bv1 + e1
y2 <- bv2 + e2
# sdy <- c(sd(y1), sd(y2))
# y1 <- (y1 - mean(y1)) / sd(y1)
# y2 <- (y2 - mean(y2)) / sd(y2)

cat("Without standardization:\n")
cat("var(y1) = ", var(y1), "\n")
cat("var(y2) = ", var(y2), "\n")
cat("cov(y1,y2) = ", cov(y1, y2), "\n")
cat("cor(y1,y2) = ", cor(y1, y2), "\n")

# cat("final varg1 =", var((bv1 - mean(y1)) / sd(y1)), "\n")
# cat("final varg2 =", var((bv2 - mean(y1)) / sd(y1)), "\n")

# save phenotypes
write.table(y1, sprintf("%s/y1.txt", outPath), row.names = F, col.names = F, quote = F)
write.table(y2, sprintf("%s/y2.txt", outPath), row.names = F, col.names = F, quote = F)

##############################################
##############################################
# save QTL effects(after standardization)
##############################################
##############################################
# if (pleio_percent != 0) {
#     #step1. find pleo QTL
#     qtl_c1 = QTLeffects[pleioQTLpos[[1]] , c("T1c1", "T2c1")]
#     #step2. scale
#     qtl_c1[,"T1c1"] = qtl_c1[,"T1c1"] / sdbv[1] * sqrt(h21) / sdy[1]
#     qtl_c1[,"T2c1"] = qtl_c1[,"T2c1"] / sdbv[2] * sqrt(h22) / sdy[2]
#     #step3. calculate cov
#     var_c1t1=var(qtl_c1[,"T1c1"])
#     var_c1t2=var(qtl_c1[,"T2c1"])
#     cov_c1t12=cov(qtl_c1[,"T1c1"], qtl_c1[,"T2c1"])
#     cov_c1 = c(var_c1t1, cov_c1t12, cov_c1t12, var_c1t2)
#     write.table(qtl_c1, sprintf("%s/scaled_pleioQTLeff_c1.txt", outPath), row.names = F, col.names = F, quote = F)
#     write.table(cov_c1, sprintf("%s/scaled_pleioQTLeff_cov_c1.txt", outPath), row.names = F, col.names = F, quote = F)

#     #step1. find pleo QTL
#     qtl_c2=QTLeffects[pleioQTLpos[[2]] , c("T1c2", "T2c2")]
#     #step2. scale
#     qtl_c2[,"T1c2"] = qtl_c2[,"T1c2"] / sdbv[1] * sqrt(h21) / sdy[1]
#     qtl_c2[,"T2c2"] = qtl_c2[,"T2c2"] / sdbv[2] * sqrt(h22) / sdy[2]
#     #step3. calculate cov
#     var_c2t1=var(qtl_c2[,"T1c2"])
#     var_c2t2=var(qtl_c2[,"T2c2"])
#     cov_c2t12=cov(qtl_c2[,"T1c2"], qtl_c2[,"T2c2"])
#     cov_c2 = c(var_c2t1, cov_c2t12, cov_c2t12, var_c2t2)
#     write.table(qtl_c2, sprintf("%s/scaled_pleioQTLeff_c2.txt", outPath), row.names = F, col.names = F, quote = F)
#     write.table(cov_c2, sprintf("%s/scaled_pleioQTLeff_cov_c2.txt", outPath), row.names = F, col.names = F, quote = F)
# }

#save scale
# write.table(sdbv, sprintf("%s/sdbv.txt", outPath), row.names = F, col.names = F, quote = F)
# write.table(sdy, sprintf("%s/sdy.txt", outPath), row.names = F, col.names = F, quote = F)


# save plink phenotypes 
phentrait1File <- paste0(outPath,"Trait1.gcta.phen.plink")
fwrite(data.frame(FID = bvs$FID, IID = bvs$IID, y = y1), file = phentrait1File, col.names = FALSE, sep = "\t")

phentrait2File <- paste0(outPath, "Trait2.gcta.phen.plink")
fwrite(data.frame(FID = bvs$FID, IID = bvs$IID, y = y2), file = phentrait2File, col.names = FALSE, sep = "\t")


##############################################
# 8. Generate summary data
##############################################
# Trait 1 
plinkGWAS <- paste0(
    plink19, " --bfile ", rawDataPath, rawBimFileSuffix,
    " --keep ", selected_ind_file,
    " --pheno ", phentrait1File,
    " --linear ",
    " --ci 0.95 ",
    " --out ", phentrait1File, ".ci"
)
system(plinkGWAS)

# Trait 2
plinkGWAS <- paste0(
    plink19, " --bfile ", rawDataPath, rawBimFileSuffix,
    " --keep ", selected_ind_file,
    " --pheno ", phentrait2File,
    " --linear ",
    " --ci 0.95 ",
    " --out ", phentrait2File, ".ci"
)
system(plinkGWAS)

cat("Check Point 4 \n")
##############################################
## Re-format result into ma format

headerPlink <- c("SNP", "A1", "A2", "MAF", "BETA", "SE", "P", "N")
headerMa <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
freqFile <- sprintf("%s.frq", selected_ind_out)

bivTraitFile1 <- paste0(phentrait1File, ".ci", ".assoc.linear")
bivTraitFile2 <- paste0(phentrait2File, ".ci", ".assoc.linear")

bivTrait1 <- fread(bivTraitFile1)
bivTrait2 <- fread(bivTraitFile2)
freq <- fread(freqFile)


bivTraitAll1 <- merge(bivTrait1, freq, by = c("SNP", "A1", "CHR"))
bivTraitAll2 <- merge(bivTrait2, freq, by = c("SNP", "A1", "CHR"))

bivTraitAll1$N <- sample_size
bivTraitAll2$N <- sample_size

bivTraitMa1 <- bivTraitAll1[, ..headerPlink]
bivTraitMa2 <- bivTraitAll2[, ..headerPlink]

names(bivTraitMa1) <- headerMa
names(bivTraitMa2) <- headerMa

fwrite(bivTraitMa1,
    file = paste0(bivTraitFile1, ".ma"),
    col.names = TRUE, sep = "\t"
)

fwrite(bivTraitMa2,
    file = paste0(bivTraitFile2, ".ma"),
    col.names = TRUE, sep = "\t"
)

##############################################
# delete unused files (to save memory)
##############################################

# delete other larger files
files_to_delete2 = paste0(outPath, c(
    "bv4T1c1.gcta.phen", "bv4T2c1.gcta.phen", "bv4T1c2.gcta.phen", "bv4T2c2.gcta.phen",
    "Trait1.gcta.phen.plink.ci.assoc.linear", "Trait2.gcta.phen.plink.ci.assoc.linear",
    "y1.txt", "y2.txt", "Trait1.gcta.phen.plink", "Trait2.gcta.phen.plink"
))

file.remove(files_to_delete2)
file.remove(selected_ind_file)
file.remove(freqFile)




print("Done!")


########################################## Test##########################################

# is_positive_definite <- function(mat) {
#     eigenvalues <- eigen(mat, only.values = TRUE)$values
#     all(eigenvalues > 0)
# }

# compute_correlation <- function(covmat) {
#     d12 <- covmat[1, 2]
#     d12 / sqrt(covmat[1, 1] * covmat[2, 2])
# }

# pleio_percent = 1
# proportion_h2_explained = c(0.5, 0.5)
# pis = c(0.01, 0.01)

# h21 = 0.5
# h22 = 0.2

# nSNPc <- c(8633, 87149)
# nQTL <- floor(nSNPc * pis)

# npleioQTL <- floor(nQTL * pleio_percent) # number of pleiotropic QTL in each annotation group
# cat("number of pleiotropic QTL in each annotation group: ", paste(npleioQTL, collapse = ", "), "\n")
# if (pleio_percent < 1.0) {
#     nnonpleioQTL <- nQTL - npleioQTL # number of non-pleiotropic QTL in each annotation group
#     cat("number of non-pleiotropic QTL in each annotation group: ", paste(nnonpleioQTL, collapse = ", "), "\n")
# } else {
#     nnonpleioQTL <- c(0, 0)
# }

# # only half affect trait 1
# nnonpleioQTLT1 <- floor(nnonpleioQTL / 2)
# # only half affect trait 2
# nnonpleioQTLT2 <- nnonpleioQTL - nnonpleioQTLT1

# QTLcovmat <- list()
# QTL_cor <- c(0.8, 0.2)
# genetic_cor <- QTL_cor * pleio_percent /(pleio_percent + (1 - pleio_percent) / 2)

# for (c in 1:2) {
#     genetic_var_11 <- h21 * proportion_h2_explained[c]
#     genetic_var_22 <- h22 * proportion_h2_explained[c]
#     genetic_cov_12 <- sqrt(genetic_var_11 * genetic_var_22) * genetic_cor[c]
#     nQTL11 <- npleioQTL[c] + nnonpleioQTLT1[c]
#     nQTL22 <- npleioQTL[c] + nnonpleioQTLT2[c]
#     nQTL12 <- nQTL21 <- npleioQTL[c]
#     QTLcovmat[[c]] <- matrix(c(genetic_var_11 / nQTL11, genetic_cov_12 / nQTL12, genetic_cov_12 / nQTL21, genetic_var_22 / nQTL22), nrow = 2, ncol = 2)
#     cat("QTL covariance matrix for category ", c, ": \n")
#     print(QTLcovmat[[c]])
#     cat("The QTL correlation is: ", compute_correlation(QTLcovmat[[c]]), "\n")
#     cat("The QTL covariance is positive or not?", is_positive_definite(QTLcovmat[[c]]), "\n")

#     cat("genetic covariance matrix for category ", c, ": \n")
#     genetic_var = matrix(c(genetic_var_11, genetic_cov_12, genetic_cov_12, genetic_var_22), nrow = 2, ncol = 2)
#     print(genetic_var)
#     cat("The genetic correlation is: ", compute_correlation(genetic_var), "\n")
# }




