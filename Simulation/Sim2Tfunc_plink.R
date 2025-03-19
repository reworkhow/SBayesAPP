#Simulation with annotation
# 3.4.1. R: generate bivarte-annotation
# ch1 only

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
  make_option(c("--annotation_size"), action="store", default=0.1, type='numeric',
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
annotation_size <- opt$annotation_size
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
QTLcor <- c(0.8, 0.2) # QTL correlation bewteen two traits for the two category
Rescor <- 0.1
nTrait <- 2

repNum <- 1 # number of replicates in gcta


rawDataPath <- paste0("/home/uqjzeng1/wd/proj/TJ/chr1_real_genotype_ukb/")
outPath <- paste0("/home/uqjzeng1/wd/proj/TJ/simu_chr1_output_new/")
rawBimFileSuffix <- paste0("ukbV3_eur_unrel_info0.6_maf0.01_hwe1e10_hm3_chr1")
outPath <- paste0(outPath, "h2_trait1.",h21,".h2_trait2.",h22,"_pleioPercent", pleio_percent, "_sampleSize", sample_size, "_annotationSize", annotation_size,"_seed",seed, "/")

###############################################
### Above are parameters to be changed
###############################################
# output parameters
cat("seed: ", seed, "\n")
cat("sample size: ", sample_size, "\n")
cat("annotation size: ", annotation_size, "\n")
cat("pleiotropic percent: ", pleio_percent, "\n")
cat("h2 for trait 1: ", h21, "\n")
cat("h2 for trait 2: ", h22, "\n")
cat("number of annotation groups: ", nCategory, "\n")
cat("proportion of non-zero effects in each annotation group: ", paste(pis, collapse = ", "), "\n")
cat("QTL correlation in each annotation category (for pleiotropic QTLs): ", paste(QTLcor, collapse = ", "), "\n")
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
nInd <- nrow(fam) #348501
nMarker <- nrow(bim) #723021
cat("number of Inds in raw file: ", nInd, "\n")
cat("number of SNPs in raw file: ", nMarker, "\n")

# Generate random indices for selected individuals without replacement
random_ind_index <- sort(sample(nInd, sample_size, replace = FALSE))
# selected individuals
selected_ind <- fam[random_ind_index, c("FID", "IID")]

selected_ind_file <- paste0(outPath, rawBimFileSuffix, "-selectedInd.txt")
fwrite(data.frame(selected_ind), file = selected_ind_file, col.names = FALSE, sep = "\t")
selected_ind_out <- sprintf("%s%s.%dInd", outPath, rawBimFileSuffix, sample_size)

plinkInd <- paste0(
    plink19, " --bfile ", rawDataPath, rawBimFileSuffix,
    " --keep ", selected_ind_file,
    " --make-bed ",
    "--threads 100 ",
    " --out ", selected_ind_out
)
system(plinkInd)


##  Calculate allele frequency
plinkFreq <- paste0(
    plink19, " --bfile ", selected_ind_out,
    " --freq ",
    " --out ", selected_ind_out
)
system(plinkFreq)

##############################################
# 2. Annotation size (10%,50%) ⇒ #SNPs in each annotation group
# select SNPs into annotation group based on annotation size parameters
##############################################
# Calculate the number of SNP to select for the category 1
nSNPc1 <- ceiling(nMarker * annotation_size)
nSNPc2 <- nMarker - nSNPc1
# Generate random SNP indices without replacement
random_snp_c1_index <- sort(sample(nMarker, nSNPc1, replace = FALSE))
# Select the random SNPs for category 1 using the generated indices
SNPsc1 <- bim$SNP[random_snp_c1_index]
SNPsc2 <- setdiff(bim$SNP, SNPsc1)
cat("number of SNPs selected in category 1: ", nSNPc1, "\n")
cat("number of SNPs selected in category 2: ", nSNPc2, "\n")

# save selected SNPs
write.table(SNPsc1, paste0(outPath, "SNPc1", ".txt"), row.names = F, col.names = F, quote = F)
write.table(SNPsc2, paste0(outPath, "SNPc2", ".txt"), row.names = F, col.names = F, quote = F)

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
}

##############################################
# 5. QTL effect correlation
##############################################

QTLmean <- c(0, 0) # mean of QTL effect for trait 1 and trait 2
QTLcovmat <- list( # marker covariance matrix for category 1 and category 2 (for npleioQTL)
    QTLcovmatc1 = matrix(c(1, QTLcor[1], QTLcor[1], 1), nrow = 2, ncol = 2, byrow = TRUE),
    QTLcovmatc2 = matrix(c(1, QTLcor[2], QTLcor[2], 1), nrow = 2, ncol = 2, byrow = TRUE)
)

QTLpos <- list(
    QTLposc1 = sample(SNPsc1, nQTL[1], replace = FALSE), # QTL positions for category 1
    QTLposc2 = sample(SNPsc2, nQTL[2], replace = FALSE) # QTL positions for category 2
) # QTL positions for category 1 and category 2

# pleiotropic QTL positions for category 1 and category 2
if (pleio_percent < 1.0 && pleio_percent!=0) {
    pleioQTLpos <- list(
        QTLposc1 = sample(QTLpos[["QTLposc1"]], npleioQTL[1], replace = FALSE),
        QTLposc2 = sample(QTLpos[["QTLposc2"]], npleioQTL[2], replace = FALSE)
    ) # pleiotropic QTL positions for category 1 and category 2
} else if (pleio_percent == 1.0) {
    pleioQTLpos <- QTLpos
} else if (pleio_percent == 0) {
    pleioQTLpos <- list(
        QTLposc1 = 0,
        QTLposc2 = 0
    )
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
        QTLeffc1 = mvrnorm(npleioQTL[1], QTLmean, QTLcovmat$QTLcovmatc1),
        QTLeffc2 = mvrnorm(npleioQTL[2], QTLmean, QTLcovmat$QTLcovmatc2)
    )
}

cat("Check Point 1 \n")
# show correlation between pleiotropic QTL effects
if (pleio_percent != 0) {
    cat("cor(QTLeffc1,QTLeffc2) for category 1: ", cor(pleioQTLeff$QTLeffc1[, 1], pleioQTLeff$QTLeffc1[, 2]), "\n")
    cat("cor(QTLeffc1,QTLeffc2) for category 2: ", cor(pleioQTLeff$QTLeffc2[, 1], pleioQTLeff$QTLeffc2[, 2]), "\n")
}

if (pleio_percent < 1.0) {
    # non-pleiotropic QTL positions for category 1 and category 2
    nonpleioQTLposc1 <- setdiff(QTLpos[["QTLposc1"]], pleioQTLpos[["QTLposc1"]])
    nonpleioQTLposc2 <- setdiff(QTLpos[["QTLposc2"]], pleioQTLpos[["QTLposc2"]])
    # only half affect trait 1
    nnonpleioQTLT1 <- floor(nnonpleioQTL / 2)
    # only half affect trait 2
    nnonpleioQTLT2 <- nnonpleioQTL - nnonpleioQTLT1
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
            T1 = rnorm(nnonpleioQTLT1[1], 0, 1),
            T2 = rnorm(nnonpleioQTLT2[1], 0, 1)
        ),
        QTLeffc2 = list(
            T1 = rnorm(nnonpleioQTLT1[2], 0, 1),
            T2 = rnorm(nnonpleioQTLT2[2], 0, 1)
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
    plink19, " --bfile ", selected_ind_out,
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
write.table(gcov_c1, sprintf("%s/bvc1_var_cov.txt", outPath), row.names = F, col.names = F, quote = F)

gcov_c2 = matrix(c(var(bvs$T1c2), cov(bvs$T1c2, bvs$T2c2), cov(bvs$T1c2, bvs$T2c2), var(bvs$T2c2)), nrow = 2, ncol = 2) 
write.table(gcov_c2, sprintf("%s/bvc2_var_cov.txt", outPath), row.names = F, col.names = F, quote = F)

# total breeding values for trait 1
bv1 = bvs$T1c1 + bvs$T1c2
# total breeding values for trait 2
bv2 = bvs$T2c1 + bvs$T2c2
sdbv <- c(sd(bv1), sd(bv2))
# Standardize breeding value
bv1 <- bv1 / sd(bv1) * sqrt(h21)
bv2 <- bv2 / sd(bv2) * sqrt(h22)
cat("After standardization:\n")
cat("var(g1) = ", var(bv1), "\n")
cat("var(g2) = ", var(bv2), "\n")
cat("cov(bv1,bv2) = ", cov(bv1, bv2), "\n")
cat("cor(bv1,bv2) = ", cor(bv1, bv2), "\n")

# save breeding value
# write.table(bv1, sprintf("%s/bv1.txt", outPath), row.names = F, col.names = F, quote = F)
# write.table(bv2, sprintf("%s/bv2.txt", outPath), row.names = F, col.names = F, quote = F)
gcov=matrix(c(var(bv1), cov(bv1, bv2), cov(bv1, bv2), var(bv2)), nrow = 2, ncol = 2)
write.table(gcov, sprintf("%s/bv_var_cov.txt", outPath), row.names = F, col.names = F, quote = F)


##############################################
# 7. Residual correlation
##############################################
res_mean <- c(0, 0) # mean of residual effect for trait 1 and trait 2
res_covmat <- matrix(c(1, Rescor, Rescor, 1), nrow = 2, ncol = 2, byrow = TRUE) # residual covariance matrix
eAll <- mvrnorm(sample_size, res_mean, res_covmat) # residual effect
e1 <- eAll[, 1] # residual effect for trait 1
e2 <- eAll[, 2] # residual effect for trait 2
e1 <- e1 / sd(e1) * sqrt(1 - h21)
e2 <- e2 / sd(e2) * sqrt(1 - h22)
cat("After standardization:\n")
cat("var(e1) = ", var(e1), "\n")
cat("var(e2) = ", var(e2), "\n")
cat("cov(e1,e2) = ", cov(e1, e2), "\n")
cat("cor(e1,e2) = ", cor(e1, e2), "\n")

# Standardize phenotypes
y1 <- bv1 + e1
y2 <- bv2 + e2
sdy <- c(sd(y1), sd(y2))
y1 <- (y1 - mean(y1)) / sd(y1)
y2 <- (y2 - mean(y2)) / sd(y2)

cat("After standardization:\n")
cat("var(y1) = ", var(y1), "\n")
cat("var(y2) = ", var(y2), "\n")
cat("cov(y1,y2) = ", cov(y1, y2), "\n")
cat("cor(y1,y2) = ", cor(y1, y2), "\n")

cat("final varg1 =", var((bv1 - mean(y1)) / sd(y1)), "\n")
cat("final varg2 =", var((bv2 - mean(y1)) / sd(y1)), "\n")

# save phenotypes
write.table(y1, sprintf("%s/y1.txt", outPath), row.names = F, col.names = F, quote = F)
write.table(y2, sprintf("%s/y2.txt", outPath), row.names = F, col.names = F, quote = F)

##############################################
##############################################
# save QTL effects(after standardization)
##############################################
##############################################
if (pleio_percent != 0) {
    #step1. find pleo QTL
    qtl_c1 = QTLeffects[pleioQTLpos[[1]] , c("T1c1", "T2c1")]
    #step2. scale
    qtl_c1[,"T1c1"] = qtl_c1[,"T1c1"] / sdbv[1] * sqrt(h21) / sdy[1]
    qtl_c1[,"T2c1"] = qtl_c1[,"T2c1"] / sdbv[2] * sqrt(h22) / sdy[2]
    #step3. calculate cov
    var_c1t1=var(qtl_c1[,"T1c1"])
    var_c1t2=var(qtl_c1[,"T2c1"])
    cov_c1t12=cov(qtl_c1[,"T1c1"], qtl_c1[,"T2c1"])
    cov_c1 = c(var_c1t1, cov_c1t12, cov_c1t12, var_c1t2)
    write.table(qtl_c1, sprintf("%s/scaled_pleioQTLeff_c1.txt", outPath), row.names = F, col.names = F, quote = F)
    write.table(cov_c1, sprintf("%s/scaled_pleioQTLeff_cov_c1.txt", outPath), row.names = F, col.names = F, quote = F)

    #step1. find pleo QTL
    qtl_c2=QTLeffects[pleioQTLpos[[2]] , c("T1c2", "T2c2")]
    #step2. scale
    qtl_c2[,"T1c2"] = qtl_c2[,"T1c2"] / sdbv[1] * sqrt(h21) / sdy[1]
    qtl_c2[,"T2c2"] = qtl_c2[,"T2c2"] / sdbv[2] * sqrt(h22) / sdy[2]
    #step3. calculate cov
    var_c2t1=var(qtl_c2[,"T1c2"])
    var_c2t2=var(qtl_c2[,"T2c2"])
    cov_c2t12=cov(qtl_c2[,"T1c2"], qtl_c2[,"T2c2"])
    cov_c2 = c(var_c2t1, cov_c2t12, cov_c2t12, var_c2t2)
    write.table(qtl_c2, sprintf("%s/scaled_pleioQTLeff_c2.txt", outPath), row.names = F, col.names = F, quote = F)
    write.table(cov_c2, sprintf("%s/scaled_pleioQTLeff_cov_c2.txt", outPath), row.names = F, col.names = F, quote = F)
}

#save scale
write.table(sdbv, sprintf("%s/sdbv.txt", outPath), row.names = F, col.names = F, quote = F)
write.table(sdy, sprintf("%s/sdy.txt", outPath), row.names = F, col.names = F, quote = F)


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
    plink19, " --bfile ", selected_ind_out,
    " --pheno ", phentrait1File,
    " --linear ",
    " --ci 0.95 ",
    " --out ", phentrait1File, ".ci"
)
system(plinkGWAS)

# Trait 2
plinkGWAS <- paste0(
    plink19, " --bfile ", selected_ind_out,
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
#delete bfiles (to save memory)
##############################################
# Create a pattern to match the file types
pattern <-  "\\.bed$|\\.bim$|\\.fam$|\\.frq$"
# List all files in the directory that match the pattern
files_to_delete <- list.files(path = outPath, pattern = pattern, full.names = TRUE)

cat("delete large files: \n")
print(files_to_delete)

# Delete the files
file.remove(files_to_delete)

# delete other larger files
files_to_delete2 = paste0(outPath,c("bv4T1c1.gcta.phen","bv4T2c1.gcta.phen", "bv4T1c2.gcta.phen", "bv4T2c2.gcta.phen",
                                    "Trait1.gcta.phen.plink.ci.assoc.linear","Trait2.gcta.phen.plink.ci.assoc.linear",
                                    "y1.txt","y2.txt",
                                    "Trait1.gcta.phen.plink","Trait2.gcta.phen.plink",
                                    "ukbV3_eur_unrel_info0.6_maf0.01_hwe1e10_hm3_chr1-selectedInd.txt"))
file.remove(files_to_delete2)

print("Done!")