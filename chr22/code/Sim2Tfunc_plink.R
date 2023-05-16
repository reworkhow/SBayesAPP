options(scipen = 999)
library(MASS)
library(data.table)

###############################################
### Below are parameters to be changed
###############################################

MainPath <- "/home/carol666/MTSBayesC/"
gcta64 <- "/home/carol666/packages/gcta"
plink19 <- "/home/carol666/packages/plink"

# parameters
seed <- 123 # set seed for reproducibility
sample_size <- 500
annotation_size <- 0.1 # 10% of SNPs in category 1
pleio_percent <- 1.0
h21 <- 0.5 # h2 for trait 1
h22 <- 0.1 # h2 for trait 2
nCategory <- 2 # number of annotation groups
pis <- rep(0.01, nCategory) # proportion of non-zero effects in each annotation group
QTLcor <- c(0.8, 0.2) # QTL correlation bewteen two traits for the two category
Rescor <- 0.1
nTrait <- 2

repNum <- 1 # number of replicates in gcta


rawDataPath <- paste0("/home/carol666/MTSBayesC/data/hapmapData/")
outPath <- paste0("/home/carol666/MTSBayesC/data/Sim2T/Sim2Tfunc/exeOutput/")
rawBimFileSuffix <- paste0("g1000_eur.chr22.hapmap3")

###############################################
### Above are parameters to be changed
###############################################

outPath <- paste0(outPath, "pleioPercent", pleio_percent, "_sampleSize", sample_size, "_annotationSize", annotation_size, "/")

# output parameters
message("seed: ", seed)
message("sample size: ", sample_size)
message("annotation size: ", annotation_size)
message("pleiotropic percent: ", pleio_percent)
message("h2 for trait 1: ", h21)
message("h2 for trait 2: ", h22)
message("number of annotation groups: ", nCategory)
message("proportion of non-zero effects in each annotation group: ", paste(pis, collapse = ", "))
message("QTL correlation: ", paste(QTLcor, collapse = ", "))
message("Residual correlation: ", Rescor)
message("folder to save simulation results: ", outPath)
dir.create(outPath)


rawBimFile <- paste0(rawDataPath, rawBimFileSuffix, ".bim")
rawFamFile <- paste0(rawDataPath, rawBimFileSuffix, ".fam")
bim <- fread(rawBimFile)
names(bim) <- c("CHR", "SNP", "POS", "CM", "A1", "A2")
fam <- fread(rawFamFile)
names(fam) <- c("FID", "IID", "PID", "MID","SEX", "PHENOTYPE")

set.seed(seed)
##############################################
# 1. Sample Size (10K or 100K) ⇒ #individuals
# select individuals based on sample size parameters
##############################################
nInd <- nrow(fam)
nMarker <- nrow(bim)
# Generate random indices for selected individuals without replacement
random_ind_index <- sort(sample(nInd, sample_size, replace = FALSE))
# selected individuals
selected_ind <- fam[random_ind_index, c("FID", "IID")]

selected_ind_file <- paste0(outPath, rawBimFileSuffix, "-selectedInd.txt")
fwrite(data.frame(selected_ind), file = selected_ind_file, col.names = FALSE, sep = "\t")
selected_ind_out <- sprintf("%s/%s.%dInd", outPath, rawBimFileSuffix, sample_size)

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
# 2. Annotation size (1%/10%) ⇒ #SNPs in each annotation group
# select SNPs into annotation group based on annotation size parameters
##############################################
# Calculate the number of SNP to select for the category 1
nSNPc1 <- ceiling(nMarker * annotation_size)
# Generate random SNP indices without replacement
random_snp_c1_index <- sort(sample(nMarker, nSNPc1, replace = FALSE))
# Select the random SNPs for category 1 using the generated indices
SNPsc1 <- bim$SNP[random_snp_c1_index]
SNPsc2 <- setdiff(bim$SNP, SNPsc1)
message("number of SNPs selected in category 1: ", length(SNPsc1))
message("number of SNPs selected in category 2: ", length(SNPsc2))

# save selected SNPs
write.table(SNPsc1, paste0(outPath, "SNPsc1", ".txt"), row.names = F, col.names = F, quote = F)
write.table(SNPsc2, paste0(outPath, "SNPsc2", ".txt"), row.names = F, col.names = F, quote = F)

##############################################
# 3. Pi ⇒ number of QTL in each annotation group
# select QTL into annotation group based on Pi parameters
##############################################
# number of QTLs in each category
nSNPc2 <- nMarker - nSNPc1
nSNPc <- c(nSNPc1, nSNPc2)
nQTL <- floor(nSNPc * pis)
message("number of QTL in each annotation group: ", paste(nQTL, collapse = ", "))

##############################################
# 4. %pleiotropic QTL for 2 traits
# select pleiotropic QTL based on %pleiotropic parameters
##############################################
npleioQTL <- floor(nQTL * pleio_percent) # number of pleiotropic QTL in each annotation group
message("number of pleiotropic QTL in each annotation group: ", paste(npleioQTL, collapse = ", "))
if (pleio_percent < 1.0) {
    nnonpleioQTL <- nQTL - npleioQTL # number of non-pleiotropic QTL in each annotation group
    message("number of non-pleiotropic QTL in each annotation group: ", paste(nnonpleioQTL, collapse = ", "))
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
if (pleio_percent < 1.0) {
    pleioQTLpos <- list(
        QTLposc1 = sample(QTLpos[["QTLposc1"]], npleioQTL[1], replace = FALSE),
        QTLposc2 = sample(QTLpos[["QTLposc2"]], npleioQTL[2], replace = FALSE)
    ) # pleiotropic QTL positions for category 1 and category 2
} else if (pleio_percent == 1.0) {
    pleioQTLpos <- QTLpos
}
# save pleiotropic QTL positions
for (i in 1:nCategory) {
    write.table(pleioQTLpos[[i]], sprintf("%s/pleioQTLposc%d.txt", outPath, i), row.names = F, col.names = F, quote = F)
}

# pleiotropic QTL effects for category 1 and category 2
pleioQTLeff <- list(
    QTLeffc1 = mvrnorm(npleioQTL[1], QTLmean, QTLcovmat$QTLcovmatc1),
    QTLeffc2 = mvrnorm(npleioQTL[2], QTLmean, QTLcovmat$QTLcovmatc2)
)

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
QTLeffects[pleioQTLpos$QTLposc1, "T1c1"] = pleioQTLeff$QTLeffc1[, 1]
# T1c2
QTLeffects[pleioQTLpos$QTLposc2, "T1c2"] = pleioQTLeff$QTLeffc2[, 1]
# T2c1
QTLeffects[pleioQTLpos$QTLposc1, "T2c1"] = pleioQTLeff$QTLeffc1[, 2]
# T2c2
QTLeffects[pleioQTLpos$QTLposc2, "T2c2"] = pleioQTLeff$QTLeffc2[, 2]
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

message(sprintf("cor(bv1,bv2) for category 1: %f", cor(bvs$T1c1, bvs$T2c1)))
message(sprintf("cor(bv1,bv2) for category 2: %f", cor(bvs$T1c2, bvs$T2c2)))

# total breeding values for trait 1
bv1 = bvs$T1c1 + bvs$T1c2
# total breeding values for trait 2
bv2 = bvs$T2c1 + bvs$T2c2
sdbv <- c(sd(bv1), sd(bv2))
# Standardize breeding value
bv1 <- bv1 / sd(bv1) * sqrt(h21)
bv2 <- bv2 / sd(bv2) * sqrt(h22)
message("After standardization:")
message("varg1 = ", var(bv1))
message("varg2 = ", var(bv2))
message("cov(bv1,bv2) = ", cov(bv1, bv2))
message("cor(bv1,bv2) = ", cor(bv1, bv2))

# save breeding value
write.table(bv1, sprintf("%s/bv1.txt", outPath), row.names = F, col.names = F, quote = F)
write.table(bv2, sprintf("%s/bv2.txt", outPath), row.names = F, col.names = F, quote = F)

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
message("After standardization:")
message("vare1 = ", var(e1))
message("vare2 = ", var(e2))
message("cov(e1,e2) = ", cov(e1, e2))
message("cor(e1,e2) = ", cor(e1, e2))

# Standardize phenotypes
y1 <- bv1 + e1
y2 <- bv2 + e2
sdy <- c(sd(y1), sd(y2))
y1 <- (y1 - mean(y1)) / sd(y1)
y2 <- (y2 - mean(y2)) / sd(y2)

message("After standardization:")
message("vary1 = ", var(y1))
message("vary2 = ", var(y2))
message("final varg1 =", var((bv1 - mean(y1)) / sd(y1)))
message("final varg2 =", var((bv2 - mean(y1)) / sd(y1)))

# save phenotypes
write.table(y1, sprintf("%s/y1.txt", outPath), row.names = F, col.names = F, quote = F)
write.table(y2, sprintf("%s/y2.txt", outPath), row.names = F, col.names = F, quote = F)

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

bivTraitAll1$N <- 300000
bivTraitAll2$N <- 300000

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

print("Done!")