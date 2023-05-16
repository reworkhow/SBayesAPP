# input is the genotypic covariate matrix
library(data.table)
library(MASS)
# parameters
geno_path <- "/Users/apple/Desktop/JiayiQu/UCD_PhD/MTSBayesC/data/" # path to genotype matrix
data_path <- "/Users/apple/Desktop/JiayiQu/UCD_PhD/MTSBayesC/data/Sim2Tfunc/" # path to save simulation parameters
seed <- 123 # set seed for reproducibility
sample_size <- 503
annotation_size <- 0.1 # 10% of SNPs in category 1
pleio_percent <- 0.5
h21 <- 0.5 # h2 for trait 1
h22 <- 0.1 # h2 for trait 2
nCategory <- 2 # number of annotation groups
pis <- rep(0.1, nCategory) # proportion of non-zero effects in each annotation group
QTLcor <- c(0.8, 0.2) # QTL correlation bewteen two traits for the two category
Rescor <- 0.1
folder <- paste0(data_path, "pleioPercent", pleio_percent, "_sampleSize", sample_size, "_annotationSize", annotation_size,"/") # folder to save simulation results
nTrait = 2

# output parameters
message("sample size: ", sample_size)
message("annotation size: ", annotation_size)
message("pleiotropic percent: ", pleio_percent)
message("h2 for trait 1: ", h21)
message("h2 for trait 2: ", h22)
message("number of annotation groups: ", nCategory)
message("proportion of non-zero effects in each annotation group: ", paste(pis, collapse = ", "))
message("QTL correlation: ", paste(QTLcor, collapse = ", "))
message("Residual correlation: ", Rescor)
message("folder to save simulation results: ", folder)
dir.create(folder)

# read in genotype matrix
geno = data.frame(fread(paste0(geno_path, "geno_n503_p18588_realSnpID.QC.csv"), header = T, sep = ","))# the first column is ID, and the rest of the colums are gene content for each SNP
rownames(geno) = geno[, 1] # IDs
geno <- geno[, -1] # remove the first column
nInd <- dim(geno)[1] # total number of individuals
nMarker <- dim(geno)[2] # total number of SNPs

##############################################
# 1. Sample Size (10K or 100K) ⇒ #individuals
# select individuals based on sample size parameters
##############################################

set.seed(seed) 
# Generate random row indices without replacement
random_indices <- sort(sample(nInd, sample_size, replace = FALSE))
# Select the random rows using the generated indices
geno <- geno[random_indices, ]
# save selected individuals 
write.table(rownames(geno), paste0(folder, "sampledInd_seed", seed, ".txt"), row.names = F, col.names = F, quote = F)
message("number of individuals selected: ", dim(geno)[1])

##############################################
# 2. Annotation size (1%/10%) ⇒ #SNPs in each annotation group
# select SNPs into annotation group based on annotation size parameters
##############################################
# Calculate the number of SNP to select for the category 1
nSNPc1 <- ceiling(nMarker * annotation_size)
# Generate random column indices without replacement
random_indices <- sort(sample(nMarker, nSNPc1, replace = FALSE))
# Select the random columns for category 1 using the generated indices
# Select the remaining columns for category 2
# save in a list
geno_list <- list(geno_c1 = geno[, random_indices], geno_c2 = geno[, -random_indices])
SNPsc1 <- colnames(geno_list$geno_c1) # SNP names in category 1
SNPsc2 <- colnames(geno_list$geno_c2) # SNP names in category 2
# save selected SNPs
write.table(SNPsc1, paste0(folder, "SNPsc1_seed", seed, ".txt"), row.names = F, col.names = F, quote = F)
write.table(SNPsc2, paste0(folder, "SNPsc2_seed", seed, ".txt"), row.names = F, col.names = F, quote = F)
message("number of SNPs selected in category 1: ", length(SNPsc1))
message("number of SNPs selected in category 2: ", length(SNPsc2))

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
    nnonpleioQTL = nQTL - npleioQTL # number of non-pleiotropic QTL in each annotation group
    message("number of non-pleiotropic QTL in each annotation group: ", paste(nnonpleioQTL, collapse = ", "))
}

##############################################
# 5. QTL effect correlation
##############################################
# # variance of QTL effect for trait 1
# QTLvar1 <- h21 / round(0.5 * sum(nnonpleioQTL) + sum(npleioQTL))
# # variance of QTL effect for trait 2
# QTLvar2 <- h22 / round(0.5 * sum(nnonpleioQTL) + sum(npleioQTL))
# # covariance of QTL effect between trait 1 and trait 2 for the two categories
# QTLcov <- QTLcor * sqrt(QTLvar1) * sqrt(QTLvar2)

QTLmean <- c(0, 0) # mean of QTL effect for trait 1 and trait 2
QTLcovmat <- list( # marker covariance matrix for category 1 and category 2 (for npleioQTL)
    QTLcovmatc1 = matrix(c(1, QTLcor[1], QTLcor[1], 1), nrow = 2, ncol = 2, byrow = TRUE),
    QTLcovmatc2 = matrix(c(1, QTLcor[2], QTLcor[2], 1), nrow = 2, ncol = 2, byrow = TRUE)
)

QTLpos = list(
    QTLposc1 = sample(SNPsc1, nQTL[1], replace = FALSE), # QTL positions for category 1
    QTLposc2 = sample(SNPsc2, nQTL[2], replace = FALSE) # QTL positions for category 2
)# QTL positions for category 1 and category 2

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
for (i in 1:nCategory){
    write.table(pleioQTLpos[[i]], sprintf("%s/pleioQTLposc%d_seed%d.txt", folder, i, seed), row.names = F, col.names = F, quote = F)
}

#pleiotropic QTL effects for category 1 and category 2
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
            write.table(nonpleioQTLpos[[i]][j], sprintf("%s/nonpleioQTLposc%dt%d_seed%d.txt", folder, i, j, seed), row.names = F, col.names = F, quote = F)
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

# Simulate breeding value
# Trait 1
bv1 <- rep(0, sample_size)
bv2 <- rep(0, sample_size)
# pleiotropic QTL
for (i in 1:nCategory) {
    bv1 = bv1 + as.matrix(geno_list[[i]][, pleioQTLpos[[i]]]) %*% pleioQTLeff[[i]][, 1]
    if (pleio_percent < 1.0) {
        bv1 = bv1 + as.matrix(geno_list[[i]][, unlist(nonpleioQTLpos[[i]]["T1"])]) %*% unlist(nonpleioQTLeff[[i]]["T1"])
    }
    bv2 = bv2 + as.matrix(geno_list[[i]][, pleioQTLpos[[i]]]) %*% pleioQTLeff[[i]][, 2]
    if (pleio_percent < 1.0) {
        bv2 = bv2 + as.matrix(geno_list[[i]][, unlist(nonpleioQTLpos[[i]]["T2"])]) %*% unlist(nonpleioQTLeff[[i]]["T2"])
    }
    message(sprintf("cor(bv1,bv2) for %d category: %f", i, cor(bv1, bv2)))
}
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
write.table(bv1, sprintf("%s/bv1.txt", folder), row.names = F, col.names = F, quote = F)
write.table(bv2, sprintf("%s/bv2.txt", folder), row.names = F, col.names = F, quote = F)

##############################################
# 5. Residual correlation
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
write.table(y1, sprintf("%s/y1.txt", folder), row.names = F, col.names = F, quote = F)
write.table(y2, sprintf("%s/y2.txt", folder), row.names = F, col.names = F, quote = F)

# save QTL effects(after standardization)
for (t in 1:nTrait){
    for (c in 1:nCategory){
        if (t == 1){ # trait 1
            pleioQTLeff[[c]][, t] = (pleioQTLeff[[c]][, t] / sdbv[1] * sqrt(h21)) / sdy[1]
        }
        else{ # trait 2
            pleioQTLeff[[c]][, t] = (pleioQTLeff[[c]][, t] / sdbv[2] * sqrt(h22)) / sdy[2]
        }
        write.table(pleioQTLeff[[c]][, t], sprintf("%s/pleioQTLeffc%dt%d_seed%d.txt", folder, c, t, seed), row.names = F, col.names = F, quote = F)
    }
}

if (pleio_percent < 1.0) {
    for (t in 1:nTrait){
        for (c in 1:nCategory){
            if (t == 1){ # trait 1
                nonpleioQTLeff[[c]]["T1"][[1]] = (unlist(nonpleioQTLeff[[c]]["T1"]) / sdbv[1] * sqrt(h21)) / sdy[1]
            }
            else{ # trait 2
                nonpleioQTLeff[[c]]["T2"][[1]] = (unlist(nonpleioQTLeff[[c]]["T2"]) / sdbv[2] * sqrt(h22)) / sdy[2]
            }
            write.table(nonpleioQTLeff[[c]][t], sprintf("%s/nonpleioQTLeffc%dt%d_seed%d.txt", folder, c, t, seed), row.names = F, col.names = F, quote = F)
        }
    }
}


