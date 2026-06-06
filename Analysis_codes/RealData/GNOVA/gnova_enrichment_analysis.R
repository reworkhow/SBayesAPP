library(data.table)
library(ggplot2)
library(dplyr)
#library(ggtext)

# Base paths
path <- "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/GNOVA_RealDataAnalysis/"
data_dir <- "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/"

# Read SNP counts
annot_file <- paste0(data_dir, "cell_type_annot_human_total162.filtered.with_rest.annot")
annot_df <- fread(annot_file)
nLoci <- colSums(annot_df[, -1])
nLoci_df <- data.frame(annot_name = names(nLoci), nLoci = as.numeric(nLoci))
ntotalsnp <- sum(nLoci)

# Read relatedness file once
#relatedness_to_T2D <- fread(sprintf("%s/relatedness_to_T2D.csv", data_dir), header = TRUE)[, 1:2]
#relatedness <- fread(sprintf("%s/relatedness_to_lung_humantotal.csv", data_dir), header = TRUE)[,c(1,3)]

# Enrichment function
compute_enrichment_se_df <- function(gnova_gcov_df, ntotal_snp) {
  gcov <- gnova_gcov_df$rho_corrected
  se_gcov <- gnova_gcov_df$rho_corrected_se
  nloci <- gnova_gcov_df$nLoci
  #S <- sum(gcov) # Sum of all genetic covariance estimates
  S <- abs(sum(gcov)) # Sum of all genetic covariance estimates -> absolute values

  enrichment <- numeric(nrow(gnova_gcov_df))
  se_enrichment <- numeric(nrow(gnova_gcov_df))

  for (index_c in 1:nrow(gnova_gcov_df)) {
    Rc <- gcov[index_c] / S
    var_Rc_c <- ((1 / S - gcov[index_c] / S^2)^2 * se_gcov[index_c]^2)
    var_Rc_others <- sum(((gcov[-index_c] / S^2)^2) * se_gcov[-index_c]^2)
    var_Rc <- var_Rc_c + var_Rc_others
    se_Rc <- sqrt(var_Rc)
    se_enrichment[index_c] <- (ntotal_snp / nloci[index_c]) * se_Rc
  }

  gnova_gcov_df$enrichment_se <- se_enrichment
  gnova_gcov_df$enrichment <- (gcov / S) * (ntotal_snp / gnova_gcov_df$nLoci)
  return(gnova_gcov_df)
}
# Create output directory if not exist
plot_dir <- file.path(path, "plots", "human_total162")
dir.create(plot_dir, showWarnings = FALSE)

# Get all GNOVA files
gnova_files <- list.files(path, pattern = "\\.cell_type_annot_human_total162.filtered.with_rest.txt$", full.names = TRUE)

# Storage objects
gcor_total_df <- data.frame(
  trait_pair = character(),
  gcor_total = numeric(),
  stringsAsFactors = FALSE
)

all_gnova_gcov_enrichment <- list()


# Loop over files
for (file in gnova_files) {
  trait_pair <- sub("\\.cell_type_annot_human_total162\\.filtered\\.with_rest\\.txt$", "", basename(file))
  message("Processing trait pairs: ", trait_pair)

  gnova_df <- fread(file)

  sum_h21 = sum(gnova_df$h2_1)
  sum_h22 = sum(gnova_df$h2_2)
  sum_rho = sum(gnova_df$rho_corrected)
  total_gcor = sum_rho/sqrt(sum_h21*sum_h22)

  # save gcor_total
  gcor_total_df <- rbind(
    gcor_total_df,
    data.frame(trait_pair = trait_pair, gcor_total = total_gcor, stringsAsFactors = FALSE)
  )
  
  gnova_gcov_df <- gnova_df[, c("annot_name", "rho_corrected", "pvalue_corrected")]

  # z-score and SE
  gnova_gcov_df$zscore_corrected <- qnorm(1 - (gnova_gcov_df$pvalue_corrected / 2))
  gnova_gcov_df$rho_corrected_se <- gnova_gcov_df$rho_corrected / gnova_gcov_df$zscore_corrected

  # Merge SNP counts
  nLoci_df$annot_name <- chartr("-", ".",   nLoci_df$annot_name)
  gnova_gcov_df <- merge(nLoci_df, gnova_gcov_df, by = "annot_name")
  names(gnova_gcov_df)[1] <- "Annotation"

  # Compute enrichment
  gnova_gcov_df <- compute_enrichment_se_df(gnova_gcov_df, ntotalsnp)

  # Merge relatedness
  gnova_gcov_df$Annotation <- gsub("\\.bed", "", gnova_gcov_df$Annotation)
  #gnova_gcov_enrichment <- merge(relatedness_to_T2D, gnova_gcov_df, by = "Annotation")
  #gnova_gcov_enrichment <- merge(relatedness, gnova_gcov_df, by = "Annotation")
  gnova_gcov_enrichment <- gnova_gcov_df

  # Order
  gnova_gcov_enrichment <- gnova_gcov_enrichment[order(gnova_gcov_enrichment$enrichment, decreasing = TRUE), ]
  orderedAnnotation <- gnova_gcov_enrichment$Annotation
  # order Annotation by orderedAnnotation
  gnova_gcov_enrichment$Annotation <- factor(gnova_gcov_enrichment$Annotation, levels = orderedAnnotation)

  # Add relatedness
  # Confidence intervals
  gnova_gcov_enrichment$lb <- gnova_gcov_enrichment$enrichment - gnova_gcov_enrichment$enrichment_se
  gnova_gcov_enrichment$ub <- gnova_gcov_enrichment$enrichment + gnova_gcov_enrichment$enrichment_se

  gnova_gcov_enrichment <- gnova_gcov_enrichment %>%
      mutate(Enrichment_by_1SD = case_when(
          enrichment > 0 & lb > 1 ~ "Enrichment",
          enrichment > 0 & ub < 1 ~ "Depletion",
          enrichment < 0 & ub < -1 ~ "Enrichment",
          enrichment < 0 & lb > -1 ~ "Depletion",
          TRUE ~ "Non-significant"
      ))
  
  # save the enrichment data
  write.table(gnova_gcov_enrichment, file = sprintf("%s/gnova_enrichment_%s.txt", plot_dir, trait_pair), sep = "\t", row.names = FALSE, quote = FALSE)

  # tag trait pair
  gnova_gcov_enrichment$trait_pair <- trait_pair

  # store
  all_gnova_gcov_enrichment[[trait_pair]] <- gnova_gcov_enrichment

  # # Color annotations
  # gnova_gcov_enrichment$Annotation_colored <- paste0(
  #   "<span style='color:",
  #   dplyr::recode(gnova_gcov_enrichment$Relatedness,
  #     "Pancreatic" = "#457a4d",
  #     "Direct"     = "#4DAF4A",
  #     "Indirect"   = "#43CD80",
  #     "No"         = "#928C8C"
  #   ),
  #   "'>", gnova_gcov_enrichment$Annotation, "</span>"
  # )
  # gnova_gcov_enrichment$Annotation_colored <- factor(
  #   gnova_gcov_enrichment$Annotation_colored,
  #   levels = gnova_gcov_enrichment$Annotation_colored[match(orderedAnnotation, gnova_gcov_enrichment$Annotation)]
  # )

  # # Plot
  # p1 <- ggplot(gnova_gcov_enrichment, aes(x = Annotation, y = enrichment, fill = Annotation)) +
  #   geom_bar(stat = "identity") +
  #   geom_errorbar(aes(ymin = enrichment - enrichment_se, ymax = enrichment + enrichment_se), width = 0.2) +
  #   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  #   geom_hline(yintercept = -1, linetype = "dashed", color = "black") +
  #   # scale_fill_manual(
  #   #   values = c(
  #   #     "Pancreatic" = "#457a4d",
  #   #     "Direct"     = "#4DAF4A",
  #   #     "Indirect"   = "#43CD80",
  #   #     "No"         = "#928C8C"
  #   #   )
  #   # ) +
  #   labs(
  #     title = sprintf("Coh2 Enrichment (GNOVA, corrected) - %s", trait_pair),
  #     y = "Enrichment",
  #     x = "Annotation"
  #   ) +
  #   theme(
  #     axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 6),
  #     legend.position = "bottom",
  #     legend.text = element_text(size = 5)
  #   )

  # # Save
  # ggsave(sprintf("%s/gnova_rho_enrichment_%s_lungAnnot_0kb_test.png", plot_dir, trait_pair), plot = p1, width = 10, height = 10, dpi = 300)
}
gcor_total_df
# Combine all enrichment dfs
all_gnova_gcov_enrichment_df <- do.call(rbind, all_gnova_gcov_enrichment)
head(all_gnova_gcov_enrichment_df)
dim(all_gnova_gcov_enrichment_df) 
# save gnova enrichment
write.table(all_gnova_gcov_enrichment_df, file = sprintf("%s/gnova_enrichment_all_trait_pairs.txt", plot_dir), sep = "\t", row.names = FALSE, quote = FALSE)
# save gcor_total
write.table(gcor_total_df, file = sprintf("%s/gnova_gcor_total_all_trait_pairs.txt", plot_dir), sep = "\t", row.names = FALSE, quote = FALSE)

# get CI
gnova_gcov_enrichment = fread(sprintf("%s/plots/gnova_enrichment_SCZ_EA_brainAnnot_0kb.txt", path))
gnova_gcov_enrichment <- fread("/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/GNOVA_RealDataAnalysis/plots/human_total162/gnova_enrichment_T2D_FG.txt")
head(gnova_gcov_enrichment)

gnova_gcov_enrichment$CIlb_90 <- gnova_gcov_enrichment$enrichment - qnorm(0.95) * gnova_gcov_enrichment$enrichment_se
gnova_gcov_enrichment$CIub_90 <- gnova_gcov_enrichment$enrichment + qnorm(0.95) * gnova_gcov_enrichment$enrichment_se
gnova_gcov_enrichment[gnova_gcov_enrichment$CIlb_90 > 1,]

## check enrichment zscore for BMI (156 annotations)
gnova_df = fread("/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/GNOVA_RealDataAnalysis/plots/human_total162/gnova_enrichment_T2D_BMI.txt")
gnova_origin = fread("/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/GNOVA_RealDataAnalysis/SCZ_EA.brain.TDEP_0kb_SBayesRC.wRest.txt")
gnova_origin$gcor <- gnova_origin$rho_corrected/sqrt(gnova_origin$h2_1*gnova_origin$h2_2)
# order gnova_origin by pvalue_corrected
head(gnova_origin, 20)
write.csv(gnova_origin, file = "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/GNOVA_RealDataAnalysis/SCZ_EA.brain.TDEP_0kb_SBayesRC.wRest.addGcor.txt", row.names = FALSE)

gnova_df$enrichment_zscore <- gnova_df$enrichment / gnova_df$enrichment_se
# order by abs zscore
gnova_df <- gnova_df[order(abs(gnova_df$pvalue_corrected)), ]
head(gnova_df, 20)

# compute h2 enrichment
file = sprintf("%s/SCZ_EA.brain.TDEP_0kb_SBayesRC.wRest.txt",path)
trait_pair <- sub("\\.brain.TDEP_0kb_SBayesRC.wRest.txt$", "", basename(file))
message("Processing trait pairs: ", trait_pair)

compute_h2_enrichment_df <- function(gnova_h2_df, ntotal_snp) {
  h21 <- gnova_h2_df$h2_1
  h22 <- gnova_h2_df$h2_2
  nloci <- gnova_h2_df$nLoci
  S_h21 <- sum(h21)
  S_h22 <- sum(h22)

  gnova_h2_df$enrichment_h21 <- (h21 / S_h21) * (ntotal_snp / gnova_h2_df$nLoci)
  gnova_h2_df$enrichment_h22 <- (h22 / S_h22) * (ntotal_snp / gnova_h2_df$nLoci)
  return(gnova_h2_df)
}

gnova_df <- fread(file)
gnova_h2_df <- gnova_df[, c("annot_name", "h2_1", "h2_2")]

# Merge SNP counts
nLoci_df$annot_name <- chartr("-", ".",   nLoci_df$annot_name)
gnova_h2_df <- merge(nLoci_df, gnova_h2_df, by = "annot_name")
names(gnova_h2_df)[1] <- "Annotation"

# Compute enrichment
gnova_h2_df <- compute_h2_enrichment_df(gnova_h2_df, ntotalsnp)
write.table(gnova_h2_df, file = sprintf("%s/gnova_h2_enrichment_%s_brainAnnot_0kb.txt", plot_dir, trait_pair), sep = "\t", row.names = FALSE, quote = FALSE)

library(data.table)
sbayesapp_df = fread("/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/analysis/real_data/brainAnnot_0kb/MTSBayesCC_3K_corrR_tuned_v3_estSigma/corBtwChains_SCZ_EA/Group1_10chains_enr/gcor_gcov_results_Group1_SCZ_EA.csv")
range(sbayesapp_df$gcor)
sbayesapp_df_EA = sbayesapp_df[sbayesapp_df$Trait == "EA",]
sbayesapp_df_EA[sbayesapp_df_EA$Significance == "Enrichment",]
sbayesapp_df_EA[sbayesapp_df_EA$Mean < -1,]
