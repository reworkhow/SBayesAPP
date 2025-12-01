#!/bin/bash
#SBATCH --job-name=hdl_batch
#SBATCH --output=hdl_batch_log.out
#SBATCH --error=hdl_batch_log.err
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=2
#SBATCH --partition=batch,guest

module load R

base_dir="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/"
data_dir="${base_dir}/GWAS_summary_statistics/ImputedSumstat/"
hdl_script="/mnt/nrdstor/zhao/jyqqu/HDL/HDL.data.wrangling.R"
LD_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/UKB_imputed_SVD_eigen99_extraction"
hdl_result_dir="${base_dir}/HDL_results"
sumstat_log_dir="${hdl_result_dir}/sumstat_format_log_dir"
analysis_log_dir="${hdl_result_dir}/analysis_log_dir"

mkdir -p "$sumstat_log_dir" "$analysis_log_dir"

trait1="LC"
trait2="CigDay"

out_dir="${hdl_result_dir}/${trait1}_${trait2}"
mkdir -p "$out_dir"

# Format Trait1
Rscript "$hdl_script" \
    gwas.file="${data_dir}/${trait1}.imputed.ma" \
    LD.path="$LD_path" \
    SNP="SNP" A1="A1" A2="A2" N="N" b="b" se="se" \
    output.file="${out_dir}/trait1" \
    log.file="${sumstat_log_dir}/trait1"

# Format Trait2
Rscript "$hdl_script" \
    gwas.file="${data_dir}/${trait2}.imputed.ma" \
    LD.path="$LD_path" \
    SNP="SNP" A1="A1" A2="A2" N="N" b="b" se="se" \
    output.file="${out_dir}/trait2" \
    log.file="${sumstat_log_dir}/trait2"

# Run HDL + Z-score calc and output analysis log
Rscript - <<EOF > "${analysis_log_dir}/trait1_trait2.log" 2>&1
library(HDL)
trait1_path <- "${out_dir}/trait1.hdl.rds"
trait2_path <- "${out_dir}/trait2.hdl.rds"

gwas1 <- readRDS(trait1_path)
gwas2 <- readRDS(trait2_path)
gwas1\$Z <- gwas1\$b / gwas1\$se
gwas2\$Z <- gwas2\$b / gwas2\$se
saveRDS(gwas1, trait1_path)
saveRDS(gwas2, trait2_path)

res.HDL <- HDL.rg.parallel(gwas1, gwas2, "${LD_path}", numCores = 2)
saveRDS(res.HDL, file = "${out_dir}/res.HDL.rds")

est_df <- as.data.frame(res.HDL\$estimates.df)
write.csv(est_df, file = "${out_dir}/hdl_estimates.csv", row.names = TRUE)

print(sprintf("Finished folder: %s", "${out_dir}"))
EOF

