#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=24:00:00
#SBATCH --job-name=gnova
#SBATCH --partition=zhao,batch,guest
#SBATCH --chdir=/common/zhao/jyqqu/MTSBayesCC/data/real_data/
#SBATCH --mail-user=jyqqu@ucdavis.edu
#SBATCH --mail-type=FAIL

# Load necessary modules (if any specific software modules are required)
module load anaconda

# Activate the environment
conda activate gnova


dir="/common/zhao/jyqqu/MTSBayesCC/data/real_data/GWAS_summary_statistics/ImputedSumstat/"
trait=$1

echo "Processing Trait: $trait"

    # Define the paths to the files
    b1_file="${dir}/T2D.imputed.ma.munge.sumstats"
    b2_file="${dir}/${trait}.imputed.convertsign.ma.munge.sumstats"

    # Check if the files are compressed (.gz) and unzip them
    if [ -f "${b1_file}.gz" ]; then
        echo "Unzipping ${b1_file}.gz"
        gunzip "${b1_file}.gz"
    fi

    if [ -f "${b2_file}.gz" ]; then
        echo "Unzipping ${b2_file}.gz"
        gunzip "${b2_file}.gz"
    fi

    # Run GNOVA for each folder
    python /common/zhao/jyqqu/GNOVA/gnova.py \
    "$b1_file" \
    "$b2_file" \
    --bfile /common/zhao/jyqqu/MTSBayesCC/data/CELLECT/data/ldsc/1000G_EUR_Phase3_plink_1M/1000G.EUR.QC.merged.filtered \
    --annot /common/zhao/jyqqu/MTSBayesCC/data/real_data/cell_type_annot_human_total_oneAnnotperSNP_BetaCellPriori.nosnpinfo.reordered4gnova.txt \
    --out /common/zhao/jyqqu/MTSBayesCC/data/real_data/GNOVA_RealDataAnalysis/T2D_${trait}.cell_type_annot_human_total_oneAnnotperSNP_BetaCellPriori.maFile.txt
