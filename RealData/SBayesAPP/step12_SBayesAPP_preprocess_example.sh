#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=24:00:00
#SBATCH --job-name=sbapp_preprocess
#SBATCH --partition=batch,guest
#SBATCH --output=/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/SBayesAPP_step12_test.log
#SBATCH --mail-user=jyqqu@ucdavis.edu
#SBATCH --mail-type=FAIL


export OPENBLAS_NUM_THREADS=1
source ~/.bashrc
export JULIA_LOAD_PATH=/mnt/nrdstor/zhao/jyqqu/mtsbayescc:

# Define path to script
script_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/codes/SBayesAPP_preprocess.jl"

# Run the Julia script with arguments
julia $script_path \
  --LD_info_path "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/eigen_data_whole_genome/" \
  --out "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.0.5.h2_trait2.0.2_pleioPercent0.5_sampleSize300000_seed1/" \
  --trait1_file "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.0.5.h2_trait2.0.2_pleioPercent0.5_sampleSize300000_seed1/Trait1.gcta.phen.plink.ci.assoc.linear.ma" \
  --trait2_file "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.0.5.h2_trait2.0.2_pleioPercent0.5_sampleSize300000_seed1/Trait2.gcta.phen.plink.ci.assoc.linear.ma" \
  --LDinfo_file "snp.info" \
  --annot_file "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/example_analysis/annotation_df.txt" \
  --annot_dict_name "anno_matrix_dict" \
  --nrank 2 \
  --nblock 49
