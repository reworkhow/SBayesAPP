#!/bin/bash
###########################
scriptPath="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/codes/"
preprocessScript="${scriptPath}/SBayesAPP_preprocess.jl"

logDir="${scriptPath}/../log_step12_real_data"
mkdir -p "$logDir"

trait_pairs=("T2D_small:FG")

# Where the imputed .ma inputs are stored
imputed_dir="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/GWAS_summary_statistics/ImputedSumstat"

# Fixed annotation file + dict
annot_file="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/cell_type_annot_human_total_unoverlap.txt"
annot_dict="anno_matrix_cell_type_human_total_unoverlap_dict"

# LD info
ldm_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/eigen_data_whole_genome/"
LDinfo_file="snp.info"

# Model params
nrank=20
nblock=591

# Julia env
JULIA_ENV="/mnt/nrdstor/zhao/jyqqu/mtsbayescc"

# ======= loop over pairs =======
for pair in "${trait_pairs[@]}"; do
    IFS=":" read -r trait1 trait2 <<< "$pair"

	echo "trait 1: $trait1, trait 2: $trait2"

    # outputs live here
    out_dir="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/${trait1}_${trait2}"
    mkdir -p "${out_dir}"

    # inputs from ImputedSumstat
    trait1_file="${imputed_dir}/${trait1}.imputed.ma"
    trait2_file="${imputed_dir}/${trait2}.imputed.convertsign.ma"

    # check files exist
    missing=0
    [[ -f "$trait1_file" ]] || { echo "Missing: $trait1_file"; missing=1; }
    [[ -f "$trait2_file" ]] || { echo "Missing: $trait2_file"; missing=1; }
    if [[ $missing -eq 1 ]]; then
        echo "Skip ${trait1}:${trait2} due to missing inputs."
        continue
    fi

    jobName="sbpp_${trait1}_${trait2}"
    logFile="${logDir}/SBayesAPP_preprocess_${trait1}_${trait2}.log"
    jobFile="${logDir}/submit_${trait1}_${trait2}.sbatch"

    cat > "$jobFile" <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --job-name=${jobName}
#SBATCH --partition=batch,guest
#SBATCH --output=${logFile}
#SBATCH --mail-user=jyqqu@ucdavis.edu
#SBATCH --mail-type=FAIL

export OPENBLAS_NUM_THREADS=1
source ~/.bashrc
export JULIA_LOAD_PATH=${JULIA_ENV}:
    
julia "${preprocessScript}" \\
  --LD_info_path "${ldm_path}" \\
  --out "${out_dir}" \\
  --trait1_file "${trait1_file}" \\
  --trait2_file "${trait2_file}" \\
  --LDinfo_file "${LDinfo_file}" \\
  --annot_file "${annot_file}" \\
  --annot_dict_name "${annot_dict}" \\
  --nrank ${nrank} \\
  --nblock ${nblock}
EOF
    echo "Submitting ${jobName}"
    sbatch "$jobFile"
done
