#!/bin/bash

parent_dir="/common/zhao/jyqqu/MTSBayesCC/data/real_data/GWAS_summary_statistics/ImputedSumstat/"

module load anaconda
conda activate gnova

for file in "$parent_dir"/*.imputed.ma; do
    # Get the base name without extension
    base_name=$(basename "$file")

    # Define expected output file
    out_file="${file}.munge.sumstats"

    # Skip if output already exists
    if [[ -f "$out_file" ]]; then
        echo "Skipping $base_name — already munged."
        continue
    fi

    echo "Processing: $base_name"

    /common/zhao/jyqqu/GNOVA/munge_sumstats.py \
        --sumstats "$file" \
        --out "${file}.munge"
done
