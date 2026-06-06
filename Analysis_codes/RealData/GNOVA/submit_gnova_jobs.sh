#!/bin/bash

# Trait list
traits=("AF" "CHOL" "DEP" "EA" "IBD" "INS" "PC" "PD" "RA" "RBC" "SBP" "SCZ")

# GNOVA sbatch script path
gnova_script="gnova_submit.sh"

# Loop through traits and submit jobs
for trait in "${traits[@]}"; do
    echo "Submitting GNOVA job for $trait"
    sbatch $gnova_script $trait
done

