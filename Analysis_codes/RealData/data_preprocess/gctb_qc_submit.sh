#!/bin/bash -l

# "GCST004748_buildGRCh37_LungCancer.tsv.ma"
# "GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.ma"
ma_file_names=("Xue_et_al_T2D_META_Nat_Commun_2018.ma")
# "LC"
# "CigDay"
short_names=("T2D_small")

# Get the number of traits
num_traits=${#ma_file_names[@]}

# Loop through each trait by index
for (( i=0; i<num_traits; i++ )); do
    ma_file_name=${ma_file_names[$i]}
    short_name=${short_names[$i]}
    
    # Define the job name dynamically based on the trait
    job_name="GCTB_${short_name}"
    
    # Submit the job
    sbatch <<EOT
#!/bin/bash 
#SBATCH --partition=batch,guest
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name=$job_name

bash gctb_qc.sh $ma_file_name $short_name
EOT

    echo "Submitted job for trait: $short_name"
done
