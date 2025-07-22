#submit MTSBayesCC job
###########################
scriptPath="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/analysis/real_data/"
logPath="${scriptPath}/log_folder_step3_estcorrR_tuned_meanR4blks_v3_gwas_snp_pthres_Height/"
mkdir -p ${logPath}

# scriptName="step3_MTSBayesCC_SNPorder_shuffleTraitOrder_estGscaleWholeChain.jl"
scriptName="step3_MTSBayesCC_SNPorder_shuffleTraitOrder_estGscale_estcorrR_tuned_meanR4blks_v3.jl"

niter=3000
outfreq=50
thin=50
nrank=20

# annot_file="../cell_type_annot_human_total_unoverlap.txt"
# annot_dict="anno_matrix_cell_type_human_total_unoverlap_dict"

annot_file="../gwas_snp_pthres_5e-5_Height.txt"
annot_dict="anno_matrix_gwas_snp_pthres_5e-5_Height_dict"

ST_folder="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/analysis/real_data/STSBayesC_vareThres/"

fixed_hyperparameters="false"

params=("Height")
#params=("BMI" "Height" "DEP" "EA" "INS" "AF" "IBD" "PC" "CHOL" "PD" "RA" "RBC" "SBP" "SCZ")
#params=("BMI" "Height" "CHOL")
chainlength="3K"
is_continue="false"

# chainlength="10K"
# previous_chainlength="$(( ${chainlength%K} - 10 ))K"

# if [ "$chainlength" == "10K" ]; then
#     is_continue="false"
# else
#     is_continue="true"
# fi

for param in "${params[@]}"; do
	echo "Trait: $param"
	for seed in 123
	do
		echo "Seed: $seed"
		starting_value_dir="XXX"
		secondary_starting_value_dir="XXX"
		ST_path="$ST_folder/T2D_$param/"
		preprocess_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/T2D_$param/"
		# v1 -> diagonal estGscale
		# v2 -> estGscale computed from alpha'alpha/m
		# v3 -> estGscale computed from beta'beta/m
		analysis_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/analysis/real_data/gwas_snp_pthres_Height/MTSBayesCC_${chainlength}_corrR_tuned_meanR4blks_v3/pthres_5e-5/T2D_${param}_seed${seed}/"
		mkdir -p $analysis_path

		#for chr in 22; do
        	#echo "Chr: $chr"	
		LogName="step3_MTSBayesCC_${param}_seed${seed}_${chainlength}_corrR_tuned_meanR4blks_v3_pthres_5e-5"
		printLog="${logPath}/$LogName.log"
		jobFile="${logPath}/$LogName.sh"
		echo "#!/bin/bash" > $jobFile
		echo "#SBATCH --nodes=1" >> $jobFile
		echo "#SBATCH --ntasks-per-node=$nrank" >> $jobFile
		echo "#SBATCH --mem-per-cpu=12G" >> $jobFile
		echo "#SBATCH --time=24:00:00" >> $jobFile  # Adjust time as necessary
		echo "#SBATCH --job-name=step3_MTSBayesCC_$param" >> $jobFile
		echo "#SBATCH --partition=zhao,batch,guest" >> $jobFile  # Adjust partition as necessary
		echo "#SBATCH --output=${printLog}" >> $jobFile
		echo "#SBATCH --mail-user=jyqqu@ucdavis.edu" >> $jobFile
		echo "#SBATCH --mail-type=FAIL" >> $jobFile

		echo "export OPENBLAS_NUM_THREADS=1" >> $jobFile
		echo "source ~/.bashrc" >> $jobFile
		echo "export JULIA_LOAD_PATH=/mnt/nrdstor/zhao/jyqqu/mtsbayescc:" >> $jobFile

		echo "srun julia --project=/mnt/nrdstor/zhao/jyqqu/mtsbayescc \\" >> $jobFile
		echo "    ${scriptPath}/${scriptName} \\" >> $jobFile
		echo "    $preprocess_path \\" >> $jobFile
		echo "    $analysis_path \\" >> $jobFile
		echo "    $niter \\" >> $jobFile
		echo "    $seed \\" >> $jobFile
		echo "    $nrank \\" >> $jobFile
		echo "    $annot_file \\" >> $jobFile
		echo "    $annot_dict \\" >> $jobFile
		echo "    $outfreq \\" >> $jobFile
		echo "    $starting_value_dir \\" >> $jobFile
		echo "    $secondary_starting_value_dir \\" >> $jobFile
		echo "    $ST_path \\" >> $jobFile
		echo "    $thin \\" >> $jobFile
		echo "    $fixed_hyperparameters \\" >> $jobFile
		echo "    $is_continue \\" >> $jobFile
		echo "    $chr > ${printLog}" >> $jobFile
		
		sbatch $jobFile  
		#done
	done		
done
