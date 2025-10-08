#submit MTSBayesCC job
###########################
scriptPath="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/analysis/real_data/"
logPath="${scriptPath}/log_folder_step3_estcorrR_tuned_v3_estSigma/"
mkdir -p ${logPath}

scriptName="step3_MTSBayesCC_SNPorder_shuffleTraitOrder_estGscale_estcorrR_tuned_v3_estSigma.jl"

niter=3000
outfreq=50
thin=50
nrank=20
chainlength="3K"
# Generate 30 random unique integers between 1 and 100000
seeds=($(shuf -i 1-100000 -n 10))

#annot_file="../cell_type_annot_human_total162.filtered.with_rest.annot"
annot_file="../brain.TDEP_0kb_SBayesRC.wRest.txt"

#annot_dict="anno_matrix_cell_type_human_total_unoverlap_dict"
#annot_dict="anno_matrix_cell_type_annot_human_total162_wRest_dict"
annot_dict="anno_matrix_brain_wRest_0kb_dict"


ST_folder="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/analysis/real_data/STSBayesC_vareThres/"

fixed_hyperparameters="false"
is_continue="false"
estimate_pi="true"

# read in the sample size N statistics for the data
csv_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/T2D_N_statistics.csv"  # <-- adjust if needed
# Build an associative array: MEAN["T2D_AF"]=1030840 (etc.)
declare -A MEAN
while IFS=, read -r Folder Std Mean Ratio; do
  [[ "$Folder" == "Folder" || -z "$Folder" ]] && continue
  # strip possible quotes
  Folder=${Folder//\"/}
  Mean=${Mean//\"/}
  MEAN["$Folder"]="$Mean"
done < "$csv_path"

#trait_pairs=("T2D:AF" "T2D:BMI" "T2D:CHOL" "T2D:DEP" "T2D:EA" "T2D:Height" "T2D:IBD" "T2D:INS" "T2D:PC" "T2D:PD" "T2D:RA" "T2D:RBC" "T2D:SBP" "T2D:SCZ")
trait_pairs=("SCZ:EA")
#N1=933970
# N1=85716 # LC
# N2=325709 # CigDay
N1=130571 # SCZ
# N2=265741 # IQ
# N1=573184 #T2D_small
# N2=88972 #FG
N2=766345 # EA


for pair in "${trait_pairs[@]}"; do
	echo "Processing pair: $pair"
	IFS=":" read -r trait1 trait2 <<< "$pair"
	folder_name="${trait1}_${trait2}"
	preprocess_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/$folder_name/"
	echo "${seeds[@]}" > "${preprocess_path}/seeds_brain_wRest_0kb2.txt"

	for seed in "${seeds[@]}"; do
		echo "Seed: $seed"
		starting_value_dir="XXX"
		secondary_starting_value_dir="XXX"
		ST_path="$ST_folder/$folder_name/"
		
		# v3 -> estGscale computed from beta'beta/m
		analysis_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/analysis/real_data/brainAnnot_0kb/MTSBayesCC_${chainlength}_corrR_tuned_v3_estSigma/${folder_name}_seed${seed}/"
		mkdir -p $analysis_path

		LogName="SBayesAPP_${folder_name}_seed${seed}_${chainlength}_brainAnnot_0kb"
		printLog="${logPath}/$LogName.log"
		jobFile="${logPath}/$LogName.sh"
		{
			echo "#!/bin/bash" 
			echo "#SBATCH --nodes=1" 
			echo "#SBATCH --ntasks-per-node=$nrank" 
			echo "#SBATCH --mem-per-cpu=15G" 
			echo "#SBATCH --time=48:00:00"  # Adjust time as necessary
			echo "#SBATCH --job-name=SBayesAPP_${folder_name}_seed${seed}_brainAnnot_0kb"
			echo "#SBATCH --partition=batch,guest"  # Adjust partition as necessary
			echo "#SBATCH --output=${printLog}" 
			echo "#SBATCH --mail-user=jyqqu@ucdavis.edu" 
			echo "#SBATCH --mail-type=FAIL" 

			echo "export OPENBLAS_NUM_THREADS=1" 
			echo "source ~/.bashrc" 
			echo "export JULIA_LOAD_PATH=/mnt/nrdstor/zhao/jyqqu/mtsbayescc:" 
			echo "unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU" 

			echo "srun -n $nrank --cpus-per-task=1 --cpu-bind=cores \\" 
			echo "julia --project=/mnt/nrdstor/zhao/jyqqu/mtsbayescc \\" 
			echo "    ${scriptPath}/${scriptName} \\" 
			echo "    $preprocess_path \\" 
			echo "    $analysis_path \\" 
			echo "    $niter \\" 
			echo "    $seed \\" 
			echo "    $nrank \\" 
			echo "    $annot_file \\" 
			echo "    $annot_dict \\" 
			echo "    $outfreq \\" 
			echo "    $starting_value_dir \\" 
			echo "    $secondary_starting_value_dir \\" 
			echo "    $ST_path \\" 
			echo "    $thin \\" 
			echo "    $N1 \\" 
			echo "    $N2 \\" 
			echo "    $estimate_pi \\" 
			echo "    $fixed_hyperparameters \\" 
			echo "    $is_continue \\" 
			echo "    $chr > ${printLog}"
		} > "$jobFile"

		sbatch $jobFile  
		#done
	done		
done
