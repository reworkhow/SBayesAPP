###########################
scriptPath="/common/zhao/jyqqu/MTSBayesCC/codes/"
logPath="${scriptPath}/../log_step1_real_data"
mkdir -p ${logPath}

# Define the base directory and parameters
#params=("AF" "BMI" "CHOL" "DEP" "EA" "Height" "IBD" "INS" "PC" "PD" "RA" "RBC" "SBP" "SCZ")
params=("AD")
Rcode="step1_standardize_pheno_var_GCTBLDinput.R"
gwas_path="/common/zhao/jyqqu/MTSBayesCC/data/real_data/GWAS_summary_statistics/ImputedSumstat/"
ldm_path="/common/zhao/jyqqu/MTSBayesCC/data/eigen_data_whole_genome/"

# Loop over the parameters and submit a job for each
for param in "${params[@]}"; do
	echo $param
	preprocess_path="/common/zhao/jyqqu/MTSBayesCC/data/real_data/T2D_$param/"
	mkdir -p $preprocess_path
	printLog="${logPath}/T2D_$param.log"
        jobFile="${logPath}/T2D_$param.sh"

	echo "#!/bin/bash" > $jobFile
	echo "#SBATCH --nodes=1" >> $jobFile
        echo "#SBATCH --ntasks-per-node=1" >> $jobFile
        echo "#SBATCH --cpus-per-task=2" >> $jobFile
        echo "#SBATCH --mem-per-cpu=4500" >> $jobFile
        echo "#SBATCH --time=24:00:00" >> $jobFile # Adjust time as necessary
	echo "#SBATCH --job-name=step1_$param" >> $jobFile
	echo "#SBATCH --partition=zhao,batch" >> $jobFile # Adjust partition as necessary
        echo "#SBATCH --output=${printLog}" >> $jobFile
        echo "#SBATCH --mail-user=jyqqu@ucdavis.edu" >> $jobFile
        echo "#SBATCH --mail-type=FAIL" >> $jobFile
	echo "module load R" >> $jobFile
	echo "Rscript ${scriptPath}/${Rcode} $gwas_path $ldm_path $preprocess_path T2D.imputed.ma $param.imputed.ma snp.info > ${printLog}" >> $jobFile
	sbatch $jobFile
done
