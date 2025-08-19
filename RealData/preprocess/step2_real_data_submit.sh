###########################
scriptPath="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/codes/"
logPath="${scriptPath}/../log_step2_real_data"
mkdir -p ${logPath}

scriptName="step2_data_dictionary_GCTBLDinput.jl"
#params=("AD" "AF" "BMI" "CHOL" "DEP" "EA" "Height" "IBD" "INS" "PC" "PD" "RA" "RBC" "SBP" "SCZ")
params=("Height")

trait="Height"
nrank=20
ldm_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/eigen_data_whole_genome/" 
annot_file="../gwas_snp_pthres_5e-5_$trait.txt"
annot_dict="anno_matrix_gwas_snp_pthres_5e-5_${trait}_dict"

for param in "${params[@]}"; do
	echo $param
	data_path="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/T2D_$param/"
	printLog="${logPath}/T2D_$param.log"
        jobFile="${logPath}/T2D_$param.sh"

        echo "#!/bin/bash" > $jobFile
        echo "#SBATCH --nodes=1" >> $jobFile
        echo "#SBATCH --ntasks-per-node=1" >> $jobFile
        echo "#SBATCH --mem=80G" >> $jobFile
        echo "#SBATCH --time=24:00:00" >> $jobFile  # Adjust time as necessary
        echo "#SBATCH --job-name=step2_$param" >> $jobFile
        echo "#SBATCH --partition=zhao,batch,guest" >> $jobFile  # Adjust partition as necessary
        echo "#SBATCH --output=${printLog}" >> $jobFile
        echo "#SBATCH --mail-user=jyqqu@ucdavis.edu" >> $jobFile
        echo "#SBATCH --mail-type=FAIL" >> $jobFile
        
        echo "export OPENBLAS_NUM_THREADS=1" >> $jobFile
	echo "source ~/.bashrc" >> $jobFile
        echo "julia --project=/mnt/nrdstor/zhao/jyqqu/mtsbayescc ${scriptPath}/${scriptName} $data_path $nrank $annot_file $ldm_path $annot_dict > ${printLog}" >> $jobFile
        sbatch $jobFile
done
