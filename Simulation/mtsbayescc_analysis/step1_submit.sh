###########################
scriptPath="/common/zhao/jyqqu/MTSBayesCC/codes/"
logPath="${scriptPath}/../log_step1"
mkdir -p ${logPath}

scriptName="step1_standardize_pheno_var_Sim2Tfunc.R"

h2_pairs=("0.5 0.2")

for h2 in "${h2_pairs[@]}"
do
    read h21 h22 <<< "$h2"

    for samSize in 300000
    do
        #for annoSize in 0.5
        #do
            for plePer in 0.1 0.5 0.8 1
            do
                for seed in {1..20}
                do
                    #qsubName="h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${plePer}_sampleSize${samSize}_annotationSize${annoSize}_seed${seed}"
                    qsubName="h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${plePer}_sampleSize${samSize}_seed${seed}"
		    printLog="${logPath}/${qsubName}.log"
                    jobFile="${logPath}/${qsubName}.sh"

                    echo "#!/bin/bash" > $jobFile
                    echo "#SBATCH --nodes=1" >> $jobFile
                    echo "#SBATCH --ntasks-per-node=1" >> $jobFile
                    echo "#SBATCH --cpus-per-task=2" >> $jobFile
                    echo "#SBATCH --mem-per-cpu=4500" >> $jobFile
                    echo "#SBATCH --time=24:00:00" >> $jobFile  # Adjust time as necessary
                    echo "#SBATCH --job-name=${qsubName}" >> $jobFile
                    echo "#SBATCH --partition=zhao,batch" >> $jobFile  # Adjust partition as necessary
                    echo "#SBATCH --output=${printLog}" >> $jobFile
                    echo "#SBATCH --mail-user=jyqqu@ucdavis.edu" >> $jobFile
                    echo "#SBATCH --mail-type=FAIL" >> $jobFile

		    echo "module load R" >> $jobFile
                    echo "Rscript ${scriptPath}/${scriptName} \
                          --seed ${seed} \
                          --sample_size ${samSize} \
                          --h21 ${h21} \
                          --h22 ${h22} \
                          --pleio_percent ${plePer}  > ${printLog}" >> $jobFile
			# --annotation_size ${annoSize} \
                    sbatch $jobFile
                done
            done
        #done
    done
done
