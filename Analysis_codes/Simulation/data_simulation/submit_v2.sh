###########################
homePath="/home/uqjzeng1/wd/proj/TJ"
scriptPath="/home/uqjzeng1/wd/proj/TJ/sim_chr1_output_v2/"
logPath="${scriptPath}/log"
mkdir -p ${logPath}

scriptName="Sim2Tfunc_plink_chr1_memory_efficient_v2.R"
chmod 751 ${scriptPath}/${scriptName}


h2_pairs=("0.5 0.2")

for h2 in "${h2_pairs[@]}"
do
    read h21 h22 <<< "$h2"

    for samSize in 300000
    do
            for plePer in 0.1 0.5 0.8 1
            do
                for seed in {1..20}
                do
                    qsubName="h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${plePer}_sampleSize${samSize}_seed${seed}"
                    printLog="${logPath}/${qsubName}.log"
                    jobFile="${logPath}/${qsubName}.sh"

                    echo "#!/bin/bash" > $jobFile
                    echo "#SBATCH --nodes=1" >> $jobFile
                    echo "#SBATCH --ntasks-per-node=1" >> $jobFile
                    echo "#SBATCH --cpus-per-task=1" >> $jobFile
                    echo "#SBATCH --mem=10G" >> $jobFile
                    echo "#SBATCH --time=48:00:00" >> $jobFile  # Adjust time as necessary
                    echo "#SBATCH --job-name=${qsubName}" >> $jobFile
                    echo "#SBATCH --partition=general" >> $jobFile  # Adjust partition as necessary
                    echo "#SBATCH --account=a_mcrae" >> $jobFile  # Adjust account as necessary
                    echo "#SBATCH --output=${printLog}" >> $jobFile

                    echo "module load r" >> $jobFile
                    
                    echo "/sw/auto/rocky8c/epyc3/software/R/4.3.3-gfbf-2023a/bin/Rscript ${scriptPath}/${scriptName} \
                          --seed ${seed} \
                          --sample_size ${samSize} \
                          --h21 ${h21} \
                          --h22 ${h22} \
                          --pleio_percent ${plePer} > ${printLog}" >> $jobFile

                    sbatch $jobFile
                done
            done
        
    done
done

