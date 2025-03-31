#submit MTSBayesCC job
###########################

source ~/.bashrc 

scriptPath="/common/zhao/jyqqu/MTSBayesCC/analysis/sim_chr1_output_v2/MTS/MTSBayesCC_SNPorder_ShuffleTraits_estGscale_STinput2/"
ST_folder="/common/zhao/jyqqu/MTSBayesCC/analysis/sim_chr1_output_v2/STS/STSBayesC/h2_0.1_pi0.9/"
logPath="${scriptPath}/log_folder/"
mkdir -p ${logPath}

scriptName="step3_MTSBayesCC_Sim2Tfunc_SNPorder_shuffleTraitOrder_estGscale_STinput.jl"

niter=1000
outfreq=10

nrank=2
nmarker=95782


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
                    # if [ ! -f "/common/zhao/tianjing/simu_ch1_res/mt_s_bayesc_c_priorsimu/${qsubName}/estGcor_total.txt" ]; then
                        echo "${qsubName}"
                        echo "#!/bin/bash" > $jobFile
                        echo "#SBATCH --nodes=1" >> $jobFile
                        echo "#SBATCH --ntasks-per-node=2" >> $jobFile
                        echo "#SBATCH --cpus-per-task=1" >> $jobFile
                        echo "#SBATCH --mem-per-cpu=4500" >> $jobFile
                        echo "#SBATCH --time=168:00:00" >> $jobFile  # Adjust time as necessary
                        echo "#SBATCH --job-name=${qsubName}" >> $jobFile
                        echo "#SBATCH --partition=zhao,batch" >> $jobFile  # Adjust partition as necessary
                        echo "#SBATCH --output=${printLog}" >> $jobFile
                        echo "#SBATCH --mail-user=jyqqu@ucdavis.edu" >> $jobFile
                        echo "#SBATCH --mail-type=FAIL" >> $jobFile

			echo "module load julia" >> $jobFile
                        echo "srun julia --project=/common/zhao/jyqqu/mtsbayescc ${scriptPath}/${scriptName} \
		      	      --seed ${seed} \
                              --sample_size ${samSize} \
                              --h21 ${h21} \
                              --h22 ${h22} \
                              --pleio_percent ${plePer} \
                              --niter ${niter} \
                              --nrank ${nrank} \
                              --nmarker ${nmarker} \
                              --outfreq ${outfreq} \
			      --analysis_folder ${scriptPath} --ST_folder ${ST_folder} > ${printLog}" >> $jobFile

			   #  --annotation_size ${annoSize} \
                        sbatch $jobFile
                    # fi
                done
            done
        #done
    done
done
