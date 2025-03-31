#!/bin/bash

module load anaconda
# Activate the environment
conda activate gnova


h2_pairs=("0.5 0.2")

for h2 in "${h2_pairs[@]}"
do
    read h21 h22 <<< "$h2"
    for samSize in 300000;
    do
        # for annoSize in 0.2 0.5; 
        # do
            for plePer in 0.1 0.5 0.8 1 
            do
                for seed in {2..20}
                do
                    for trait in Trait1 Trait2 
                    do
                        # Define the base path for input and output
                        basePath=/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${plePer}_sampleSize${samSize}_seed${seed}/${trait}.gcta.phen.plink.ci.assoc.linear
                        echo $basePath
                        # Define the output path for the munge file
                        mungePath=${basePath}.ma.munge

                        # Run the Python script to generate the .gz file
                        /common/zhao/tianjing/GNOVA/munge_sumstats.py \
                        --sumstats ${basePath}.ma \
                        --out $mungePath

                        # Unzip the generated .gz file
                        gunzip ${mungePath}.sumstats.gz
                    done
                done
            done
        #done
    done
done
