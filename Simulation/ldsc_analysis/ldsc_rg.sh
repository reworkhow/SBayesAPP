#!/bin/bash

h2_pairs=("0.5 0.2")

for h2 in "${h2_pairs[@]}"
do
    read h21 h22 <<< "$h2"

    # Loop over different conditions
    for samSize in 300000
    do
        # for annoSize in 0.5
        # do
            for plePer in 0.1 0.5 0.8 1
            do
                for seed in {1..20}
                do
                    # Generate the output directory path
                    output_dir="/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${plePer}_sampleSize${samSize}_seed${seed}"
                    cd $output_dir
                    mkdir -p ./ldsc/totalrg/

                    # Generate the batch script for the current configuration
                    cat << EOF > ldsc_${samSize}_${plePer}_${seed}.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4500
#SBATCH --time=24:00:00
#SBATCH --job-name=ldsc_${samSize}_${plePer}_${seed}
#SBATCH --partition=zhao,batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jyqqu@ucdavis.edu

module load anaconda
conda activate ldsc

echo "Processing in directory: $output_dir"

/home/zhao/jyqqu/ldsc/munge_sumstats.py \\
--sumstats Trait1.gcta.phen.plink.ci.assoc.linear.ma  \\
--out ./ldsc/totalrg/trait1

/home/zhao/jyqqu/ldsc/munge_sumstats.py \\
--sumstats Trait2.gcta.phen.plink.ci.assoc.linear.ma \\
--out ./ldsc/totalrg/trait2

/home/zhao/jyqqu/ldsc/ldsc.py \\
--rg ./ldsc/totalrg/trait1.sumstats.gz,./ldsc/totalrg/trait2.sumstats.gz \\
--ref-ld-chr /common/zhao/tianjing/software/eur_w_ld_chr/ \\
--w-ld-chr /common/zhao/tianjing/software/eur_w_ld_chr/ \\
--out ./ldsc/totalrg/ldsc_res

EOF

                    # Submit the batch job
                    sbatch ldsc_${samSize}_${plePer}_${seed}.sh

                done
            done
        #done
    done
done
