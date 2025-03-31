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
                for seed in {2..20}
                do
                    # Generate the output directory path
                    output_dir="/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${plePer}_sampleSize${samSize}_seed${seed}"
                    cd $output_dir

                    # Generate the batch script for the current configuration
                    cat << EOF > gnova_${samSize}_${plePer}_${seed}.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4500
#SBATCH --time=24:00:00
#SBATCH --job-name=gnova_${samSize}_${plePer}_${seed}
#SBATCH --partition=zhao,batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jyqqu@ucdavis.edu

module load anaconda
conda activate gnova

echo "Processing in directory: $output_dir"

python /common/zhao/jyqqu/GNOVA/gnova.py \\
Trait1.gcta.phen.plink.ci.assoc.linear.ma.munge.sumstats \\
Trait2.gcta.phen.plink.ci.assoc.linear.ma.munge.sumstats \\
--bfile /common/zhao/jyqqu/MTSBayesCC/data/bfiles/1000G.EUR.QC.1.filtered \\
--annot anno.txt \\
--out gnova.txt
EOF

                    # Submit the batch job
                    sbatch gnova_${samSize}_${plePer}_${seed}.sh

                done
            done
        #done
    done
done



