#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4500
#SBATCH --time=00:50:00
#SBATCH --job-name=LDSC
#SBATCH --partition=zhao,batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jyqqu@ucdavis.edu

plink_folder="/common/zhao/jyqqu/MTSBayesCC/data/bfiles/sldsc/"

module load anaconda
conda activate ldsc

for pleio in 0.1 0.5 0.8 1 
do
        echo "pleio $pleio"
        for samSize in 300000
        do
                echo "samSize $samSize"
                for hsq in "0.5-0.2"
                do
                        echo "hsq $hsq"
                        h21=$(echo "$hsq" | cut -d'-' -f1)
                        h22=$(echo "$hsq" | cut -d'-' -f2)
                        for rep in {2..20}
                        do
                                echo "rep $rep"
                                analysis_folder="h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${pleio}_sampleSize${samSize}_seed${rep}"
                                annot_folder="/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/${analysis_folder}/sldsc_annot/"

                                
                                ## Step 2: Computing LD scores with an annot file
                                echo "Computing LD scores with the annot file SimData.1.annot.gz"
                                python /home/zhao/jyqqu/ldsc/ldsc.py \
                                --l2 \
                                --bfile ${plink_folder}/1000G.EUR.QC.merged.1 \
                                --ld-wind-cm 1 \
                                --annot ${annot_folder}/SimData.1.annot.gz \
                                --thin-annot \
                                --out ${annot_folder}/SimData.1
                                
                        done
                done
        done
done
