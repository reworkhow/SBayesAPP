#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=00:20:00
#SBATCH --job-name=SLDSC
#SBATCH --partition=zhao,batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jyqqu@ucdavis.edu

module load anaconda
conda activate ldsc

h21="0.5"
h22="0.2"
pleio="0.1"
samSize="300000"
rep="1"

plink_folder="/common/zhao/jyqqu/MTSBayesCC/data/bfiles/"
analysis_folder="h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${pleio}_sampleSize${samSize}_seed${rep}"
annot_folder="/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/${analysis_folder}/sldsc_annot/"
cd $annot_folder

python /home/zhao/jyqqu/ldsc/munge_sumstats.py \
--sumstats ../Trait1.gcta.phen.plink.ci.assoc.linear.ma \
--out ./trait1

python /home/zhao/jyqqu/ldsc/ldsc.py \
	--h2 ./trait1.sumstats.gz \
	--ref-ld-chr ./SimData.1. \
	--w-ld-chr $plink_folder/weights.hm3_noMHC.1. \
	--out ./ldsc_example_trait1
