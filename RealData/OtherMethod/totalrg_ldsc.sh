#!/bin/bash
#mkdir -p $ldsc_dir/ldsc/totalrg_gctb_impute/

gwas_dir="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/GWAS_summary_statistics/ImputedSumstat/"
ldsc_dir="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/ldsc/"
# List of directories you want to process
trait1="SCZ"
trait2="IQ"

#gwas_files=("AF.imputed.ma" "CHOL.imputed.ma" "FG.imputed.convertsign.ma" "Height.imputed.ma" "PD.imputed.ma")
#traits=("AD" "AF" "BMI" "CHOL" "DEP" "EA" "FG" "Height" "IBD" "INS" "PC" "PD" "RA" "RBC" "SBP" "SCZ")

module load anaconda
conda activate ldsc

# Step 2: Run munge_sumstats.py with formatted b1 and b2 files
/home/zhao/jyqqu/ldsc/munge_sumstats.py \
--sumstats $gwas_dir/${trait1}.imputed.ma \
--out $ldsc_dir/totalrg_gctb_impute/${trait1}

/home/zhao/jyqqu/ldsc/munge_sumstats.py \
--sumstats $gwas_dir/${trait2}.imputed.ma \
--out $ldsc_dir/totalrg_gctb_impute/${trait2}

# Step 3: Run ldsc.py with the formatted files
/home/zhao/jyqqu/ldsc/ldsc.py \
--rg $ldsc_dir/totalrg_gctb_impute/${trait1}.sumstats.gz,$ldsc_dir/totalrg_gctb_impute/${trait2}.sumstats.gz \
--ref-ld-chr /mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/ldsc_eur_w_ld_chr/ \
--w-ld-chr /mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/ldsc_eur_w_ld_chr/ \
--out $ldsc_dir/totalrg_gctb_impute/${trait1}_${trait2}

echo "Finished processing trait1 $trait1 and trait2 $trait2"

conda deactivate
