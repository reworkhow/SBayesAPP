#!/bin/bash

ma_file_name=$1
short_name=$2
output="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/real_data/GWAS_summary_statistics/ImputedSumstat"

/home/zhao/jyqqu/gctb_2.5.2_Linux/gctb --ldm-eigen /mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/eigen_data_whole_genome/ \
--gwas-summary $output/../$ma_file_name \
--impute-summary \
--out $output/$short_name \
--thread 4
