Rcode="/common/zhao/jyqqu/MTSBayesCC/codes/annot.gz4sldsc.R"

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
                        for rep in {1..20}
                        do
                                echo "rep $rep"
                                #qsubName="h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${pleio}_sampleSize${samSize}_annotationSize0.5_seed${rep}"
                                qsubName="h2_trait1.${h21}.h2_trait2.${h22}_pleioPercent${pleio}_sampleSize${samSize}_seed${rep}"
				file="sbatch_files/annot_$qsubName.sh"
                                cat ST_header > $file
                                echo "#SBATCH -J ${qsubName}" >> $file
                                echo module load R >> $file
                                echo Rscript $Rcode --seed $rep --sample_size $samSize --pleio_percent $pleio --h21 ${h21} --h22 ${h22} >> $file
                                chmod u+x $file
                                sbatch $file
                        done
                done
        done
done
