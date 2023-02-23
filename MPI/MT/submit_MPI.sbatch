#!/bin/bash -l
#SBATCH --job-name=MPI
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2500
#SBATCH --time=00:10:00
#SBATCH --partition=high
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu

module load julia

srun julia MTBayesC-func-LDblocks.mpi.jl
