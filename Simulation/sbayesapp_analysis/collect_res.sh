#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4500
#SBATCH --time=99:00:00
#SBATCH --job-name=collect_res
#SBATCH --partition=zhao,batch
#SBATCH --mail-user=jyqqu@ucdavis.edu
#SBATCH --mail-type=FAIL

module load julia
julia --project=/common/zhao/jyqqu/mtsbayescc collect_res_total_gcor_h2.jl
julia --project=/common/zhao/jyqqu/mtsbayescc collect_res_anno_stratified_gcor_h2.jl
