#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ======= Input Arguments =======
DATA_PATH="$REPO_ROOT/example/SBayesAPP_input_chr1/"
ANALYSIS_PATH="$REPO_ROOT/example/SBayesAPP_res_chr1/"
mkdir -p "$ANALYSIS_PATH"
NITER=1000
SEED=42
ANNOT_FILE="annotation_df.txt"
ANNOT_DICT="anno_matrix_dict"
OUTFREQ=100
STARTING_VALUE_DIR="XXX" # only needed if using starting values or extend chain length
GSCALE_VALUE_DIR="XXX" # only needed if using fixed Gscale values during continuation
ST_PATH="$REPO_ROOT/example/ST_res/"
THIN=50
N1=300000
N2=300000
N_CON=0
IS_CONTINUE="false"  # or "true" if continuing from a previous run
ESTIMATE_VARE="true"
ESTIMATE_VARA="true"
ESTIMATE_PI="true"
ESTIMATE_GSCALE="true"
ESTGSCALE_ITER=500
REPORT_PLEIOTROPIC_QTL_EFFECT_MATRIX="false"
OUTPUT_MCMC_DELTA="false"



# ======= Julia Environment Setup (optional) =======
if [[ -f /cvmfs/hpc.ucdavis.edu/sw/conda/root/etc/profile.d/conda.sh ]]; then
  set +u
  source /cvmfs/hpc.ucdavis.edu/sw/conda/root/etc/profile.d/conda.sh
  conda activate nnmm_twas || true
  set -u
fi

JULIA_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
echo "Launching Julia with ${JULIA_THREADS} thread(s)"

# ======= Run the Julia Script =======
julia --project="$REPO_ROOT" --threads "$JULIA_THREADS" "$REPO_ROOT/scripts/run_nonmpi.jl" \
  --data_path "$DATA_PATH" \
  --analysis_path "$ANALYSIS_PATH" \
  --n_iter "$NITER" \
  --seed "$SEED" \
  --annot_file "$ANNOT_FILE" \
  --annot_dict "$ANNOT_DICT" \
  --out_freq "$OUTFREQ" \
  --starting_value_dir "$STARTING_VALUE_DIR" \
  --gscale_value_dir "$GSCALE_VALUE_DIR" \
  --st_path "$ST_PATH" \
  --thin "$THIN" \
  --n1 "$N1" \
  --n2 "$N2" \
  --n_con "$N_CON" \
  --is_continue "$IS_CONTINUE" \
  --estimate_vare "$ESTIMATE_VARE" \
  --estimate_vara "$ESTIMATE_VARA" \
  --estimate_pi "$ESTIMATE_PI" \
  --estimate_gscale "$ESTIMATE_GSCALE" \
  --estgscale_iter "$ESTGSCALE_ITER" \
  --report_pleiotropic_qtl_effect_matrix "$REPORT_PLEIOTROPIC_QTL_EFFECT_MATRIX" \
  --output_mcmc_delta "$OUTPUT_MCMC_DELTA"
