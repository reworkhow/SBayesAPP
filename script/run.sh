#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ======= Input Arguments =======
DATA_PATH="$REPO_ROOT/example/SBayesAPP_input_first10blks/"
ANALYSIS_PATH="$REPO_ROOT/example/SBayesAPP_res_first10blks/"
mkdir -p "$ANALYSIS_PATH"
NITER=1000
SEED=42
NRANK=1
ANNOT_FILE="annotation_df.txt"
ANNOT_DICT="anno_matrix_dict"
OUTFREQ=100
STARTING_VALUE_DIR="XXX" # only needed if using starting values or extend chain length
GSCALE_VALUE_DIR="XXX" # only needed if using fixed Gscale values during continuation
ST_PATH="$REPO_ROOT/example/ST_res/"
THIN=50
IS_CONTINUE="false"  # or "true" if continuing from a previous run

# ======= Julia Environment Setup (optional) =======
if [[ -f /cvmfs/hpc.ucdavis.edu/sw/conda/root/etc/profile.d/conda.sh ]]; then
  source /cvmfs/hpc.ucdavis.edu/sw/conda/root/etc/profile.d/conda.sh
  conda activate nnmm_twas || true
fi

# ======= Run the Julia Script =======
julia --project="$REPO_ROOT" "$REPO_ROOT/scripts/run_example.jl"
