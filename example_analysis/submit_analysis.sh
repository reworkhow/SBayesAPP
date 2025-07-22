#!/bin/bash

# ======= Input Arguments =======
DATA_PATH="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/example_analysis/SBayesAPP_input_first10blks/"
ANALYSIS_PATH="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/example_analysis/SBayesAPP_res_first10blks/"
mkdir -p "$ANALYSIS_PATH"
NITER=1000
SEED=42
NRANK=1
ANNOT_FILE="annotation_df.txt"
ANNOT_DICT="anno_matrix_dict"
OUTFREQ=100
STARTING_VALUE_DIR="XXX" # only needed if using starting values or extend chain length
SECONDARY_STARTING_VALUE_DIR="XXX" # only needed if using prefixed Gscale or extend chain length
ST_PATH="/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/example_analysis/ST_res/"
THIN=50
IS_CONTINUE="false"  # or "true" if continuing from a previous run

# ======= Julia Environment Setup (optional) =======
# module load julia/1.10.2  # if needed on HPC
source ~/.bashrc          # if Julia path is set there

# ======= Run the Julia Script =======
julia SBayesAPP_nonMPI.jl \
  "$DATA_PATH" \
  "$ANALYSIS_PATH" \
  "$NITER" \
  "$SEED" \
  "$NRANK" \
  "$ANNOT_FILE" \
  "$ANNOT_DICT" \
  "$OUTFREQ" \
  "$STARTING_VALUE_DIR" \
  "$SECONDARY_STARTING_VALUE_DIR" \
  "$ST_PATH" \
  "$THIN" \
  "$IS_CONTINUE"
