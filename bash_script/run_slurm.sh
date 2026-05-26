#!/bin/bash
#SBATCH --exclude=bm[19-27]
#SBATCH --job-name=SBayesAPP_first10blks
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=3-00:00:00
#SBATCH --partition=bmh
#SBATCH --account=qtlchenggrp
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jyqqu@ucdavis.edu

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"

if [[ -f "$SUBMIT_DIR/run.sh" ]]; then
	SCRIPT_DIR="$SUBMIT_DIR"
elif [[ -f "$SUBMIT_DIR/script/run.sh" ]]; then
	SCRIPT_DIR="$SUBMIT_DIR/script"
else
	echo "Could not locate run.sh from submit directory: $SUBMIT_DIR" >&2
	exit 1
fi

REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$REPO_ROOT"
bash "$SCRIPT_DIR/run.sh"