#!/bin/bash
#SBATCH --exclude=bm[19-27]
#SBATCH --job-name=SBayesAPP_chr1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH --time=3-00:00:00
#SBATCH --partition=bmh
#SBATCH --account=qtlchenggrp
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jyqqu@ucdavis.edu

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$REPO_ROOT"
bash "$REPO_ROOT/script/run.sh"