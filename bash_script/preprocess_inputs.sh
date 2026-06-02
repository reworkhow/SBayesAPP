#!/bin/bash
#SBATCH --exclude=bm[19-27]
#SBATCH --job-name=SBayesAPP_preprocess
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --partition=bmh
#SBATCH --account=qtlchenggrp
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jyqqu@ucdavis.edu

set -euo pipefail

usage() {
	cat <<'EOF'
Usage:
  1. Edit the input variables in this script.
  2. Submit with:
     sbatch bash_script/preprocess_inputs.sh

Example:
  out_path=/path/to/output
  ld_info_path=/path/to/LD_data/eigen_data_whole_genome
  trait1_file=/path/to/trait1.imputed.ma
  trait2_file=/path/to/trait2.imputed.ma
  ldinfo_file=snp.info
  annot_file=/path/to/annotation_df.txt
  sbatch bash_script/preprocess_inputs.sh

Notes:
  - Adjust the #SBATCH resource lines in this file if your job needs different resources.
  - LDinfo_file can be either an absolute path or a filename relative to LD_info_path.
EOF
}

out_path="/path/to/output"
ld_info_path="/path/to/LD_data/eigen_data_whole_genome"
trait1_file="/path/to/trait1.imputed.ma"
trait2_file="/path/to/trait2.imputed.ma"
ldinfo_file="snp.info"
annot_file="/path/to/annotation_df.txt"
annot_dict_name="anno_matrix_dict"
readable_files_dir="readableFiles"
nblock=""

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
	usage
	exit 0
fi

for required_var in out_path ld_info_path trait1_file trait2_file ldinfo_file annot_file; do
	if [[ -z "${!required_var}" || "${!required_var}" == /path/to/* ]]; then
		echo "Please set ${required_var} near the top of bash_script/preprocess_inputs.sh before submitting." >&2
		usage >&2
		exit 1
	fi
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$REPO_ROOT"

if ! command -v julia >/dev/null 2>&1; then
	echo "julia is not available on PATH" >&2
	exit 1
fi

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

cmd=(
	julia --project=. scripts/preprocess_inputs.jl
	--out "$out_path"
	--LD_info_path "$ld_info_path"
	--trait1_file "$trait1_file"
	--trait2_file "$trait2_file"
	--LDinfo_file "$ldinfo_file"
	--annot_file "$annot_file"
	--annot_dict_name "$annot_dict_name"
	--readable_files_dir "$readable_files_dir"
)

if [[ -n "$nblock" ]]; then
	cmd+=(--nblock "$nblock")
fi

printf 'Running command:'
	printf ' %q' "${cmd[@]}"
	printf '\n'

srun "${cmd[@]}"