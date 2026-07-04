#!/bin/bash
#SBATCH --job-name=TENT5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=logs/__%j.log
#SBATCH --error=logs/__%j.log
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" && -s "${SLURM_SUBMIT_DIR}/Scripts/00-Wrapper.sh" ]]; then
    BASE_DIR="$(cd "$SLURM_SUBMIT_DIR" && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi
SCRIPT_DIR="${BASE_DIR}/Scripts"
PIPELINE_OUTPUT_PREFIX="${PIPELINE_OUTPUT_PREFIX:-SLIT}"

cd "$BASE_DIR"

export PIPELINE_OUTPUT_PREFIX="TENT5"
export PIPELINE_QUERY_NAME="TENT5A-B-C-D"
export PIPELINE_OUTGROUP_PATTERN="Echinodermata"

INPUT_FASTA="${SCRIPT_DIR}/${PIPELINE_OUTPUT_PREFIX}-homologs.fasta"
if [[ ! -s "$INPUT_FASTA" ]]; then
    echo "Error: alignment input FASTA not found or empty: $INPUT_FASTA" >&2
    exit 1
fi

bash "${SCRIPT_DIR}/00-Wrapper.sh"   --from align   --to compare_tree   --aligner prank   --force-local
