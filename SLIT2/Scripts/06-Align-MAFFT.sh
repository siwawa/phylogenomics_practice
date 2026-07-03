#!/bin/bash
#SBATCH --job-name=PhyloMafft
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --output=logs/mafft_%j.log
#SBATCH --error=logs/mafft_%j.log
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate alignments

echo "Job started on $(hostname) at $(date)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PIPELINE_OUTPUT_PREFIX="${PIPELINE_OUTPUT_PREFIX:-SLIT}"

if [[ -z "${PIPELINE_INPUT_FASTA:-}" && -n "${SLIT_INPUT_FASTA:-}" ]]; then
    PIPELINE_INPUT_FASTA="$SLIT_INPUT_FASTA"
fi

INPUT_FILE="${PIPELINE_INPUT_FASTA:-${SCRIPT_DIR}/${PIPELINE_OUTPUT_PREFIX}-homologs.fasta}"
OUTPUT_DIR="${BASE_DIR}/Alignments"
OUTPUT_FILE="${OUTPUT_DIR}/${PIPELINE_OUTPUT_PREFIX}_aligned.fasta"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "$OUTPUT_DIR"

if [[ ! -s "$INPUT_FILE" ]]; then
    echo "Error: input file not found or empty: $INPUT_FILE"
    exit 1
fi

echo "Starting ${PIPELINE_OUTPUT_PREFIX} alignment with MAFFT --localpair..."
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"

mafft --localpair --maxiterate 1000 --thread "$THREADS" "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Alignment complete. Output: $OUTPUT_FILE"
echo "Job ended at $(date)"
