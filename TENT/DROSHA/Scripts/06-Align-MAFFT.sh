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

if [[ -n "${PIPELINE_BASE_DIR:-}" && -d "${PIPELINE_BASE_DIR}/Scripts" ]]; then
    BASE_DIR="$(cd "$PIPELINE_BASE_DIR" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/Scripts" ]]; then
    BASE_DIR="$(cd "$SLURM_SUBMIT_DIR" && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi
SCRIPT_DIR="${BASE_DIR}/Scripts"
PIPELINE_OUTPUT_PREFIX="${PIPELINE_OUTPUT_PREFIX:-SLIT}"

INPUT_FILE="${SCRIPT_DIR}/${PIPELINE_OUTPUT_PREFIX}-homologs.fasta"
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
