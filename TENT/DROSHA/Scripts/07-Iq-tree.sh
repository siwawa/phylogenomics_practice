#!/bin/bash
#SBATCH --job-name=PhyloTree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --output=logs/iqtree_%j.log
#SBATCH --error=logs/iqtree_%j.log
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate iqtree

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
INPUT_FASTA="${PIPELINE_ALIGNMENT_FASTA:-${BASE_DIR}/Alignments/${PIPELINE_OUTPUT_PREFIX}_aligned.fasta}"
OUTPUT_DIR="${BASE_DIR}/Tree"
OUT_PREFIX="${OUTPUT_DIR}/${PIPELINE_OUTPUT_PREFIX}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "$OUTPUT_DIR"

if [[ ! -s "$INPUT_FASTA" ]]; then
    echo "Error: alignment FASTA not found or empty: $INPUT_FASTA"
    exit 1
fi

if compgen -G "${OUT_PREFIX}.*" > /dev/null; then
    echo "Old tree files found for ${OUT_PREFIX}. Cleaning up..."
    rm -f "${OUT_PREFIX}".*
fi

echo "Starting IQ-TREE analysis for ${PIPELINE_OUTPUT_PREFIX}..."
echo "Input: $INPUT_FASTA"
echo "Output prefix: $OUT_PREFIX"

iqtree3 -s "$INPUT_FASTA" \
        -B 5000 \
        --prefix "$OUT_PREFIX" \
        -T AUTO \
        --threads-max "$THREADS"

echo "IQ-TREE analysis completed at $(date)."
