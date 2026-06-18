#!/bin/bash
#SBATCH --job-name=SLIT_Prnk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --output=/rna/liha/phylogenomics_practice/SLIT2/logs/prank_%j.log
#SBATCH --error=/rna/liha/phylogenomics_practice/SLIT2/logs/prank_%j.log

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate alignments

echo "Job started on $(hostname) at $(date)"

BASE_DIR="/rna/liha/phylogenomics_practice/SLIT2"
INPUT_FILE="${BASE_DIR}/Scripts/SLIT-homologs.fasta"
OUTPUT_DIR="${BASE_DIR}/Alignments"
OUTPUT_PREFIX="${OUTPUT_DIR}/SLIT_prank"
OUTPUT_FILE="${OUTPUT_DIR}/SLIT_aligned.fasta"

mkdir -p "$OUTPUT_DIR"

if [[ ! -s "$INPUT_FILE" ]]; then
    echo "Error: input file not found or empty: $INPUT_FILE"
    exit 1
fi

echo "Starting SLIT alignment with PRANK..."
echo "Input: $INPUT_FILE"
echo "Output prefix: $OUTPUT_PREFIX"

rm -f "${OUTPUT_PREFIX}".*

# No -t guide tree is supplied here. PRANK will estimate its own guide tree
# instead of using a fixed topology.
prank \
    -d="$INPUT_FILE" \
    -o="$OUTPUT_PREFIX" \
    -protein \
    -quiet

if [[ ! -s "${OUTPUT_PREFIX}.best.fas" ]]; then
    echo "Error: PRANK did not create expected output: ${OUTPUT_PREFIX}.best.fas"
    exit 1
fi

cp -f "${OUTPUT_PREFIX}.best.fas" "$OUTPUT_FILE"

echo "Alignment completed successfully! Output saved to: $OUTPUT_FILE"
echo "Job ended at $(date)"
