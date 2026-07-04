#!/bin/bash
#SBATCH --job-name=PhyloRun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=logs/__%j.log
#SBATCH --error=logs/__%j.log
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PIPELINE_OUTPUT_PREFIX="${PIPELINE_OUTPUT_PREFIX:-SLIT}"

cd "$BASE_DIR"
mkdir -p logs

export PIPELINE_INPUT_FASTA="${PIPELINE_INPUT_FASTA:-${SCRIPT_DIR}/${PIPELINE_OUTPUT_PREFIX}-homologs-removed-divergent.fasta}"

bash "${SCRIPT_DIR}/00-Wrapper.sh" \
  --from align \
  --to compare_tree \
  --aligner prank \
  --force-local
