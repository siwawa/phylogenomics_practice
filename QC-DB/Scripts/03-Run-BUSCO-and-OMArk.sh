#!/bin/bash
#SBATCH --job-name=proteome_qc
#SBATCH --output=/rna/liha/phylogenomics_practice/QC-DB/Scripts/logs/qc_%A_%a.out
#SBATCH --error=/rna/liha/phylogenomics_practice/QC-DB/Scripts/logs/qc_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --nodelist=pichu
#SBATCH --array=1-70%5
#
# Runs 03-BUSCO-and-OMArk.sh once per taxid, reading taxid + selected
# accession from selected_assemblies.tsv (one line per taxid, output of
# 02-Select-assemblies.sh). "%5" throttles Slurm to at most 5 concurrently
# running array tasks.
#
# NOTE: the #SBATCH --array=1-70%5 above is just a fallback default.
# Always override it at submission time to match the actual row count:
#
#   N=$(($(wc -l < selected_assemblies.tsv) - 1))   # minus 1 for the header
#   sbatch --array=1-${N}%5 04-Submit-array.sbatch

set -euo pipefail

SCRIPT_DIR=/rna/liha/phylogenomics_practice/QC-DB/Scripts
ASSEMBLIES_FILE=${SCRIPT_DIR}/selected_assemblies.tsv

mkdir -p "${SCRIPT_DIR}/logs"

# +1 to skip the header line, so array index 1 => data row 1
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${ASSEMBLIES_FILE}")

TAXID=$(echo "$LINE" | cut -f1)
TAXONOMY=$(echo "$LINE" | cut -f2)
ACCESSION=$(echo "$LINE" | cut -f3)

if [[ -z "$TAXID" || -z "$ACCESSION" ]]; then
    echo "[ERROR] No taxid/accession found for array index ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

bash "${SCRIPT_DIR}/03-BUSCO-and-OMArk.sh" "$TAXID" "$ACCESSION" "$TAXONOMY"
