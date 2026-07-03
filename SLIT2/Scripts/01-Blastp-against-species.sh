#!/bin/bash
#SBATCH --job-name=PhyloBlastp
#SBATCH --output=Blast/logs/%A_%a.out
#SBATCH --error=Blast/logs/%A_%a.err
#SBATCH --array=2-47%3
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PIPELINE_QUERY_NAME="${PIPELINE_QUERY_NAME:-SLIT1-2-3}"
PIPELINE_MAIL_TO="${PIPELINE_MAIL_TO:-}"
PIPELINE_MAIL_SCRIPT="${PIPELINE_MAIL_SCRIPT:-}"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate NCBI-download

SPECIES_FILE="${SCRIPT_DIR}/species.txt"
QUERY_DIR="${BASE_DIR}/Blast/Query"
OUTPUT_DIR="${BASE_DIR}/Blast/Results"
DONE_DIR="${BASE_DIR}/Blast/Done"

mkdir -p "$OUTPUT_DIR" "$DONE_DIR"

ROW="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_FILE")"
if [[ -z "$ROW" ]]; then
    echo "No species row found for array task ${SLURM_ARRAY_TASK_ID}."
    exit 1
fi

TARGET_CLADE="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $1}')"
TAXID="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $7}')"
SPECIES="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $8}')"
SPECIES_LABEL="$(printf '%s\n' "$SPECIES" | awk '{ if (NF >= 2) print $1, $2; else print $0 }')"

SPECIES_SAFE="${TARGET_CLADE}_$(echo "$SPECIES_LABEL" | tr ' ' '_')"
SPECIES_DIR="${OUTPUT_DIR}/${SPECIES_SAFE}"
DONE_FILE="${DONE_DIR}/${SPECIES_SAFE}.done"
QUERY="${QUERY_DIR}/${PIPELINE_QUERY_NAME}.fasta"
RESULT_FILE="${SPECIES_DIR}/${TARGET_CLADE}_${PIPELINE_QUERY_NAME}.txt"

mkdir -p "$SPECIES_DIR"

if [[ ! -s "$QUERY" ]]; then
    echo "Error: query FASTA not found or empty: $QUERY"
    exit 1
fi

if [[ -f "$DONE_FILE" ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} | ${SPECIES_SAFE} already completed. Skipping BLAST."
else
    echo "Array task ${SLURM_ARRAY_TASK_ID} | taxid: ${TAXID} | species: ${SPECIES_SAFE}"
    echo "Running BLASTP: ${PIPELINE_QUERY_NAME} x ${SPECIES_SAFE}"

    # Output format: QueryID, TargetID, PID, Alignment length, query coverage, E-value, Bitscore, Subject title
    blastp \
        -query "$QUERY" \
        -db refseq_protein \
        -remote \
        -entrez_query "txid${TAXID}[Organism]" \
        -out "$RESULT_FILE" \
        -outfmt "6 qseqid saccver pident length qcovs evalue bitscore stitle" \
        -max_target_seqs 10

    echo "Done: ${PIPELINE_QUERY_NAME} x ${SPECIES_SAFE}"
    touch "$DONE_FILE"
fi

send_finished_mail() {
    [[ -n "$PIPELINE_MAIL_TO" && -s "$PIPELINE_MAIL_SCRIPT" ]] || return 0

    local body mail_status
    body="BLASTp search finished for ${PIPELINE_QUERY_NAME} in ${BASE_DIR}."

    set +e
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate mail
    python3 "$PIPELINE_MAIL_SCRIPT" \
        --to "$PIPELINE_MAIL_TO" \
        --subject "BLASTp search finished: ${PIPELINE_QUERY_NAME}" \
        --body "$body"
    mail_status=$?
    conda deactivate
    set -e

    if [[ "$mail_status" -ne 0 ]]; then
        echo "Warning: mail notification failed with status ${mail_status}."
    fi
}

N_DONE="$(find "$DONE_DIR" -name "*.done" | wc -l | tr -d ' ')"
N_TOTAL="$(awk 'NR > 1 && NF { n++ } END { print n + 0 }' "$SPECIES_FILE")"

if [[ "$N_TOTAL" -gt 0 && "$N_DONE" -ge "$N_TOTAL" ]]; then
    send_finished_mail
    echo "Done at $(date)"
fi
