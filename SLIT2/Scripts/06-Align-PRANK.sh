#!/bin/bash
#SBATCH --job-name=PhyloPrank
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --output=logs/prank_%j.log
#SBATCH --error=logs/prank_%j.log
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate alignments

echo "Job started on $(hostname) at $(date)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PIPELINE_OUTPUT_PREFIX="${PIPELINE_OUTPUT_PREFIX:-SLIT}"
PIPELINE_MAIL_TO="${PIPELINE_MAIL_TO:-}"
PIPELINE_MAIL_SCRIPT="${PIPELINE_MAIL_SCRIPT:-}"

if [[ -z "${PIPELINE_INPUT_FASTA:-}" && -n "${SLIT_INPUT_FASTA:-}" ]]; then
    PIPELINE_INPUT_FASTA="$SLIT_INPUT_FASTA"
fi

INPUT_FILE="${PIPELINE_INPUT_FASTA:-${SCRIPT_DIR}/${PIPELINE_OUTPUT_PREFIX}-homologs.fasta}"
OUTPUT_DIR="${BASE_DIR}/Alignments"
OUTPUT_PREFIX="${OUTPUT_DIR}/${PIPELINE_OUTPUT_PREFIX}_prank"
OUTPUT_FILE="${OUTPUT_DIR}/${PIPELINE_OUTPUT_PREFIX}_aligned.fasta"
PRANK_INPUT_FILE="${OUTPUT_PREFIX}_input.fasta"

mkdir -p "$OUTPUT_DIR"

if [[ ! -s "$INPUT_FILE" ]]; then
    echo "Error: input file not found or empty: $INPUT_FILE"
    exit 1
fi

echo "Starting ${PIPELINE_OUTPUT_PREFIX} alignment with PRANK..."
echo "Input: $INPUT_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo "PRANK-safe input: $PRANK_INPUT_FILE"

rm -f "${OUTPUT_PREFIX}".*

awk '
    /^>/ {
        header = substr($0, 2)
        gsub(/[^A-Za-z0-9_.-]+/, "_", header)
        print ">" header
        next
    }
    { print }
' "$INPUT_FILE" > "$PRANK_INPUT_FILE"

send_mail() {
    local code="$?" subject log_file body mail_status
    trap - EXIT
    [[ "$code" -eq 0 ]] && subject="finished" || subject="failed"
    log_file="${BASE_DIR}/logs/prank_${SLURM_JOB_ID:-manual}.log"
    body="$(printf 'PRANK %s\nJob: %s %s\nInput: %s\nOutput: %s\n\nRecent log:\n%s\n' "$subject" "${SLURM_JOB_NAME:-manual}" "${SLURM_JOB_ID:-manual}" "$INPUT_FILE" "$OUTPUT_FILE" "$(tail -n 40 "$log_file" 2>/dev/null || echo "Log file not found: $log_file")")"

    if [[ -n "$PIPELINE_MAIL_TO" && -s "$PIPELINE_MAIL_SCRIPT" ]]; then
        set +e
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate mail
        python3 "$PIPELINE_MAIL_SCRIPT" \
            --to "$PIPELINE_MAIL_TO" \
            --subject "PRANK alignment ${subject}: ${SLURM_JOB_ID:-manual}" \
            --body "$body"
        mail_status=$?
        conda deactivate
        set -e

        if [[ "$mail_status" -ne 0 ]]; then
            echo "Warning: mail notification failed with status ${mail_status}."
        fi
    fi

    exit "$code"
}
trap send_mail EXIT

prank \
    -d="$PRANK_INPUT_FILE" \
    -o="$OUTPUT_PREFIX" \
    -iterate="30" \
    -protein \
    -quiet

if [[ ! -s "${OUTPUT_PREFIX}.best.fas" ]]; then
    echo "Error: PRANK did not create expected output: ${OUTPUT_PREFIX}.best.fas"
    exit 1
fi

cp -f "${OUTPUT_PREFIX}.best.fas" "$OUTPUT_FILE"

echo "Alignment completed successfully! Output saved to: $OUTPUT_FILE"
echo "Job ended at $(date)"
