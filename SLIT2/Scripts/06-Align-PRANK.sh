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
PRANK_INPUT_FILE="${OUTPUT_PREFIX}_input.fasta"

mkdir -p "$OUTPUT_DIR"

if [[ ! -s "$INPUT_FILE" ]]; then
    echo "Error: input file not found or empty: $INPUT_FILE"
    exit 1
fi

echo "Starting SLIT alignment with PRANK..."
echo "Input: $INPUT_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo "PRANK-safe input: $PRANK_INPUT_FILE"

rm -f "${OUTPUT_PREFIX}".*

# Convert FASTA headers to single-token labels before PRANK reads them.
awk '
    /^>/ {
        header = substr($0, 2)
        gsub(/[^A-Za-z0-9_.-]+/, "_", header)
        print ">" header
        next
    }
    { print }
' "$INPUT_FILE" > "$PRANK_INPUT_FILE"

# Send mail if PRANK fails or suceeds in whatever reason. 
send_mail() {
    local code="$?" subject log_file body
    log_file="${BASE_DIR}/logs/prank_${SLURM_JOB_ID}.log"
    [[ "$code" -eq 0 ]] && subject="finished" || subject="failed"
    body="$(printf 'PRANK %s\nJob: %s %s\nInput: %s\nOutput: %s\n\nRecent log:\n%s\n' "$subject" "$SLURM_JOB_NAME" "$SLURM_JOB_ID" "$INPUT_FILE" "$OUTPUT_FILE" "$(tail -n 40 "$log_file" 2>/dev/null || echo "Log file not found: $log_file")")"
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate mail
    python3 /rna/liha/tools/send_mail/mail.py --to "siwawa@snu.ac.kr" --subject "PRANK alignment ${subject}: ${SLURM_JOB_ID}" --body "$body"
    conda deactivate
    exit "$code"
}
trap send_mail EXIT


# No -t guide tree is supplied here. PRANK will estimate its own guide tree
# instead of using a fixed topology.
prank \
    -d="$PRANK_INPUT_FILE" \
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



if [[ "$N_DONE" -eq "$N_TOTAL" ]]; then
    # Send mail if finished 
    MAIL_ENV="/rna/liha/tools/send_mail"
    MAIL_SCRIPT="${MAIL_ENV}/mail.py" 
    BODY=$(head -n 20 ${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log 2>/dev/null || echo "Log file not found") 

    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate mail 
    python3 "$MAIL_SCRIPT" \
        --to "siwawa@snu.ac.kr" \
        --subject "BLASTp search Finished" \
        --body "$BODY"

    echo "Done at $(date)"
    conda deactivate 
fi