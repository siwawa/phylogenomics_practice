#!/usr/bin/env bash
#SBATCH --job-name=TENT2_BLASTp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --exclusive
#SBATCH --nodelist=guanine
#SBATCH --partition=intel
#SBATCH --output=logs/TENT2_BLASTp_%j.log
#SBATCH --error=logs/TENT2_BLASTp_%j.log


set -euo pipefail

BASE_DIR="/rna/liha/phylogenomics_practice/TENT/TENT2"
QUERY="${BASE_DIR}/TENT2.fasta"
DB_DIR="/bank/ncbi/gnathostome/refseq_protein"
DATABASE="${DB_DIR}/refseq_protein"
OUTPUT_DIR="${BASE_DIR}/BLASTp-results"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate NCBI-download

export BLASTDB="$DB_DIR"
mkdir -p "$OUTPUT_DIR"



OUTFMT="6 qseqid saccver staxids sscinames pident length qcovs evalue bitscore stitle"
MAX_TARGET_SEQS=50000

blastp \
    -query "$QUERY" \
    -db "$DATABASE" \
    -evalue 0.1 \
    -num_threads "${SLURM_CPUS_PER_TASK:-1}" \
    -max_target_seqs "$MAX_TARGET_SEQS" \
    -outfmt "$OUTFMT" \
    -out "${OUTPUT_DIR}/TENT2.txt"





# Send mail if finished. 
PIPELINE_QUERY_NAME="TENT2"
PIPELINE_MAIL_TO="${PIPELINE_MAIL_TO:-"siwawa@snu.ac.kr"}"
PIPELINE_MAIL_SCRIPT="${PIPELINE_MAIL_SCRIPT:-/rna/liha/tools/send_mail/mail.py}"

# Temporary function to send mail notification when the pipeline starts and finished. 
send_finished_mail() {
    [[ -n "$PIPELINE_MAIL_TO" && -s "$PIPELINE_MAIL_SCRIPT" ]] || return 0

    local body mail_status
    body="BLASTp search finished for ${PIPELINE_QUERY_NAME} in ${BASE_DIR}."

    set +e
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate mail
    python3 "$PIPELINE_MAIL_SCRIPT" \
        --to "$PIPELINE_MAIL_TO" \
        --subject "BLASTp search: ${PIPELINE_QUERY_NAME}" \
        --body "$body"
    mail_status=$?
    conda deactivate
    set -e

    if [[ "$mail_status" -ne 0 ]]; then
        echo "Warning: mail notification failed with status ${mail_status}."
    fi
} 


send_finished_mail 
echo "Done at $(date)"
