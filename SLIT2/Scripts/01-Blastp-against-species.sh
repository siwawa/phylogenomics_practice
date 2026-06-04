#!/bin/bash
#SBATCH --job-name=SLIT-Blastp
#SBATCH --output=/rna/liha/phylogenomics_practice/SLIT2/Blast/logs/%A_%a.out
#SBATCH --error=/rna/liha/phylogenomics_practice/SLIT2/Blast/logs/%A_%a.err
#SBATCH --array=1-46%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
set -euo pipefail


source ~/miniconda3/etc/profile.d/conda.sh
conda activate NCBI-download 

SPECIES_FILE="/rna/liha/phylogenomics_practice/SLIT2/Scripts/species.txt"
QUERY_DIR="/rna/liha/phylogenomics_practice/SLIT2/Blast/Query"
OUTPUT_DIR="/rna/liha/phylogenomics_practice/SLIT2/Blast/Results"
DONE_DIR="/rna/liha/phylogenomics_practice/SLIT2/Blast/Done"

mkdir -p "$OUTPUT_DIR" "$DONE_DIR"

TAXID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_FILE" | awk '{print $8}')
SPECIES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_FILE" | awk '{print $9, $10}')
SPECIES_SAFE=$(echo "$SPECIES" | tr ' ' '_')
SPECIES_DIR="${OUTPUT_DIR}/${SPECIES_SAFE}"
DONE_FILE="${DONE_DIR}/${SPECIES_SAFE}.done"

mkdir -p "$SPECIES_DIR"


if [[ -f "$DONE_FILE" ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} | ${SPECIES_SAFE} already completed. Skipping BLAST."
else
    echo "Array task ${SLURM_ARRAY_TASK_ID} | taxid: ${SPECIES_SAFE}"

    for QUERY_PTN in SLIT1 SLIT2 SLIT3; do
        QUERY="${QUERY_DIR}/${QUERY_PTN}.fasta"

        echo "Running BLASTP: ${QUERY_PTN} x ${SPECIES_SAFE}"

        # Output format: 
        # QueryID, TargetID, PID, Alignment length(include internal gaps), E-value, Bitscore, Subject title
        blastp \
            -query "$QUERY" \
            -db nr_cluster_seq \
            -remote \
            -entrez_query "txid${TAXID}[Organism]" \
            -out "${SPECIES_DIR}/${QUERY_PTN}.txt" \
            -outfmt "6 qseqid saccver pident length qcovs evalue bitscore stitle" \
            -max_target_seqs 10

        echo "Done: ${QUERY_PTN} x ${SPECIES_SAFE}"
        sleep 45
    done

    touch "$DONE_FILE"
fi


N_DONE=$(find "$DONE_DIR" -name "*.done" | wc -l)
N_TOTAL=$(wc -l < "$SPECIES_FILE")

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






