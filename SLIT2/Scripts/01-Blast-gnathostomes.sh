#!/bin/bash
#SBATCH --job-name=slit_blast
#SBATCH --output=/rna/liha/phylogenomics_practice/SLIT2/Blast/Logs/%A_%a.out
#SBATCH --error=/rna/liha/phylogenomics_practice/SLIT2/Blast/Logs/%A_%a.err
#SBATCH --array=1-25
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate NCBI-download 

SPECIES_FILE="/rna/liha/phylogenomics_practice/SLIT2/Scripts/species.txt"
QUERY_DIR="/rna/liha/phylogenomics_practice/SLIT2/Blast/Query"
OUTPUT_DIR="/rna/liha/phylogenomics_practice/SLIT2/Blast/Results"
DONE_DIR="/rna/liha/phylogenomics_practice/SLIT2/Blast/Done"

mkdir -p "$OUTPUT_DIR" "$DONE_DIR"

TAXID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_FILE" | awk '{print $8}')
SPECIES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_FILE" | awk '{print $10, $11}')
SPECIES_SAFE=$(echo "$SPECIES" | tr ' ' '_')


echo "Array task ${SLURM_ARRAY_TASK_ID} | taxid: ${TAXID}"

for QUERY_PTN in SLIT1 SLIT2 SLIT3; do
    QUERY="${QUERY_DIR}/${QUERY_PTN}.fasta"

    echo "Running BLASTP: ${QUERY_PTN} x txid${TAXID}"

    blastp \
        -query "$QUERY" \
        -db nr \
        -remote \
        -entrez_query "txid${TAXID}[Organism]" \
        -out "${OUTPUT_DIR}/blast_${QUERY_PTN}_${TAXID}.txt" \
        -outfmt 6 \
        -max_target_seqs 10

    echo "Done: ${QUERY_PTN} x txid${TAXID}"
done

touch "${DONE_DIR}/${SPECIES_SAFE}.done"




if [[ "$SLURM_ARRAY_TASK_ID" == "$SLURM_ARRAY_TASK_MAX" ]]; then
    # Send mail if finished 
    MAIL_ENV="/rna/liha/tools/send_mail"
    MAIL_SCRIPT="${MAIL_ENV}/mail.py" 
    BODY=$(head -n 20 ${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log 2>/dev/null || echo "Log file not found") 

    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate mail 
    python3 "$MAIL_SCRIPT" \
        --to "siwawa@snu.ac.kr" \
        --subject "[ALL-BY_ALL COMPARISON] Finished" \
        --body "$BODY"

    echo "Done at $(date)"
    conda deactivate 
fi






