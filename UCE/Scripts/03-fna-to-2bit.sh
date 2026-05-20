#!/bin/bash
#SBATCH --job-name=fna_to_2bit
#SBATCH --nodes=1
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=128
#SBATCH --mem=128G                 
#SBATCH --output=/rna/liha/phylogenomics_practice/UCE/harvest_UCE/logs/faToTwoBit_%j.out
#SBATCH --error=/rna/liha/phylogenomics_practice/UCE/harvest_UCE/logs/faToTwoBit_%j.err
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate phyluce
echo "Job started on $(hostname) at $(date)" 
 

TSV_FILE="/rna/liha/phylogenomics_practice/UCE/harvest_UCE/00-species/species.tsv"
OUT_BASE="/rna/liha/phylogenomics_practice/UCE/harvest_UCE/03-two-bit"
IN_BASE="/bank/ncbi/gnathostome" 
MAX_JOBS=$SLURM_CPUS_PER_TASK  


do_convert() {
    # xargs로부터 전달받은 탭 구분 데이터를 변수로 분리
    IFS=$'\t' read -r Target_Clade Accession Organism_Name <<< "$1"
    
    FORMATTED_NAME=$(echo "$Organism_Name" | tr -d '\r' | tr ' ' '_')
    IN_DIR="${IN_BASE}/${Target_Clade}/ncbi_dataset/data/${Accession}"
    FNA_FILE=$(find "$IN_DIR" -maxdepth 1 -name "*_genomic.fna" | head -n 1)
    OUT_FILE="${OUT_BASE}/${FORMATTED_NAME}/${FORMATTED_NAME}.2bit"

    if [[ -f "$OUT_FILE" ]]; then
        return # 이미 있으면 종료
    fi

    mkdir -p "$(dirname "$OUT_FILE")"
    echo "Converting: $FORMATTED_NAME ..."
    faToTwoBit "$FNA_FILE" "$OUT_FILE"
}

export -f do_convert
export IN_BASE OUT_BASE 


tail -n +2 "$TSV_FILE" | awk -F'\t' '{print $1"\t"$2"\t"$8}' | \
xargs -I {} -P "$SLURM_CPUS_PER_TASK" -d $'\n' bash -c 'do_convert "$@"' _ {}

echo "All tasks completed!"