#!/bin/bash
#SBATCH --job-name=run_iqtree
#SBATCH --nodes=1                 
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=128
#SBATCH --mem=128G                 
#SBATCH --output=/rna/liha/phylogenomics_practice/UCE/logs/iqtree_%j.log    
#SBATCH --error=/rna/liha/phylogenomics_practice/UCE/logs/iqtree_%j.log      

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate iqtree


BASE_DIR="/rna/liha/phylogenomics_practice/UCE"
OUTPUT_DIR="${BASE_DIR}/Tree"  
 
INPUT_FASTA="/rna/liha/phylogenomics_practice/UCE/Aligned_fasta/uce-4587_aligned_modified.fasta"
UCE_ID="uce-4587_aligned"


UCE_SPECIFIC_DIR="${OUTPUT_DIR}/${UCE_ID}"

if [ -d "$UCE_SPECIFIC_DIR" ]; then
    echo "Directory exists for ${UCE_ID}. Cleaning up old files..."
    rm -f "${UCE_SPECIFIC_DIR}"/*
else
    mkdir -p "$UCE_SPECIFIC_DIR"
fi

OUT_PREFIX="${UCE_SPECIFIC_DIR}/${UCE_ID}" 



iqtree3 -s "$INPUT_FASTA" -B 1000 --prefix "$OUT_PREFIX" -T AUTO 



echo "All IQ-TREE tasks completed successfully at $(date)!" 