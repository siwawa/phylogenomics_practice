#!/bin/bash
#SBATCH --job-name=uce_align
#SBATCH --nodes=1                 
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=128
#SBATCH --mem=128G                 
#SBATCH --output=/rna/liha/phylogenomics_practice/UCE/logs/mafft_%j.log    
#SBATCH --error=/rna/liha/phylogenomics_practice/UCE/logs/mafft_%j.log      
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate alignments 
echo "Job started on $(hostname) at $(date)" 


BASE_DIR="/rna/liha/phylogenomics_practice/UCE"
INPUT_DIR="${BASE_DIR}/Unaligned_fasta" 
OUTPUT_DIR="${BASE_DIR}/Aligned_fasta"  
STATS_FILE="${BASE_DIR}/uce_stats_top10.tsv" 
THREADS=$SLURM_CPUS_PER_TASK
echo "Starting UCE alignment with MAFFT --localpair..." 


mkdir -p $OUTPUT_DIR
SELECTED_LOCI=$(awk -F'\t' 'NR>1 {print $1}' "$STATS_FILE") 

for locus in $SELECTED_LOCI; do
    file="${INPUT_DIR}/${locus}.fasta"
    aligned_file="${OUTPUT_DIR}/${locus}_aligned.fasta"
    
    if [ -f "$file" ]; then
        if [ -f "$aligned_file" ]; then
            echo "Skipping $locus... (Already aligned)"
        else 
            echo "Processing $locus..."
            mafft --localpair --maxiterate 1000 --thread $THREADS "$file" > "$aligned_file"
        fi
    else
        echo "Warning: File for $locus not found in $INPUT_DIR, skipping."
    fi
done

echo "Alignment completed for selected UCEs."



