#!/bin/bash
#SBATCH --job-name=SLIT_Tree
#SBATCH --nodes=1                 
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=64       
#SBATCH --mem=64G                 
#SBATCH --output=/rna/liha/phylogenomics_practice/SLIT2/logs/iqtree_%j.log    
#SBATCH --error=/rna/liha/phylogenomics_practice/SLIT2/logs/iqtree_%j.log      

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate iqtree

echo "Job started on $(hostname) at $(date)" 

# 1. Set paths for the previous MAFFT alignment output.
BASE_DIR="/rna/liha/phylogenomics_practice/SLIT2"
INPUT_FASTA="${BASE_DIR}/Alignments/SLIT_alignment.fas"
OUTPUT_DIR="${BASE_DIR}/Tree"  

# 2. Create the output directory and define the IQ-TREE prefix.
mkdir -p "$OUTPUT_DIR"
OUT_PREFIX="${OUTPUT_DIR}/SLIT" 

# 3. Remove old files with the same prefix before starting.
# IQ-TREE can fail or stop when files with the same prefix already exist.
if ls ${OUT_PREFIX}.* 1> /dev/null 2>&1; then
    echo "Old tree files found for ${OUT_PREFIX}. Cleaning up..."
    rm -f ${OUT_PREFIX}.*
fi

echo "Starting IQ-TREE analysis for SLIT..." 

# 4. Run IQ-TREE.
# Use -T AUTO while capping threads at the CPU count allocated by SLURM.
iqtree3 -s "$INPUT_FASTA" \
        -B 5000 \
        --prefix "$OUT_PREFIX" \
        -T AUTO \
        --threads-max $SLURM_CPUS_PER_TASK

echo "IQ-TREE analysis completed at $(date)."