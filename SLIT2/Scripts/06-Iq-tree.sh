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

# 1. 경로 및 파일 설정 (이전 MAFFT 결과물 기준)
BASE_DIR="/rna/liha/phylogenomics_practice/SLIT2"
INPUT_FASTA="${BASE_DIR}/Alignments/SLIT_aligned.fasta"
OUTPUT_DIR="${BASE_DIR}/Tree"  

# 2. 출력 폴더 생성 및 Prefix 지정
mkdir -p "$OUTPUT_DIR"
OUT_PREFIX="${OUTPUT_DIR}/SLIT" 

# 3. 기존 파일 정리 (클린업 로직)
# IQ-TREE는 동일한 Prefix의 파일이 존재하면 에러를 뿜거나 실행을 멈춥니다.
if ls ${OUT_PREFIX}.* 1> /dev/null 2>&1; then
    echo "Old tree files found for ${OUT_PREFIX}. Cleaning up..."
    rm -f ${OUT_PREFIX}.*
fi

echo "Starting IQ-TREE analysis for SLIT..." 

# 4. IQ-TREE 실행
# -T AUTO를 사용하되, Slurm에서 할당한 코어 수를 넘지 못하도록 안전장치 추가
iqtree3 -s "$INPUT_FASTA" \
        -B 5000 \
        --prefix "$OUT_PREFIX" \
        -T AUTO \
        --threads-max $SLURM_CPUS_PER_TASK

echo "All IQ-TREE tasks completed successfully at $(date)!"