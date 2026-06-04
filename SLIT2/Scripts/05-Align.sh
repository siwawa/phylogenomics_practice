#!/bin/bash
#SBATCH --job-name=SLIT_Mafft
#SBATCH --nodes=1                 
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=32      
#SBATCH --mem=32G                 
#SBATCH --output=/rna/liha/phylogenomics_practice/SLIT2/logs/mafft_%j.log    
#SBATCH --error=/rna/liha/phylogenomics_practice/SLIT2/logs/mafft_%j.log      

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate alignments 

echo "Job started on $(hostname) at $(date)" 

# 1. 경로 설정 (방금 파이썬으로 만든 파일 기준)
BASE_DIR="/rna/liha/phylogenomics_practice/SLIT2"
INPUT_FILE="${BASE_DIR}/Scripts/Blast-filtered-longest-isoforms.fasta"

# 정렬된 결과를 저장할 폴더 및 파일명 지정
OUTPUT_DIR="${BASE_DIR}/Alignments"
OUTPUT_FILE="${OUTPUT_DIR}/SLIT_aligned.fasta"

# 2. 출력 폴더 생성 (없으면 만들기)
mkdir -p "$OUTPUT_DIR"

echo "Starting SLIT alignment with MAFFT --localpair..." 

# 3. MAFFT 실행 로직 (단일 파일이므로 반복문 불필요)
if [ -f "$INPUT_FILE" ]; then
    echo "Processing SLIT sequences..."
    
    # L-INS-i 알고리즘(--localpair) 적용 및 SLURM 코어 수 연동
    mafft --localpair --maxiterate 1000 --thread $SLURM_CPUS_PER_TASK "$INPUT_FILE" > "$OUTPUT_FILE"
    
    echo "Alignment completed successfully! Output saved to: $OUTPUT_FILE"
else
    echo "Error: Input file ($INPUT_FILE) not found!"
    exit 1
fi

echo "Job ended at $(date)"