#!/bin/bash
#SBATCH --nodelist=guanine
#SBATCH --partition=intel
#SBATCH --nodes=1                 
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=32         
#SBATCH --mem=128G                
#SBATCH --output=/rna/liha/phylogenomics_practice/UCE/harvest_UCE/logs/%x_%j.log    
#SBATCH --error=/rna/liha/phylogenomics_practice/UCE/harvest_UCE/logs/error_%x_%j.log       
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate NCBI-download
echo "Job started on $(hostname) at $(date)" 


if [ -z "$1" ]; then
    echo "Error: 타겟 Clade 이름이 입력되지 않았습니다."
    exit 1
fi


# 디렉토리 경로 변수 지정
WORK_DIR="/bank/ncbi/gnathostome"

echo "Node hostname: $(hostname)"
echo "Trying to access: $WORK_DIR"


if [ -d "$WORK_DIR" ]; then
    cd "$WORK_DIR"
    echo "Successfully entered $WORK_DIR"
else
    echo "CRITICAL ERROR: Directory not found on node $(hostname)!"
    echo "List of /bank/ncbi:"
    ls -l /bank/ncbi
    exit 1
fi
 

# 설정 변수
INPUT_TSV="/rna/liha/phylogenomics_practice/UCE/harvest_UCE/00-species/gnathostome_reference.tsv"      
TARGET_CLADE="$1"    
SAFE_TARGET_CLADE=$(echo "$TARGET_CLADE" | tr ' ' '_')

ACC_LIST="${SAFE_TARGET_CLADE}_accessions.txt" 
OUTPUT_ZIP="${SAFE_TARGET_CLADE}_data.zip"  

echo "========================================================"
echo " Extracting Accessions for: $TARGET_CLADE"
echo "========================================================"

# 1. TSV 파일에서 해당 Clade의 Accession(2번째 컬럼) 추출
# awk를 사용하여 첫번째 컬럼($1)이 정확히 TARGET_CLADE와 일치하는 경우 두번째 컬럼($2) 출력
awk -F '\t' -v target="$TARGET_CLADE" '$1 == target {print $2}' "$INPUT_TSV" > "$ACC_LIST"

# 추출된 개수 확인
COUNT=$(wc -l < "$ACC_LIST")
echo "Found $COUNT genomes for $TARGET_CLADE."

echo "========================================================"
echo " Downloading Dehydrated Zip (Metadata only)..."
echo "========================================================"

# 2. datasets 툴을 사용하여 일괄 다운로드
# --include genome: 게놈 서열 (fna)
datasets download genome accession \
    --inputfile "$ACC_LIST" \
    --include genome \
    --dehydrated \
    --filename "$OUTPUT_ZIP"

echo "--------------------------------------------------------"
echo " Unzipping structure..."



unzip -o "$OUTPUT_ZIP" -d "${SAFE_TARGET_CLADE}"
datasets rehydrate \
    --directory "${SAFE_TARGET_CLADE}" \
    --max-workers 16

echo "All tasks completed at $(date)" 


# 3. 무결성 검증. 
MD5_DIR="${SAFE_TARGET_CLADE}" 

if [ ! -d "$MD5_DIR" ]; then
    echo "Error: $MD5_DIR 디렉토리를 찾을 수 없습니다. 압축 해제에 실패했는지 확인하세요."
    exit 1
fi 

cd "$MD5_DIR" 

if [ ! -f "md5sum.txt" ]; then
    echo "No md5sum file" 
    exit 1 
fi 

if ! md5sum -c md5sum.txt --quiet; then 
    echo "download corrupted" 
    exit 1 
fi 
