#!/bin/bash
#SBATCH --nodelist=guanine
#SBATCH --partition=intel
#SBATCH --nodes=1                 
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G                
#SBATCH --output=/rna/liha/phylogenomics_practice/UCE/harvest_UCE/logs/slice_%j.log    
#SBATCH --error=/rna/liha/phylogenomics_practice/UCE/harvest_UCE/logs/slice_%j.log      

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate phyluce
echo "Job started on $(hostname) at $(date)" 

# -----------------------------------------------------------------
BASE_DIR="/rna/liha/phylogenomics_practice/UCE/harvest_UCE"
TWOBIT_DIR="${BASE_DIR}/03-two-bit"
LASTZ_DIR="${BASE_DIR}/04-lastz-output"  # lastz.clean 파일들이 있는 곳 (추가됨)
OUTPUT_DIR="${BASE_DIR}/05-fasta-output" # 잘라낸 서열들이 저장될 새 폴더
CONF_FILE="${BASE_DIR}/taxon_set.conf"   # 종합 지도가 만들어질 위치
# -----------------------------------------------------------------

cd "$BASE_DIR"

# 1. 03-two-bit 폴더를 뒤져서 자동으로 taxon_set.conf 종합 지도를 만듭니다.
echo "========================================================"
echo " Creating taxon_set.conf automatically..."
echo "[scaffolds]" > "$CONF_FILE"

TAXA_LIST=$(find "$TWOBIT_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)
for TAXON in $TAXA_LIST; do
    echo "${TAXON}:${TWOBIT_DIR}/${TAXON}/${TAXON}.2bit" >> "$CONF_FILE"
done

cat "$CONF_FILE"
echo "========================================================"

# 2. 만들어진 지도를 바탕으로 서열을 잘라냅니다.
echo " Extracting UCE loci and filtering paralogs..."

# 기존 폴더 삭제 (입력 대기열 멈춤 에러 방지)
if [ -d "$OUTPUT_DIR" ]; then
    rm -rf "$OUTPUT_DIR"
fi


phyluce_probe_slice_sequence_from_genomes \
    --conf "$CONF_FILE" \
    --lastz "$LASTZ_DIR" \
    --output "$OUTPUT_DIR" \
    --name-pattern "metazoan_probes-70.fasta_v_{}.lastz.clean" \
    --flank 500

echo "All tasks completed at $(date)"
