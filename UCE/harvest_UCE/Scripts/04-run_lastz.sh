#!/bin/bash
#SBATCH --nodes=1                 
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G                 
#SBATCH --output=/rna/liha/phylogenomics_practice/UCE/harvest_UCE/logs/lastz_%j.log    
#SBATCH --error=/rna/liha/phylogenomics_practice/UCE/harvest_UCE/logs/lastz_%j.log      

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate phyluce
echo "Job started on $(hostname) at $(date)" 


BASE_DIR="/rna/liha/phylogenomics_practice/UCE/harvest_UCE"
TWOBIT_DIR="${BASE_DIR}/03-two-bit"
PROBE_FILE="${BASE_DIR}/02-probes/metazoan_probes-70.fasta" 
OUTPUT_DIR="${BASE_DIR}/04-lastz-output"
DB_FILE="${BASE_DIR}/uce_search.sqlite"
CORES=$SLURM_CPUS_PER_TASK



cd "$BASE_DIR"

if [ -d "$OUTPUT_DIR" ]; then
    echo "[Info] Removing existing output directory: $OUTPUT_DIR"
    rm -rf "$OUTPUT_DIR"
fi

if [ -f "$DB_FILE" ]; then
    echo "[Info] Removing existing database file: $DB_FILE"
    rm -f "$DB_FILE"
fi 



# 03-two-bit 폴더에 있는 모든 .2bit 파일의 이름을 추출해서 리스트로 만듭니다.
TAXA_LIST=$(find "$TWOBIT_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | tr '\n' ' ')

if [[ -z "$TAXA_LIST" ]]; then
    echo "[Error] No .2bit files found directly in $TWOBIT_DIR! Please move them out of subfolders."
    exit 1
fi

echo "========================================================"
echo " Target Taxa: $TAXA_LIST"
echo " Cores to use: $CORES"
echo "========================================================"

echo " Running LASTZ to find UCEs..."

phyluce_probe_run_multiple_lastzs_sqlite \
    --db "$DB_FILE" \
    --output "$OUTPUT_DIR" \
    --scaffoldlist $TAXA_LIST \
    --genome-base-path "$TWOBIT_DIR" \
    --probefile "$PROBE_FILE" \
    --identity 80 \
    --coverage 80 \
    --cores "$CORES" 

echo "All tasks completed at $(date)"