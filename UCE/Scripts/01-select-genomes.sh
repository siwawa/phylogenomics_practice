#!/bin/bash

# Download genomes from NCBI 
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate NCBI-download
echo "Job started on $(hostname) at $(date)" 

# 새로 추가된 Target Clade
TARGET_LIST=(
    "Choanoflagellata"   # outgroup
    "Porifera"
    "Ctenophora"
    "Placozoa"
    "Cnidaria"
    "Xenacoelomorpha"   
    "Gnathifera"
    "Platyhelminthes"
    "Lophotrochozoa"
    "Ecdysozoa"
)

# Adopted from /rna/liha/selection_project/Whole-genome/NCBI-select.sh 
CUTOFF_DATE="2015-01-01"
OUTPUT_FILE="phylogenomic_practice.tsv"

# 헤더 작성 (기존 파일 덮어쓰기)
echo -e "Target_Clade\tAccession\tRefSeq_Status\tCoverage\tN50\tGenome_Size\tTaxID\tOrganism_Name" > "$OUTPUT_FILE"

for TAXON in "${TARGET_LIST[@]}"; do

    QUERY_TAXON="$TAXON"
    if [ "$TAXON" == "Cyclostomata" ]; then
        QUERY_TAXON="1476529"
    elif [ "$TAXON" == "Ctenophora" ]; then 
        QUERY_TAXON="10197"
    elif [ "$TAXON" == "Gnathifera" ]; then
        QUERY_TAXON="10190"
    elif [ "$TAXON" == "Cephalochordata" ]; then
        QUERY_TAXON="7737"
    fi


    SAFE_NAME=$(echo "$TAXON" | tr ' ' '_')
    echo "========================================================"
    echo " Processing: $TAXON"
    echo " Strategy: RefSeq first. If missing, take Non-Ref (Cov >= 10x)."
    echo "========================================================"
    
    # 1. 메타데이터 다운로드 및 복합 로직 처리
    datasets summary genome taxon "$QUERY_TAXON" \
        --assembly-level contig,scaffold,chromosome,complete \
        --as-json-lines 2>/dev/null \
    | jq -r -s --arg clade "$TAXON" --arg date "$CUTOFF_DATE" '
        # [Step 1] 기본 품질 필터링
        map(select(
            # Atypical 제외
            (.assembly_info.atypical.is_atypical | not) and
            (
                # reference sequence 혹은 representative sequence일 시 조건에 구애받지 않고 통과
                (.assembly_info.refseq_category == "reference genome") or
                (.assembly_info.refseq_category == "representative genome") or
                # reference sequence 혹은 representative sequence가 아닐 시 2015년 이후 + Coverage 10x 이상일 경우에만 통과 
                (
                    ((.assembly_info.release_date // "1900-01-01") >= $date) and 
                    (( (.assembly_stats.genome_coverage // "0") | sub("x";"";"i") | tonumber? // 0) >= 10)
                )
            )
        ))
        
        | sort_by(
            # 1. Reference/Representative면 무조건 최우선
            (if (.assembly_info.refseq_category == "reference genome" or 
                 .assembly_info.refseq_category == "representative genome") then 0 else 1 end),

            # 2. Reference/Representative일 시 GCF로 시작하는 파일을 선정
            (if (.assembly_info.refseq_category == "reference genome" or 
                 .assembly_info.refseq_category == "representative genome") then 
                 (if (.accession | startswith("GCF_")) then 0 else 1 end) else 1 end),
            
            # 3. 만약 Reference/Representative이 아니라면 최신의 subission을 선정 
            ((.assembly_info.release_date // "0") | sub("-";"";"g") | tonumber? * -1)
        )
          
        # [Step 3] 1종 1게놈 선택
        | group_by(.organism.tax_id)
        | map(.[0])

        # [Step 4] 출력
        | .[] | "\($clade)\t\(.accession)\t\(.assembly_info.refseq_category // "na")\t\(.assembly_stats.genome_coverage // "na")\t\(.assembly_stats.contig_n50 // "0")\t\(.assembly_stats.total_sequence_length // "0")\t\(.organism.tax_id)\t\(.organism.organism_name // .organism.common_name // "Unknown")"
    ' >> "$OUTPUT_FILE"
done

echo "All tasks completed at $(date)"