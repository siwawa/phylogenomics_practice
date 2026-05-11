#!/bin/bash

# NCBI 데이터셋 요약 정보 추출용 스크립트
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate NCBI-download
echo "Job started on $(hostname) at $(date)" 

# 출력 경로 및 설정
OUTPUT_FILE="/rna/liha/phylogenomics_practice/UCE/harvest_UCE/00-species/gnathostome_reference.tsv"
CUTOFF_DATE="2015-01-01"
TARGET_CLADE="Gnathostomes"

# 분석 대상 13종 리스트 (추천드린 외군 포함)
SPECIES_LIST=(
    "Homo sapiens"
    "Mus musculus"
    "Gallus gallus"
    "Xenopus tropicalis"
    "Tetraodon nigroviridis"
    "Lepisosteus oculatus"
    "Scyliorhinus torazame"
    "Callorhinchus milii"
    "Chrysemys picta bellii"
    "Latimeria chalumnae"
    "Polypterus senegalus"
    "Acipenser ruthenus"
    "Amia calva"
)

# 헤더 작성
echo -e "Target_Clade\tAccession\tRefSeq_Status\tCoverage\tN50\tGenome_Size\tTaxID\tOrganism_Name" > "$OUTPUT_FILE"

for SPECIES in "${SPECIES_LIST[@]}"; do
    echo "--------------------------------------------------------"
    echo " Searching metadata for: $SPECIES"
    
    # datasets summary를 통해 메타데이터 확보 후 기존 jq 로직으로 필터링
    datasets summary genome taxon "$SPECIES" \
        --assembly-level contig,scaffold,chromosome,complete \
        --as-json-lines 2>/dev/null \
    | jq -r -s --arg clade "$TARGET_CLADE" --arg date "$CUTOFF_DATE" '
        # [Step 1] 기본 품질 및 기준 필터링
        map(select(
            (.assembly_info.atypical.is_atypical | not) and
            (
                (.assembly_info.refseq_category == "reference genome") or
                (.assembly_info.refseq_category == "representative genome") or
                (
                    ((.assembly_info.release_date // "1900-01-01") >= $date) and 
                    (( (.assembly_stats.genome_coverage // "0") | sub("x";"";"i") | tonumber? // 0) >= 10)
                )
            )
        ))
        
        # [Step 2] 우선순위 정렬 (RefSeq 우선 -> GCF 우선 -> 최신순)
        | sort_by(
            (if (.assembly_info.refseq_category == "reference genome" or 
                 .assembly_info.refseq_category == "representative genome") then 0 else 1 end),
            (if (.assembly_info.refseq_category == "reference genome" or 
                 .assembly_info.refseq_category == "representative genome") then 
                 (if (.accession | startswith("GCF_")) then 0 else 1 end) else 1 end),
            ((.assembly_info.release_date // "0") | sub("-";"";"g") | tonumber? * -1)
        )
          
        # [Step 3] 가장 적합한 1개만 선택
        | .[0:1] 

        # [Step 4] TSV 포맷 출력
        | .[] | "\($clade)\t\(.accession)\t\(.assembly_info.refseq_category // "na")\t\(.assembly_stats.genome_coverage // "na")\t\(.assembly_stats.contig_n50 // "0")\t\(.assembly_stats.total_sequence_length // "0")\t\(.organism.tax_id)\t\(.organism.organism_name // .organism.common_name // "Unknown")"
    ' >> "$OUTPUT_FILE"
done

echo "--------------------------------------------------------"
echo "Metadata extraction for Gnathostomes completed!"
echo "Check output: $OUTPUT_FILE"