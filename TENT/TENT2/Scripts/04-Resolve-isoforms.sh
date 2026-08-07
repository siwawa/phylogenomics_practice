#!/usr/bin/env bash

set -euo pipefail

BASE_DIR="/rna/liha/phylogenomics_practice/TENT/TENT2"
INPUT="${BASE_DIR}/BLASTp-results/TENT2_BLASTp_results.tsv"
OUTPUT="${BASE_DIR}/BLASTp-results/TENT2_BLASTp_results_isoform_resolved.tsv"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

awk -F '\t' 'NR > 1 && $2 != "NA" {print $2}' "$INPUT" | sort -u > "$TMP_DIR/accessions.txt"
split -l 100 -d "$TMP_DIR/accessions.txt" "$TMP_DIR/accessions_"

: > "$TMP_DIR/accession_gene.tsv"

for accession_chunk in "$TMP_DIR"/accessions_*; do
    query=$(awk 'BEGIN { ORS = "" } { if (NR > 1) printf " OR "; printf "%s[ACCN]", $0 }' "$accession_chunk")

    uids=$(curl -L -sS --fail --retry 4 --retry-all-errors --retry-delay 2 -G "$EUTILS/esearch.fcgi" \
        --data-urlencode "db=protein" \
        --data-urlencode "term=$query" \
        --data-urlencode "retmax=100" \
        --data-urlencode "retmode=json" | jq -r '.esearchresult.idlist[]' | paste -sd,)

    [[ -n "$uids" ]] || continue
    sleep 0.5

    curl -L -sS --fail --retry 4 --retry-all-errors --retry-delay 2 -G "$EUTILS/efetch.fcgi" \
        --data-urlencode "db=protein" \
        --data-urlencode "id=$uids" \
        --data-urlencode "rettype=gb" \
        --data-urlencode "retmode=text" |
        awk '
        /^VERSION/ { accession = $2 }
        /\/db_xref="GeneID:/ {
            gene = $0
            sub(/^.*GeneID:/, "", gene)
            sub(/".*$/, "", gene)
            genes[gene] = 1
        }
        /^\/\/$/ {
            count = 0
            for (gene in genes) {
                only_gene = gene
                count++
            }
            if (accession != "" && count == 1) print accession "\t" only_gene
            delete genes
            accession = ""
        }
        ' >> "$TMP_DIR/accession_gene.tsv"
    sleep 0.5
done

# Keep only accessions with exactly one linked NCBI GeneID.
sort -u "$TMP_DIR/accession_gene.tsv" |
    awk -F '\t' -v OFS='\t' '
    $1 != previous {
        if (count == 1) print previous, gene
        previous = $1
        gene = $2
        count = 1
        next
    }
    { count++ }
    END { if (count == 1) print previous, gene }
    ' > "$TMP_DIR/accession_gene_unique.tsv"

# Keep the first AWK input non-empty even if no accession has a GeneID.
printf '__no_accession__\t__no_gene__\n' >> "$TMP_DIR/accession_gene_unique.tsv"

awk -F '\t' -v OFS='\t' '
NR == FNR {
    gene[$1] = $2
    next
}

FNR == 1 {
    print $0, "Isoform_resolved"
    next
}

{
    row_number = FNR
    resolved = ($2 in gene)
    key = resolved ? $3 SUBSEP gene[$2] : "unresolved" SUBSEP row_number

    if (!(key in best) ||
        $9 + 0 > best_bitscore[key] ||
        ($9 + 0 == best_bitscore[key] && $8 + 0 < best_evalue[key]) ||
        ($9 + 0 == best_bitscore[key] && $8 + 0 == best_evalue[key] && $7 + 0 > best_qcovs[key]) ||
        ($9 + 0 == best_bitscore[key] && $8 + 0 == best_evalue[key] && $7 + 0 == best_qcovs[key] && $6 + 0 > best_length[key])) {
        if (key in best) keep[best[key]] = 0
        best[key] = row_number
        best_bitscore[key] = $9 + 0
        best_evalue[key] = $8 + 0
        best_qcovs[key] = $7 + 0
        best_length[key] = $6 + 0
        line[row_number] = $0
        flag[row_number] = resolved ? "TRUE" : "FALSE"
        keep[row_number] = 1
    }
}

END {
    for (row_number = 2; row_number <= FNR; row_number++)
        if (keep[row_number]) print line[row_number], flag[row_number]
}
' "$TMP_DIR/accession_gene_unique.tsv" "$INPUT" > "$OUTPUT"
