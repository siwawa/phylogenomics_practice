#!/usr/bin/env bash

set -euo pipefail

BASE_DIR="/rna/liha/phylogenomics_practice/TENT/TENT2"
INPUT="${BASE_DIR}/BLASTp-results/TENT2_BLASTp_results_isoform_resolved.tsv"
OUTPUT="${BASE_DIR}/BLASTp-results/TENT2_BLASTp_results_isoform_resolved.fasta"
BEST_OUTPUT="${BASE_DIR}/BLASTp-results/TENT2_best_forward_hit_per_species.fasta"
REFSEQ_DB="/bank/ncbi/gnathostome/refseq_protein/refseq_protein"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

source ~/miniconda3/etc/profile.d/conda.sh
conda activate NCBI-download

awk -F '\t' -v OFS='\t' '
NR > 1 && $2 != "NA" {
    species = $12
    gsub(/[[:space:]]+/, "_", species)
    print $2, $11 "_" species "_" $2
}
' "$INPUT" > "$TMP_DIR/accession_headers.tsv"

cut -f1 "$TMP_DIR/accession_headers.tsv" > "$TMP_DIR/accessions.txt"

blastdbcmd \
    -db "$REFSEQ_DB" \
    -entry_batch "$TMP_DIR/accessions.txt" \
    -target_only \
    -outfmt $'%a\t%s' \
    > "$TMP_DIR/sequences.tsv"

awk -F '\t' '
NR == FNR {
    header[$1] = $2
    next
}
$1 in header {
    print ">" header[$1]
    print $2
}
' "$TMP_DIR/accession_headers.tsv" "$TMP_DIR/sequences.tsv" > "$OUTPUT"

# Select one best forward hit per requested species.
awk -F '\t' -v OFS='\t' '
NR > 1 && $2 != "NA" {
    species_key = $12
    accession = $2
    bitscore = $9 + 0
    evalue = $8 + 0
    qcovs = $7 + 0
    pident = $5 + 0
    alignment_length = $6 + 0

    if (!(species_key in seen) ||
        bitscore > best_bitscore[species_key] ||
        (bitscore == best_bitscore[species_key] && evalue < best_evalue[species_key]) ||
        (bitscore == best_bitscore[species_key] && evalue == best_evalue[species_key] && qcovs > best_qcovs[species_key]) ||
        (bitscore == best_bitscore[species_key] && evalue == best_evalue[species_key] && qcovs == best_qcovs[species_key] && pident > best_pident[species_key]) ||
        (bitscore == best_bitscore[species_key] && evalue == best_evalue[species_key] && qcovs == best_qcovs[species_key] && pident == best_pident[species_key] && alignment_length > best_length[species_key]) ||
        (bitscore == best_bitscore[species_key] && evalue == best_evalue[species_key] && qcovs == best_qcovs[species_key] && pident == best_pident[species_key] && alignment_length == best_length[species_key] && accession < best_accession[species_key])) {
        species = $12
        gsub(/[[:space:]]+/, "_", species)
        best_accession[species_key] = accession
        best_header[species_key] = $11 "_" species "_" accession
        best_bitscore[species_key] = bitscore
        best_evalue[species_key] = evalue
        best_qcovs[species_key] = qcovs
        best_pident[species_key] = pident
        best_length[species_key] = alignment_length
        seen[species_key] = 1
    }
}
END {
    for (species_key in seen)
        print best_accession[species_key], best_header[species_key]
}
' "$INPUT" | sort -k2,2 > "$TMP_DIR/best_accession_headers.tsv"

cut -f1 "$TMP_DIR/best_accession_headers.tsv" > "$TMP_DIR/best_accessions.txt"

blastdbcmd \
    -db "$REFSEQ_DB" \
    -entry_batch "$TMP_DIR/best_accessions.txt" \
    -target_only \
    -outfmt $'%a\t%s' \
    > "$TMP_DIR/best_sequences.tsv"

awk -F '\t' '
NR == FNR {
    header[$1] = $2
    next
}
$1 in header {
    print ">" header[$1]
    print $2
}
' "$TMP_DIR/best_accession_headers.tsv" "$TMP_DIR/best_sequences.tsv" > "$BEST_OUTPUT"

echo "FASTA sequences: $(grep -c '^>' "$OUTPUT")"
echo "Output: $OUTPUT"
echo "Best-per-species FASTA sequences: $(grep -c '^>' "$BEST_OUTPUT")"
echo "Best-per-species output: $BEST_OUTPUT"
