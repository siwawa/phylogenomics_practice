#!/usr/bin/env bash
#SBATCH --job-name=TENT2_RBH
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --nodelist=pichu
#SBATCH --output=logs/TENT2_reciprocal_%j.log
#SBATCH --error=logs/TENT2_reciprocal_%j.log

set -euo pipefail

BASE_DIR="/rna/liha/phylogenomics_practice/TENT/TENT2"
FORWARD_HITS="${BASE_DIR}/BLASTp-results/TENT2_BLASTp_results_isoform_resolved.tsv"
HUMAN_TENT_FASTA="${BASE_DIR}/TENT.fasta"
REFSEQ_DB="/bank/ncbi/gnathostome/refseq_protein/refseq_protein"
HUMAN_DB="/bank/ncbi/gnathostome/HumanDB/human_prot"
OUTPUT_DIR="${BASE_DIR}/BLAST-reciprocal-hits"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

HUMAN_TENT_REFERENCE="${OUTPUT_DIR}/human_TENT_reference_blast.tsv"
HUMAN_TENT_ACCESSIONS="${OUTPUT_DIR}/human_TENT_accessions.txt"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

source ~/miniconda3/etc/profile.d/conda.sh
conda activate NCBI-download

mkdir -p "$OUTPUT_DIR"

# Find the HumanDB accessions corresponding to the known human TENT proteins.
blastp \
    -query "$HUMAN_TENT_FASTA" \
    -db "$HUMAN_DB" \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -num_threads "$THREADS" \
    -outfmt "6 qseqid saccver pident qcovs evalue bitscore stitle" \
    -out "$HUMAN_TENT_REFERENCE"

awk -F '\t' -v OFS='\t' '
$4 == 100 {
    tent_name = $1
    sub(/_Homo_sapiens$/, "", tent_name)
    print $2, tent_name
}
' "$HUMAN_TENT_REFERENCE" |
    sort -u > "$TMP_DIR/human_TENT_seed_accessions.tsv"

# Expand every seed TENT protein to all RefSeq protein isoforms from its GeneID.
: > "$TMP_DIR/human_TENT_all_accessions.tsv"
while IFS=$'\t' read -r seed_accession tent_name; do
    gene_id=$(curl -L -sS --fail --retry 4 --retry-all-errors --retry-delay 2 -G "$EUTILS/elink.fcgi" \
        --data-urlencode "dbfrom=protein" \
        --data-urlencode "db=gene" \
        --data-urlencode "id=$seed_accession" \
        --data-urlencode "linkname=protein_gene" \
        --data-urlencode "retmode=json" |
        jq -r '.linksets[0].linksetdbs[]? | select(.linkname == "protein_gene") | .links[]?' |
        sort -u)

    [[ $(printf '%s\n' "$gene_id" | sed '/^$/d' | wc -l) -eq 1 ]] || continue
    sleep 0.5

    protein_uids=$(curl -L -sS --fail --retry 4 --retry-all-errors --retry-delay 2 -G "$EUTILS/elink.fcgi" \
        --data-urlencode "dbfrom=gene" \
        --data-urlencode "db=protein" \
        --data-urlencode "id=$gene_id" \
        --data-urlencode "linkname=gene_protein_refseq" \
        --data-urlencode "retmode=json" |
        jq -r '.linksets[0].linksetdbs[]? | select(.linkname == "gene_protein_refseq") | .links[]?' |
        paste -sd,)

    [[ -n "$protein_uids" ]] || continue
    sleep 0.5

    curl -L -sS --fail --retry 4 --retry-all-errors --retry-delay 2 -G "$EUTILS/esummary.fcgi" \
        --data-urlencode "db=protein" \
        --data-urlencode "id=$protein_uids" \
        --data-urlencode "retmode=json" |
        jq -r --arg tent_name "$tent_name" '.result.uids[] as $uid | [.result[$uid].accessionversion, $tent_name] | @tsv' \
        >> "$TMP_DIR/human_TENT_all_accessions.tsv"
    sleep 0.5
done < "$TMP_DIR/human_TENT_seed_accessions.tsv"

# Retain only TENT isoform accessions that occur in the HumanDB BLAST database.
cut -f1 "$TMP_DIR/human_TENT_all_accessions.tsv" | sort -u > "$TMP_DIR/human_TENT_all_accession_ids.txt"
blastdbcmd \
    -db "$HUMAN_DB" \
    -entry_batch "$TMP_DIR/human_TENT_all_accession_ids.txt" \
    -target_only \
    -outfmt "%a" \
    > "$TMP_DIR/human_TENT_in_HumanDB.txt" 2> "$TMP_DIR/human_TENT_absent_from_HumanDB.log" || true

awk -F '\t' -v OFS='\t' '
NR == FNR { in_human_db[$1] = 1; next }
$1 in in_human_db { print }
' "$TMP_DIR/human_TENT_in_HumanDB.txt" "$TMP_DIR/human_TENT_all_accessions.tsv" |
    sort -u > "$HUMAN_TENT_ACCESSIONS"

# Search every retained forward hit separately against HumanDB.
awk -F '\t' 'NR > 1 && $2 != "NA" {print $2 "\t" $11}' "$FORWARD_HITS" |
while IFS=$'\t' read -r accession taxonomy; do
    query_name="${accession}_${taxonomy}"
    query_dir="${OUTPUT_DIR}/${query_name}"
    query_fasta="${TMP_DIR}/${query_name}.fasta"
    raw_hits="${TMP_DIR}/${query_name}_raw.tsv"
    reciprocal_hits="${query_dir}/${query_name}_reciprocal_blast_hit.txt"

    mkdir -p "$query_dir"

    blastdbcmd \
        -db "$REFSEQ_DB" \
        -entry "$accession" \
        -target_only \
        -outfmt "%f" \
        > "$query_fasta"

    blastp \
        -query "$query_fasta" \
        -db "$HUMAN_DB" \
        -evalue 0.1 \
        -max_target_seqs 100 \
        -max_hsps 1 \
        -num_threads "$THREADS" \
        -outfmt "6 qseqid saccver staxids sscinames pident length qcovs evalue bitscore stitle" \
        -out "$raw_hits"

    {
        printf 'qseqid\tsaccver\tstaxids\tsscinames\tpident\tlength\tqcovs\tevalue\tbitscore\tstitle\tis_human_TENT\thuman_TENT\n'
        awk -F '\t' -v OFS='\t' '
        NR == FNR {
            tent[$1] = ($1 in tent ? tent[$1] "," $2 : $2)
            next
        }
        { print $0, ($2 in tent ? "TRUE" : "FALSE"), ($2 in tent ? tent[$2] : "NA") }
        ' "$HUMAN_TENT_ACCESSIONS" "$raw_hits"
    } > "$reciprocal_hits"
done

echo "Human TENT reference accessions: $(wc -l < "$HUMAN_TENT_ACCESSIONS")"
echo "Reciprocal-hit directories: $OUTPUT_DIR"
