#!/usr/bin/env bash
#SBATCH --job-name=REFSEQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --nodelist=pichu
#SBATCH --output=logs/REFSEQ_%j.log
#SBATCH --error=logs/REFSEQ_%j.log



set -euo pipefail

SCRIPT_DIR=${CHECK_REFSEQ_WORK_DIR:-/rna/liha/phylogenomics_practice/TENT/TENT2/Scripts}
DB_DIR=/bank/ncbi/gnathostome/refseq_protein
DB="$DB_DIR/refseq_protein"
NAMES_DMP="$DB_DIR/taxonomy/names.dmp"
BLASTDBCMD=/rna/liha/miniconda3/envs/NCBI-download/bin/blastdbcmd

FINAL_SUMMARY="$SCRIPT_DIR/refseq_protein_all_taxid_summary.tsv"

if [[ ! -x "$BLASTDBCMD" ]]; then
    echo "ERROR: blastdbcmd is not executable: $BLASTDBCMD" >&2
    exit 1
fi

if [[ ! -r "$NAMES_DMP" ]]; then
    echo "ERROR: NCBI taxonomy names file is not readable: $NAMES_DMP" >&2
    exit 1
fi

mkdir -p "$SCRIPT_DIR"
export BLASTDB="$DB_DIR"
export LC_ALL=C

TMP_COUNTS=$(mktemp "$SCRIPT_DIR/.refseq_counts.XXXXXX")
TMP_FINAL=$(mktemp "$SCRIPT_DIR/.refseq_summary.XXXXXX")
trap 'rm -f -- "$TMP_COUNTS" "$TMP_FINAL"' EXIT

echo "[1/2] Scanning all records in refseq_protein..."

"$BLASTDBCMD" \
    -db "$DB" \
    -entry all \
    -target_only \
    -outfmt $'%T\t%a' |
awk -F '\t' '
BEGIN {
    OFS = "\t"
}

NF >= 2 && $1 ~ /^[0-9]+$/ {
    taxid = $1
    accession = $2
    prefix = "OTHER"

    # Example: XP_012345.1 -> XP_
    if (match(accession, /^[[:alpha:]]+_/))
        prefix = substr(accession, RSTART, RLENGTH)

    total[taxid]++
    count[taxid SUBSEP prefix]++
}

END {
    for (taxid in total) {
        np = count[taxid SUBSEP "NP_"] + 0
        xp = count[taxid SUBSEP "XP_"] + 0
        wp = count[taxid SUBSEP "WP_"] + 0
        yp = count[taxid SUBSEP "YP_"] + 0
        other = total[taxid] - np - xp - wp - yp
        ratio = sprintf("%.4f", xp / total[taxid])

        print taxid, total[taxid], np, xp, wp, yp, other, ratio
    }
}
' | sort -k1,1n > "$TMP_COUNTS"

echo "[2/2] Adding scientific names and writing the final summary..."

awk -F '\t' -v OFS='\t' '
NR == FNR {
    split($0, field, "\t\\|\t")
    if (field[4] ~ /^scientific name/)
        scientific_name[field[1]] = field[2]
    next
}

{
    name = scientific_name[$1]
    if (name == "") name = "(name unavailable)"

    print name, $1, "YES", $2, $3, $4, $5, $6, $7, $8
}
' "$NAMES_DMP" "$TMP_COUNTS" > "$TMP_FINAL"

{
    printf 'species\ttaxid\tin_refseq_protein\ttotal\tNP_\tXP_\tWP_\tYP_\tOTHER\tXP_ratio\n'
    cat "$TMP_FINAL"
} > "$FINAL_SUMMARY"

rm -f -- "$TMP_COUNTS" "$TMP_FINAL"
trap - EXIT

echo "Done: $FINAL_SUMMARY"
