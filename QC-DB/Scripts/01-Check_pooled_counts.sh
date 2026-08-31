#!/bin/bash
# 01-Check_pooled_counts.sh
#
# Computes pooled-protein/assembly counts for ONE row of species_info.tsv.
# Array-batched: one Slurm array task per data row (see
# tmp-run-01-array.sbatch), rather than one long serial loop over all rows
# -- a single serial run of ~70 taxa against a flaky external API can
# exceed any time limit and, without per-row checkpointing, lose all
# progress at once (see CLAUDE.md).
#
# Writes a single 15-column line (no header) to
# Results/.pooled_counts/<NNN>.tsv, where NNN is the zero-padded row index.
# 01-Collapse-pooled-counts.sh concatenates these into the final
# species_info_pooled_counts.tsv.
#
# Idempotent/resumable: if the output file already exists with a non-NA
# n_assemblies_pooled, this exits immediately without redoing any NCBI
# lookups. Resubmitting the array is how NA rows (transient IPG failures)
# get retried -- there is no separate "retry" script.
#
# Usage: bash 01-Check_pooled_counts.sh <ROW_INDEX>
#   ROW_INDEX is 1-based, matching a data row of species_info.tsv
#   (SLURM_ARRAY_TASK_ID maps directly onto it).

set -euo pipefail

ROW_INDEX="${1:?Provide the 1-based species_info.tsv data-row index as the first argument}"

INPUT="${INPUT:-/rna/liha/phylogenomics_practice/QC-DB/Scripts/species_info.tsv}"
ASSEMBLY_SUMMARY="/rna/liha/phylogenomics_practice/QC-DB/assembly_summary_refseq.txt"
BLASTDB="/bank/ncbi/gnathostome/refseq_protein/refseq_protein"
PER_TAXON_DIR="/rna/liha/phylogenomics_practice/QC-DB/Results/.pooled_counts"

BLASTDBCMD="/rna/liha/miniconda3/envs/NCBI-download/bin/blastdbcmd"
EPOST="/rna/liha/miniconda3/envs/NCBI-download/bin/epost"
EFETCH="/rna/liha/miniconda3/envs/NCBI-download/bin/efetch"

CHUNK_SIZE=500
SLEEP_BETWEEN=1
MAX_ATTEMPTS=3
RETRY_BASE_SLEEP=5

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib_ipg_retry.sh"

mkdir -p "$PER_TAXON_DIR"
OUT_FILE="${PER_TAXON_DIR}/$(printf '%03d' "$ROW_INDEX").tsv"

if [[ -s "$OUT_FILE" ]] && [[ "$(cut -f15 "$OUT_FILE")" != "NA" ]]; then
    echo "[INFO] Row ${ROW_INDEX} already done (${OUT_FILE}), skipping." >&2
    exit 0
fi

line=$(sed -n "$((ROW_INDEX + 1))p" "$INPUT")
[[ -n "$line" ]] || { echo "[ERROR] No data row ${ROW_INDEX} in ${INPUT}" >&2; exit 1; }

IFS=$'\t' read -r -a fields <<< "$line"
taxid="${fields[3]}"
in_refseq="${fields[4]}"

if [[ "$in_refseq" != "YES" || "$taxid" == "NA" ]]; then
    printf '%s\tNA\tNA\tNA\n' "$(IFS=$'\t'; echo "${fields[*]}")" > "$OUT_FILE"
    exit 0
fi

echo "[INFO] Row ${ROW_INDEX}: taxid ${taxid} (${fields[2]})" >&2

WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

n_assemblies=$(awk -F'\t' -v taxid="$taxid" '$6 == taxid {n++} END {print n + 0}' "$ASSEMBLY_SUMMARY")

protein_file="$WORK_DIR/proteins.txt"
"$BLASTDBCMD" \
    -db "$BLASTDB" \
    -taxids "$taxid" \
    -no_taxid_expansion \
    -target_only \
    -outfmt '%a' \
    | sort -u > "$protein_file"

n_proteins=$(wc -l < "$protein_file" | tr -d ' ')

if [[ "$n_proteins" -eq 0 ]]; then
    printf '%s\t%s\t0\t0\n' "$(IFS=$'\t'; echo "${fields[*]}")" "$n_assemblies" > "$OUT_FILE"
    exit 0
fi

assembly_file="$WORK_DIR/assemblies.txt"
: > "$assembly_file"

split -l "$CHUNK_SIZE" "$protein_file" "$WORK_DIR/chunk_"
ipg_ok=TRUE

for chunk in "$WORK_DIR/chunk_"*; do
    if ! ipg_lookup_with_retry "$EPOST" "$EFETCH" "$chunk" "$assembly_file" "$MAX_ATTEMPTS" "$RETRY_BASE_SLEEP"; then
        ipg_ok=FALSE
    fi
    sleep "$SLEEP_BETWEEN"
done

if [[ "$ipg_ok" == "TRUE" ]]; then
    n_pooled_assemblies=$(sort -u "$assembly_file" | awk 'NF {n++} END {print n + 0}')
else
    n_pooled_assemblies=NA
    echo "[WARNING] IPG lookup failed for taxid ${taxid} (row ${ROW_INDEX}); n_assemblies_pooled is NA — will retry if the array is resubmitted" >&2
fi

printf '%s\t%s\t%s\t%s\n' \
    "$(IFS=$'\t'; echo "${fields[*]}")" \
    "$n_assemblies" \
    "$n_proteins" \
    "$n_pooled_assemblies" > "$OUT_FILE"
