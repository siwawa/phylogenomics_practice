#!/bin/bash
# 01-Collapse-pooled-counts.sh
#
# Concatenates the per-row output files from the 01 array
# (Results/.pooled_counts/<NNN>.tsv, one 15-column line each, no header)
# into the final species_info_pooled_counts.tsv, in original
# species_info.tsv row order.
#
# Safe to run at any point, including with the array still in-flight or
# partially failed: missing rows are just counted and reported, not
# treated as fatal by this script -- the caller (tmp-collapse-notify.sbatch)
# decides what counts as an acceptable vs. broken outcome.
#
# Usage: bash 01-Collapse-pooled-counts.sh
# Prints one line to stdout: "<total> <present> <missing> <still_na>"
# (space-separated, for `read`). Human-readable detail goes to stderr.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT="${INPUT:-${SCRIPT_DIR}/species_info.tsv}"
OUTPUT="${OUTPUT:-${SCRIPT_DIR}/species_info_pooled_counts.tsv}"
PER_TAXON_DIR="${PER_TAXON_DIR:-/rna/liha/phylogenomics_practice/QC-DB/Results/.pooled_counts}"

total=$(($(wc -l < "$INPUT") - 1))

TMP_OUT=$(mktemp)
head -n 1 "$INPUT" | awk -F'\t' -v OFS='\t' '{
    print $0, "n_assemblies_exact_taxid", "n_proteins_pooled", "n_assemblies_pooled"
}' > "$TMP_OUT"

present=0
missing=0
missing_rows=()

for i in $(seq 1 "$total"); do
    f="${PER_TAXON_DIR}/$(printf '%03d' "$i").tsv"
    if [[ -s "$f" ]]; then
        cat "$f" >> "$TMP_OUT"
        present=$((present + 1))
    else
        missing=$((missing + 1))
        missing_rows+=("$i")
    fi
done

mv "$TMP_OUT" "$OUTPUT"

still_na=$(awk -F'\t' 'NR>1 && $5=="YES" && $15=="NA"' "$OUTPUT" | wc -l)

echo "[INFO] Collapsed ${present}/${total} rows into ${OUTPUT} (${missing} missing, ${still_na} still NA)" >&2
[[ "$missing" -eq 0 ]] || echo "[INFO] Missing row indices: ${missing_rows[*]}" >&2

echo "${total} ${present} ${missing} ${still_na}"
