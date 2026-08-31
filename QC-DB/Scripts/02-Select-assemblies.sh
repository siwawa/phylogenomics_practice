#!/bin/bash
# 02-Select-assemblies.sh
#
# For each taxid (in_refseq_protein == YES, Taxonomy != Asgard_Archaea),
# looks up which assembly each pooled protein belongs to (via NCBI edirect
# IPG lookup) and selects the assembly with the MOST proteins as the
# representative accession for that taxid.
#
# Incremental/resumable: if OUTPUT already has a row for a taxid, that row
# is reused as-is and the IPG lookup is skipped for it -- only taxids not
# yet in OUTPUT are queried. This also means a taxid that failed the IPG
# lookup last run (and so never got a row written) is retried automatically
# on the next run, same as the 01 stage's per-row cache. To force a
# specific taxid to be re-selected, delete its row from OUTPUT (or delete
# the whole file to start over).
#
# Usage:
#   bash 02-Select-assemblies.sh species_info.tsv [output_tsv]
#   (output_tsv defaults to selected_assemblies.tsv next to this script)

set -euo pipefail

INPUT="${1:?Provide path to species_info.tsv as the first argument}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT="${2:-${SCRIPT_DIR}/selected_assemblies.tsv}"
BLASTDB="/bank/ncbi/gnathostome/refseq_protein/refseq_protein"

BLASTDBCMD="/rna/liha/miniconda3/envs/NCBI-download/bin/blastdbcmd"
EPOST="/rna/liha/miniconda3/envs/NCBI-download/bin/epost"
EFETCH="/rna/liha/miniconda3/envs/NCBI-download/bin/efetch"

CHUNK_SIZE=500
SLEEP_BETWEEN=1
MAX_ATTEMPTS=3
RETRY_BASE_SLEEP=5
NEEDS_REVIEW="/rna/liha/phylogenomics_practice/QC-DB/Scripts/needs_review.txt"

source "${SCRIPT_DIR}/lib_ipg_retry.sh"

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

declare -A cached_row
if [[ -s "$OUTPUT" ]]; then
    while IFS= read -r line; do
        cached_taxid="${line%%$'\t'*}"
        [[ "$cached_taxid" =~ ^[0-9]+$ ]] || continue
        cached_row["$cached_taxid"]="$line"
    done < <(tail -n +2 "$OUTPUT")
    echo "[INFO] Loaded ${#cached_row[@]} cached selections from ${OUTPUT}" >&2
fi

TMP_OUT="$WORKDIR/output.tsv"
echo -e "taxid\ttaxonomy\tselected_accession\tn_proteins_in_selected_assembly" > "$TMP_OUT"

tail -n +2 "$INPUT" | while IFS=$'\t' read -r -a fields; do
    taxonomy="${fields[0]}"
    taxid="${fields[3]}"
    in_refseq="${fields[4]}"

    if [[ "$in_refseq" != "YES" || "$taxid" == "NA" ]]; then
        continue
    fi

    if [[ "$taxonomy" == "Asgard_Archaea" ]]; then
        echo "[INFO] Skipping taxid $taxid (Asgard_Archaea, excluded from analysis)" >&2
        continue
    fi

    if [[ -n "${cached_row[$taxid]:-}" ]]; then
        echo "[INFO] taxid=$taxid already selected, reusing cached row." >&2
        echo "${cached_row[$taxid]}" >> "$TMP_OUT"
        continue
    fi

    echo "[INFO] Selecting assembly for taxid=$taxid..." >&2

    protein_file="$WORKDIR/${taxid}_proteins.txt"
    "$BLASTDBCMD" -db "$BLASTDB" -taxids "$taxid" \
        -no_taxid_expansion -target_only \
        -outfmt '%a' | sort -u > "$protein_file"

    n_proteins=$(wc -l < "$protein_file" | tr -d ' ')
    if [[ "$n_proteins" -eq 0 ]]; then
        echo "[WARNING] No proteins found for taxid $taxid; skipping" >&2
        continue
    fi

    assembly_hits="$WORKDIR/${taxid}_assembly_hits.txt"
    : > "$assembly_hits"

    split -l "$CHUNK_SIZE" "$protein_file" "$WORKDIR/${taxid}_chunk_"
    ipg_ok=TRUE
    for chunk in "$WORKDIR/${taxid}_chunk_"*; do
        if ! ipg_lookup_with_retry "$EPOST" "$EFETCH" "$chunk" "$assembly_hits" "$MAX_ATTEMPTS" "$RETRY_BASE_SLEEP"; then
            ipg_ok=FALSE
            printf '%s\t02-Select-assemblies\t%s\t%s\n' \
                "$(date -Iseconds)" "$taxid" "$(basename "$chunk")" >> "$NEEDS_REVIEW"
        fi
        sleep "$SLEEP_BETWEEN"
    done
    rm -f "$WORKDIR/${taxid}_chunk_"*

    if [[ "$ipg_ok" != "TRUE" ]]; then
        echo "[WARNING] IPG lookup incomplete for taxid $taxid (see needs_review.txt); skipping rather than picking from partial data" >&2
        continue
    fi

    if [[ ! -s "$assembly_hits" ]]; then
        echo "[WARNING] IPG lookup returned no assembly hits for taxid $taxid; skipping" >&2
        continue
    fi

    # pick the assembly accession with the most protein hits
    read -r n_hits selected < <(sort "$assembly_hits" | uniq -c | sort -rn | head -1)

    echo -e "${taxid}\t${taxonomy}\t${selected}\t${n_hits}" >> "$TMP_OUT"
done

mv -f "$TMP_OUT" "$OUTPUT"
n_selected=$(($(wc -l < "$OUTPUT") - 1))
echo "[INFO] Wrote ${n_selected} selected assemblies to ${OUTPUT}" >&2
