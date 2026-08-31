#!/bin/bash
# lib_ipg_retry.sh
#
# Shared retry helper for epost|efetch IPG lookups, used by
# 01-Check_pooled_counts.sh and 02-Select-assemblies.sh.
#
# NCBI eutils intermittently returns 502s. edirect already retries once
# internally (visible as "WARNING: FAILURE" + "SECOND ATTEMPT" in stderr),
# but can still fail twice in a row. This wraps a whole chunk lookup in
# extra attempts with linear backoff before giving up on it.

# ipg_lookup_with_retry <epost_bin> <efetch_bin> <chunk_file> <outfile> <max_attempts> <base_sleep_seconds>
#
# On success: appends efetch's IPG hits (already awk-filtered to GCA_/GCF_
# accessions) to outfile, returns 0.
# On exhausted retries: returns 1. outfile may contain a partial/empty
# result for this chunk (harmless — only ever well-formed accession lines,
# deduped by the caller via `sort -u`).
ipg_lookup_with_retry() {
    local epost_bin="$1" efetch_bin="$2" chunk="$3" outfile="$4"
    local max_attempts="$5" base_sleep="$6"
    local attempt=1

    while (( attempt <= max_attempts )); do
        if "$epost_bin" -db protein -input "$chunk" -format acc \
            | "$efetch_bin" -format ipg 2>/dev/null \
            | awk -F'\t' '$NF ~ /^GC[AF]_/ {print $NF}' \
            >> "$outfile"; then
            return 0
        fi

        echo "[WARNING] IPG lookup failed (attempt ${attempt}/${max_attempts}) for chunk $(basename "$chunk"); retrying in $((base_sleep * attempt))s..." >&2
        sleep "$((base_sleep * attempt))"
        (( attempt++ ))
    done

    echo "[ERROR] IPG lookup permanently failed for chunk $(basename "$chunk") after ${max_attempts} attempts" >&2
    return 1
}
