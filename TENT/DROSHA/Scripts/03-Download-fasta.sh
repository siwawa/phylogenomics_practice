#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HIGH_HITS_TABLE="${PIPELINE_HIGH_HITS_TABLE:-${SCRIPT_DIR}/Blast-high-scoring-hits.txt}"
HIGH_HITS_FASTA="${PIPELINE_HIGH_HITS_FASTA:-${SCRIPT_DIR}/Blast-high-scoring-hits.fasta}"

if [[ ! -s "$HIGH_HITS_TABLE" ]]; then
    echo "Error: high-scoring hit table not found or empty: $HIGH_HITS_TABLE" >&2
    exit 1
fi

ACC_LIST="$(awk 'NR > 1 { print $2 }' "$HIGH_HITS_TABLE" | paste -sd, -)"
if [[ -z "$ACC_LIST" ]]; then
    echo "Error: no accessions found in $HIGH_HITS_TABLE" >&2
    exit 1
fi

curl --fail -sS -d "db=protein&id=${ACC_LIST}&rettype=fasta&retmode=text" \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
    > "$HIGH_HITS_FASTA"
