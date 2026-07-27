#!/usr/bin/env bash
#SBATCH --job-name=TENT5_RBH
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --exclusive
#SBATCH --nodelist=guanine
#SBATCH --partition=intel
#SBATCH --output=logs/reciprocal_%j.log
#SBATCH --error=logs/reciprocal_%j.log

set -euo pipefail

# Reciprocal best-hit screen for the TENT5 forward-BLAST results.
#
# 1. Collect unique forward hits passing the configured E-value mode.
# 2. Retrieve those protein sequences from NCBI Protein.
# 3. Search them against the local human protein BLAST database.
# 4. Retain a candidate only when its top human hit is an annotated
#    TENT5A, TENT5B, TENT5C, or TENT5D protein.

if [[ -n "${PIPELINE_BASE_DIR:-}" && -d "${PIPELINE_BASE_DIR}/Scripts" ]]; then
    BASE_DIR="$(cd "$PIPELINE_BASE_DIR" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/Scripts" ]]; then
    BASE_DIR="$(cd "$SLURM_SUBMIT_DIR" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" &&
        "$(basename "$SLURM_SUBMIT_DIR")" == "Scripts" &&
        -d "${SLURM_SUBMIT_DIR}/../Blast" ]]; then
    BASE_DIR="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi

SCRIPT_DIR="${BASE_DIR}/Scripts"
RESULTS_DIR="${BASE_DIR}/Blast/Results"
DONE_DIR="${BASE_DIR}/Blast/Done"
SPECIES_FILE="${SCRIPT_DIR}/species.txt"
HUMAN_DB="${HUMAN_DB:-/bank/ncbi/gnathostome/HumanDB/human_prot}"
QUERY_NAME="${PIPELINE_QUERY_NAME:-TENT5A-B-C-D}"
EVALUE_CUTOFF="${RECIPROCAL_EVALUE_CUTOFF:-1e-10}"
NO_EVALUE_FILTER="${RECIPROCAL_NO_EVALUE_FILTER:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
NCBI_EMAIL="${NCBI_EMAIL:-}"
CONDA_ENV="${PIPELINE_BLAST_CONDA_ENV:-NCBI-download}"

if [[ "$NO_EVALUE_FILTER" == "1" ]]; then
    OUTPUT_DIR="${BASE_DIR}/Blast/Reciprocal_no_evalue_filter"
    FORWARD_HITS="${OUTPUT_DIR}/forward_hits_all.tsv"
    RB_HIGH_HITS="${SCRIPT_DIR}/Blast-high-scoring-hits-rb-no-evalue-filter.txt"
    BLAST_EVALUE="10"
    FILTER_DESCRIPTION="without an E-value filter"
else
    OUTPUT_DIR="${BASE_DIR}/Blast/Reciprocal"
    FORWARD_HITS="${OUTPUT_DIR}/forward_hits_evalue_lt_1e-10.tsv"
    RB_HIGH_HITS="${SCRIPT_DIR}/Blast-high-scoring-hits-rb.txt"
    BLAST_EVALUE="$EVALUE_CUTOFF"
    FILTER_DESCRIPTION="with E-value < ${EVALUE_CUTOFF}"
fi

FORWARD_ACCESSIONS="${OUTPUT_DIR}/forward_hit_accessions.txt"
FORWARD_FASTA="${OUTPUT_DIR}/forward_hit_sequences.fasta"
MISSING_ACCESSIONS="${OUTPUT_DIR}/forward_hit_accessions_not_retrieved.txt"
TENT5_ACCESSIONS="${OUTPUT_DIR}/human_TENT5_accessions.tsv"
RECIPROCAL_TOP_HITS="${OUTPUT_DIR}/reciprocal_top_hits_vs_human.tsv"
RECIPROCAL_TENT5_HITS="${OUTPUT_DIR}/reciprocal_TENT5_hits.tsv"
RECIPROCAL_TENT5_FASTA="${OUTPUT_DIR}/reciprocal_TENT5_hit_sequences.fasta"

log() {
    printf '[TENT5 reciprocal BLAST] %s\n' "$*"
}

die() {
    printf '[TENT5 reciprocal BLAST] ERROR: %s\n' "$*" >&2
    exit 1
}

CONDA_SH="${HOME}/miniconda3/etc/profile.d/conda.sh"
[[ -s "$CONDA_SH" ]] || die "Conda initialization script not found: ${CONDA_SH}"
# shellcheck source=/dev/null
source "$CONDA_SH"
conda activate "$CONDA_ENV"

for command_name in blastp blastdbcmd curl awk sort split; do
    command -v "$command_name" >/dev/null 2>&1 ||
        die "Required command not found: ${command_name}"
done

[[ -d "$RESULTS_DIR" ]] || die "Forward-BLAST result directory not found: ${RESULTS_DIR}"
[[ -s "$SPECIES_FILE" ]] || die "Species table not found or empty: ${SPECIES_FILE}"
[[ "$NO_EVALUE_FILTER" == "0" || "$NO_EVALUE_FILTER" == "1" ]] ||
    die "RECIPROCAL_NO_EVALUE_FILTER must be 0 or 1, not '${NO_EVALUE_FILTER}'."
[[ "$EVALUE_CUTOFF" =~ ^[0-9]+([.][0-9]+)?[eE]-[0-9]+$ ]] ||
    die "RECIPROCAL_EVALUE_CUTOFF must look like 1e-10, not '${EVALUE_CUTOFF}'."
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "Invalid thread count: ${THREADS}"

if ! blastdbcmd -db "$HUMAN_DB" -info >/dev/null 2>&1; then
    die "Human BLAST database cannot be opened: ${HUMAN_DB}"
fi

mkdir -p "$OUTPUT_DIR"
WORK_DIR="$(mktemp -d "${OUTPUT_DIR}/.reciprocal_work.XXXXXX")"

cleanup() {
    local status=$?
    trap - EXIT
    if [[ -n "${WORK_DIR:-}" && -d "$WORK_DIR" ]]; then
        rm -rf -- "$WORK_DIR"
    fi
    exit "$status"
}
trap cleanup EXIT

FORWARD_STAGE="${WORK_DIR}/forward_hits.tsv"
ACCESSION_STAGE="${WORK_DIR}/forward_accessions.txt"
FASTA_STAGE="${WORK_DIR}/forward_sequences.fasta"
FETCHED_ACCESSIONS="${WORK_DIR}/fetched_accessions.txt"
MISSING_STAGE="${WORK_DIR}/missing_accessions.txt"
TENT5_STAGE="${WORK_DIR}/human_TENT5_accessions.tsv"
RECIPROCAL_STAGE="${WORK_DIR}/reciprocal_top_hits.tsv"
ACCEPTED_STAGE="${WORK_DIR}/reciprocal_TENT5_hits.tsv"
ACCEPTED_IDS="${WORK_DIR}/accepted_accessions.txt"
ACCEPTED_FASTA_STAGE="${WORK_DIR}/reciprocal_TENT5_sequences.fasta"
RB_HIGH_HITS_STAGE="${WORK_DIR}/Blast-high-scoring-hits-rb.txt"

MISSING_DONE="${WORK_DIR}/missing_done_markers.txt"
awk -F '\t' -v done_dir="$DONE_DIR" '
    NR > 1 && NF {
        n = split($8, words, " ")
        species = words[1]
        if (n >= 2) species = species "_" words[2]
        print done_dir "/" $1 "_" species ".done"
    }
' "$SPECIES_FILE" |
    while IFS= read -r done_file; do
        [[ -f "$done_file" ]] || printf '%s\n' "$done_file"
    done > "$MISSING_DONE"

N_MISSING_DONE="$(awk 'END { print NR + 0 }' "$MISSING_DONE")"
if [[ "$N_MISSING_DONE" -gt 0 ]]; then
    sed -n '1,20p' "$MISSING_DONE" >&2
    log "WARNING: Forward BLAST is incomplete: ${N_MISSING_DONE} expected .done markers are missing. Continuing with available result files."
fi

printf 'qseqid\tsaccver\tpident\tlength\tqcovs\tevalue\tbitscore\tstitle\tsource_file\torganism\tfilepath\n' \
    > "$FORWARD_STAGE"

N_RESULT_FILES=0
while IFS= read -r -d '' result_file; do
    N_RESULT_FILES=$((N_RESULT_FILES + 1))
    source_file="$(basename "$result_file")"
    organism="$(basename "$(dirname "$result_file")")"
    awk -F '\t' -v OFS='\t' -v cutoff="$EVALUE_CUTOFF" \
        -v no_filter="$NO_EVALUE_FILTER" \
        -v source="$source_file" -v organism="$organism" -v filepath="$result_file" \
        'NF >= 8 && (no_filter == 1 || ($6 + 0) < cutoff) {
            print $1, $2, $3, $4, $5, $6, $7, $8, source, organism, filepath
        }' "$result_file" >> "$FORWARD_STAGE"
done < <(
    find "$RESULTS_DIR" -type f -name "*_${QUERY_NAME}.txt" -size +0c -print0
)

[[ "$N_RESULT_FILES" -gt 0 ]] ||
    die "No non-empty forward-BLAST files ending in _${QUERY_NAME}.txt were found."

awk -F '\t' 'NR > 1 && $2 != "" { print $2 }' "$FORWARD_STAGE" |
    sort -u > "$ACCESSION_STAGE"

N_FORWARD_ROWS="$(awk 'END { print (NR > 0 ? NR - 1 : 0) }' "$FORWARD_STAGE")"
N_ACCESSIONS="$(awk 'END { print NR + 0 }' "$ACCESSION_STAGE")"
[[ "$N_ACCESSIONS" -gt 0 ]] ||
    die "No forward hits were found ${FILTER_DESCRIPTION}."

log "Found ${N_FORWARD_ROWS} forward-hit rows ${FILTER_DESCRIPTION}, representing ${N_ACCESSIONS} unique accessions."
log "Retrieving candidate sequences from NCBI Protein in batches of 200."

mkdir -p "${WORK_DIR}/batches"
split -l 200 -d -a 4 "$ACCESSION_STAGE" "${WORK_DIR}/batches/accessions."
: > "$FASTA_STAGE"

for batch_file in "${WORK_DIR}"/batches/accessions.*; do
    batch_ids="$(paste -sd, "$batch_file")"
    curl_args=(
        --fail
        --silent
        --show-error
        --retry 5
        --retry-all-errors
        --data-urlencode "db=protein"
        --data-urlencode "id=${batch_ids}"
        --data-urlencode "rettype=fasta"
        --data-urlencode "retmode=text"
        --data-urlencode "tool=tent5_reciprocal_blast"
    )
    if [[ -n "$NCBI_EMAIL" ]]; then
        curl_args+=(--data-urlencode "email=${NCBI_EMAIL}")
    fi
    curl "${curl_args[@]}" \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
        >> "$FASTA_STAGE"
    sleep 0.4
done

[[ -s "$FASTA_STAGE" ]] || die "NCBI returned no candidate protein sequences."

awk '
    /^>/ {
        id = substr($0, 2)
        sub(/[[:space:]].*$/, "", id)
        n = split(id, fields, "|")
        if (n >= 2 && (fields[1] == "ref" || fields[1] == "gb" ||
                       fields[1] == "emb" || fields[1] == "dbj")) {
            id = fields[2]
        }
        print id
    }
' "$FASTA_STAGE" | sort -u > "$FETCHED_ACCESSIONS"

comm -23 "$ACCESSION_STAGE" "$FETCHED_ACCESSIONS" > "$MISSING_STAGE"
N_MISSING="$(awk 'END { print NR + 0 }' "$MISSING_STAGE")"
if [[ "$N_MISSING" -gt 0 ]]; then
    log "WARNING: NCBI did not return ${N_MISSING} requested accessions; see $(basename "$MISSING_ACCESSIONS")."
fi

log "Identifying all annotated human TENT5 proteins in ${HUMAN_DB}."
blastdbcmd \
    -db "$HUMAN_DB" \
    -entry all \
    -outfmt $'%a\t%t' |
    awk -F '\t' -v OFS='\t' '
        BEGIN { IGNORECASE = 1 }
        $2 ~ /^terminal nucleotidyltransferase 5[ABCD]([[:space:]]|$)/ {
            print $1, $2
        }
    ' > "$TENT5_STAGE"

N_HUMAN_TENT5="$(awk 'END { print NR + 0 }' "$TENT5_STAGE")"
[[ "$N_HUMAN_TENT5" -gt 0 ]] ||
    die "No TENT5A/B/C/D annotations were found in the human BLAST database."
log "Recognized ${N_HUMAN_TENT5} human TENT5 protein records."

log "Running reciprocal BLASTP against the human database with ${THREADS} threads."
blastp \
    -query "$FASTA_STAGE" \
    -db "$HUMAN_DB" \
    -evalue "$BLAST_EVALUE" \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -num_threads "$THREADS" \
    -out "$RECIPROCAL_STAGE" \
    -outfmt '6 qaccver saccver evalue bitscore pident length qcovs stitle'

printf 'candidate_accession\tbest_forward_evalue\treciprocal_top_accession\treciprocal_evalue\tbitscore\tpident\talignment_length\tquery_coverage\treciprocal_subject_title\n' \
    > "$ACCEPTED_STAGE"

awk -F '\t' -v OFS='\t' -v cutoff="$EVALUE_CUTOFF" \
    -v no_filter="$NO_EVALUE_FILTER" \
    -v tent5_file="$TENT5_STAGE" -v forward_file="$FORWARD_STAGE" '
    FILENAME == tent5_file {
        tent5[$1] = 1
        next
    }
    FILENAME == forward_file {
        if (FNR > 1 && (!($2 in forward_best) || ($6 + 0) < forward_best[$2])) {
            forward_best[$2] = $6 + 0
            forward_best_text[$2] = $6
        }
        next
    }
    ($2 in tent5) && (no_filter == 1 || ($3 + 0) < cutoff) {
        print $1, forward_best_text[$1], $2, $3, $4, $5, $6, $7, $8
    }
' "$TENT5_STAGE" "$FORWARD_STAGE" "$RECIPROCAL_STAGE" >> "$ACCEPTED_STAGE"

awk -F '\t' 'NR > 1 { print $1 }' "$ACCEPTED_STAGE" |
    sort -u > "$ACCEPTED_IDS"

awk -F '\t' -v OFS='\t' -v ids_file="$ACCEPTED_IDS" '
    FILENAME == ids_file {
        keep[$1] = 1
        next
    }
    FNR == 1 || ($2 in keep) {
        print
    }
' "$ACCEPTED_IDS" "$FORWARD_STAGE" > "$RB_HIGH_HITS_STAGE"

awk '
    NR == FNR {
        keep[$1] = 1
        next
    }
    /^>/ {
        id = substr($0, 2)
        sub(/[[:space:]].*$/, "", id)
        n = split(id, fields, "|")
        if (n >= 2 && (fields[1] == "ref" || fields[1] == "gb" ||
                       fields[1] == "emb" || fields[1] == "dbj")) {
            id = fields[2]
        }
        emit = (id in keep)
    }
    emit { print }
' "$ACCEPTED_IDS" "$FASTA_STAGE" > "$ACCEPTED_FASTA_STAGE"

mv -f "$FORWARD_STAGE" "$FORWARD_HITS"
mv -f "$ACCESSION_STAGE" "$FORWARD_ACCESSIONS"
mv -f "$FASTA_STAGE" "$FORWARD_FASTA"
mv -f "$MISSING_STAGE" "$MISSING_ACCESSIONS"
mv -f "$TENT5_STAGE" "$TENT5_ACCESSIONS"
mv -f "$RECIPROCAL_STAGE" "$RECIPROCAL_TOP_HITS"
mv -f "$ACCEPTED_STAGE" "$RECIPROCAL_TENT5_HITS"
mv -f "$ACCEPTED_FASTA_STAGE" "$RECIPROCAL_TENT5_FASTA"
mv -f "$RB_HIGH_HITS_STAGE" "$RB_HIGH_HITS"

N_ACCEPTED="$(awk 'END { print (NR > 0 ? NR - 1 : 0) }' "$RECIPROCAL_TENT5_HITS")"
log "Accepted ${N_ACCEPTED} reciprocal TENT5 best hits."
log "Table: ${RECIPROCAL_TENT5_HITS}"
log "FASTA: ${RECIPROCAL_TENT5_FASTA}"
log "RB-filtered forward-hit table: ${RB_HIGH_HITS}"
