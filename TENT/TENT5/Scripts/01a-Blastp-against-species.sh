#!/bin/bash
#SBATCH --job-name=PhyloBlastp
#SBATCH --output=Blast/logs/%A_%a.out
#SBATCH --error=Blast/logs/%A_%a.err
#SBATCH --array=2-708%3
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodelist=guanine
#SBATCH --partition=intel
set -euo pipefail


# Downloads one assembly proteome, searches it locally, and removes the temporary files.
usage() {
    cat <<'USAGE'
Usage: sbatch 01a-Blastp-against-species.sh [options]

For each array task, download the annotated proteome for the assembly accession
in species.txt, run BLASTP locally against that proteome, and delete the
temporary download. Output paths and columns match the original workflow.

Options:
  --tmp-root PATH  Place per-task temporary directories under PATH.
                   Default: SLURM_TMPDIR when set, otherwise /tmp.
  --force          Rerun and overwrite a result even if its .done marker exists.
  -h, --help       Show this help.

Environment equivalents:
  PIPELINE_TMP_ROOT=PATH   Override the temporary-storage root.
  PIPELINE_BLAST_FORCE=1   Rerun completed array tasks.
USAGE
}

TMP_ROOT="${PIPELINE_TMP_ROOT:-${SLURM_TMPDIR:-/tmp}}"
BLAST_FORCE="${PIPELINE_BLAST_FORCE:-0}"

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --tmp-root)
            [[ "$#" -ge 2 ]] || { echo "Error: --tmp-root requires a directory path." >&2; exit 2; }
            TMP_ROOT="$2"
            shift 2
            ;;
        --force)
            BLAST_FORCE=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ "$BLAST_FORCE" != "0" && "$BLAST_FORCE" != "1" ]]; then
    echo "Error: PIPELINE_BLAST_FORCE must be 0 or 1, not '${BLAST_FORCE}'." >&2
    exit 2
fi

if [[ -n "${PIPELINE_BASE_DIR:-}" && -d "${PIPELINE_BASE_DIR}/Scripts" ]]; then
    BASE_DIR="$(cd "$PIPELINE_BASE_DIR" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/Scripts" ]]; then
    BASE_DIR="$(cd "$SLURM_SUBMIT_DIR" && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi
SCRIPT_DIR="${BASE_DIR}/Scripts"
PIPELINE_QUERY_NAME="${PIPELINE_QUERY_NAME:-SLIT1-2-3}"
PIPELINE_MAIL_TO="${PIPELINE_MAIL_TO:-}"
PIPELINE_MAIL_SCRIPT="${PIPELINE_MAIL_SCRIPT:-}"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate NCBI-download

for REQUIRED_COMMAND in datasets unzip blastp; do
    if ! command -v "$REQUIRED_COMMAND" >/dev/null 2>&1; then
        echo "Error: required command not found in NCBI-download environment: ${REQUIRED_COMMAND}" >&2
        exit 1
    fi
done

SPECIES_FILE="${SCRIPT_DIR}/species.txt"
QUERY_DIR="${BASE_DIR}/Blast/Query"
OUTPUT_DIR="${BASE_DIR}/Blast/Results"
DONE_DIR="${BASE_DIR}/Blast/Done"

mkdir -p "$OUTPUT_DIR" "$DONE_DIR"

ROW="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_FILE")"
if [[ -z "$ROW" ]]; then
    echo "No species row found for array task ${SLURM_ARRAY_TASK_ID}."
    exit 1
fi

TARGET_CLADE="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $1}')"
ACCESSION="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $2}')"
TAXID="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $7}')"
SPECIES="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $8}')"
SPECIES_LABEL="$(printf '%s\n' "$SPECIES" | awk '{ if (NF >= 2) print $1, $2; else print $0 }')"

SPECIES_SAFE="${TARGET_CLADE}_$(echo "$SPECIES_LABEL" | tr ' ' '_')"
SPECIES_DIR="${OUTPUT_DIR}/${SPECIES_SAFE}"
DONE_FILE="${DONE_DIR}/${SPECIES_SAFE}.done"
QUERY="${QUERY_DIR}/${PIPELINE_QUERY_NAME}.fasta"
RESULT_FILE="${SPECIES_DIR}/${TARGET_CLADE}_${PIPELINE_QUERY_NAME}.txt"

mkdir -p "$SPECIES_DIR"

if [[ ! -s "$QUERY" ]]; then
    echo "Error: query FASTA not found or empty: $QUERY"
    exit 1
fi

if [[ -f "$DONE_FILE" && "$BLAST_FORCE" -eq 0 ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} | ${SPECIES_SAFE} already completed. Skipping BLAST."
else
    if [[ -z "$ACCESSION" ]]; then
        echo "Error: assembly accession is empty for array task ${SLURM_ARRAY_TASK_ID}." >&2
        exit 1
    fi
    if [[ ! -d "$TMP_ROOT" || ! -w "$TMP_ROOT" ]]; then
        echo "Error: temporary-storage root does not exist or is not writable: ${TMP_ROOT}" >&2
        exit 1
    fi

    TASK_TMP=""
    RESULT_STAGE=""
    cleanup_task_files() {
        local status=$?
        trap - EXIT
        if [[ -n "$RESULT_STAGE" && -e "$RESULT_STAGE" ]]; then
            rm -f -- "$RESULT_STAGE"
        fi
        if [[ -n "$TASK_TMP" && -d "$TASK_TMP" ]]; then
            echo "Removing temporary task directory: ${TASK_TMP}"
            rm -rf -- "$TASK_TMP"
        fi
        exit "$status"
    }
    trap cleanup_task_files EXIT

    TASK_TMP="$(mktemp -d "${TMP_ROOT%/}/tent5_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}.XXXXXX")"
    DOWNLOAD_ZIP="${TASK_TMP}/ncbi_dataset.zip"
    PACKAGE_DIR="${TASK_TMP}/package"
    TASK_RESULT="${TASK_TMP}/blast_result.txt"

    echo "Array task ${SLURM_ARRAY_TASK_ID} | accession: ${ACCESSION} | taxid: ${TAXID} | species: ${SPECIES_SAFE}"
    echo "Running BLASTP: ${PIPELINE_QUERY_NAME} x ${SPECIES_SAFE}"
    echo "Temporary task directory: ${TASK_TMP}"

    echo "Downloading annotated proteome for ${ACCESSION}"
    datasets download genome accession "$ACCESSION" \
        --include protein \
        --filename "$DOWNLOAD_ZIP" \
        --no-progressbar

    mkdir -p "$PACKAGE_DIR"
    unzip -q "$DOWNLOAD_ZIP" -d "$PACKAGE_DIR"
    PROTEIN_FASTA="$(find "$PACKAGE_DIR" -type f -name 'protein.faa' -size +0c -print -quit)"
    if [[ -z "$PROTEIN_FASTA" ]]; then
        echo "Error: NCBI supplied no annotated protein FASTA for assembly ${ACCESSION}." >&2
        exit 1
    fi
    echo "Proteome: ${PROTEIN_FASTA}"

    # Output format: QueryID, TargetID, PID, Alignment length, query coverage, E-value, Bitscore, Subject title
    blastp \
        -query "$QUERY" \
        -subject "$PROTEIN_FASTA" \
        -parse_deflines \
        -out "$TASK_RESULT" \
        -outfmt "6 qseqid saccver pident length qcovs evalue bitscore stitle" \
        -max_target_seqs 10

    # Copy through a staging file so a failed task cannot replace a valid result.
    RESULT_STAGE="${SPECIES_DIR}/.${TARGET_CLADE}_${PIPELINE_QUERY_NAME}.${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}.tmp"
    cp "$TASK_RESULT" "$RESULT_STAGE"
    mv -f "$RESULT_STAGE" "$RESULT_FILE"
    RESULT_STAGE=""

    echo "Done: ${PIPELINE_QUERY_NAME} x ${SPECIES_SAFE}"
    touch "$DONE_FILE"
fi

send_finished_mail() {
    [[ -n "$PIPELINE_MAIL_TO" && -s "$PIPELINE_MAIL_SCRIPT" ]] || return 0

    local body mail_status
    body="BLASTp search finished for ${PIPELINE_QUERY_NAME} in ${BASE_DIR}."

    set +e
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate mail
    python3 "$PIPELINE_MAIL_SCRIPT" \
        --to "$PIPELINE_MAIL_TO" \
        --subject "BLASTp search finished: ${PIPELINE_QUERY_NAME}" \
        --body "$body"
    mail_status=$?
    conda deactivate
    set -e

    if [[ "$mail_status" -ne 0 ]]; then
        echo "Warning: mail notification failed with status ${mail_status}."
    fi
}

N_DONE="$(find "$DONE_DIR" -name "*.done" | wc -l | tr -d ' ')"
N_TOTAL="$(awk 'NR > 1 && NF { n++ } END { print n + 0 }' "$SPECIES_FILE")"

if [[ "$N_TOTAL" -gt 0 && "$N_DONE" -ge "$N_TOTAL" ]]; then
    send_finished_mail
    echo "Done at $(date)"
fi
