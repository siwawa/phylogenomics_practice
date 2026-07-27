#!/bin/bash
#SBATCH --job-name=TENT5BlastRescue
#SBATCH --output=Blast/logs/rescue_%A_%a.out
#SBATCH --error=Blast/logs/rescue_%A_%a.err
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
set -euo pipefail

# Emergency remote-BLAST accelerator for the 12 Ecdysozoa and annelid targets.
# Three array workers share the targets round-robin (4 targets per worker).
# Results and .done markers intentionally match 01-Blastp-against-species.sh.

if [[ -n "${PIPELINE_BASE_DIR:-}" && -d "${PIPELINE_BASE_DIR}/Scripts" ]]; then
    BASE_DIR="$(cd "$PIPELINE_BASE_DIR" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/Scripts" ]]; then
    BASE_DIR="$(cd "$SLURM_SUBMIT_DIR" && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi

SCRIPT_DIR="${BASE_DIR}/Scripts"
SPECIES_FILE="${SCRIPT_DIR}/species.txt"
QUERY_DIR="${BASE_DIR}/Blast/Query"
OUTPUT_DIR="${BASE_DIR}/Blast/Results"
DONE_DIR="${BASE_DIR}/Blast/Done"
PIPELINE_QUERY_NAME="${PIPELINE_QUERY_NAME:-TENT5A-B-C-D}"
QUERY="${QUERY_DIR}/${PIPELINE_QUERY_NAME}.fasta"

N_WORKERS=3
EXPECTED_TARGETS=12
WORKER_ID="${SLURM_ARRAY_TASK_ID:-}"

if [[ ! "$WORKER_ID" =~ ^[0-2]$ ]]; then
    echo "Error: this script must run as SLURM array task 0-2; got '${WORKER_ID:-unset}'." >&2
    exit 2
fi
if [[ ! -s "$SPECIES_FILE" ]]; then
    echo "Error: species table not found or empty: $SPECIES_FILE" >&2
    exit 1
fi
if [[ ! -s "$QUERY" ]]; then
    echo "Error: query FASTA not found or empty: $QUERY" >&2
    exit 1
fi

source ~/miniconda3/etc/profile.d/conda.sh
conda activate NCBI-download

if ! command -v blastp >/dev/null 2>&1; then
    echo "Error: blastp is unavailable in the NCBI-download environment." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR" "$DONE_DIR"

mapfile -t TARGET_ROWS < <(
    awk -F '\t' '
        NR > 1 {
            if (($1 ~ /^Ecdysozoa_/ && $1 !~ /^Ecdysozoa_Arthropoda/) || $1 ~ /^Ecdysozoa_Arthropoda/ || $1 ~ /^Lophotrochozoa_Annelida/) {
                print
            }
        }
    ' "$SPECIES_FILE"
)

if [[ "${#TARGET_ROWS[@]}" -ne "$EXPECTED_TARGETS" ]]; then
    echo "Error: expected ${EXPECTED_TARGETS} rescue targets but selected ${#TARGET_ROWS[@]}." >&2
    echo "Refusing to run with an unexpected target set." >&2
    exit 1
fi

echo "Rescue worker ${WORKER_ID}/${N_WORKERS}: ${PIPELINE_QUERY_NAME}"
echo "Selected targets: ${#TARGET_ROWS[@]}; assignment: indexes congruent to ${WORKER_ID} modulo ${N_WORKERS}"

N_ASSIGNED=0
N_COMPLETED=0
N_SKIPPED=0
N_FAILED=0

for INDEX in "${!TARGET_ROWS[@]}"; do
    if (( INDEX % N_WORKERS != WORKER_ID )); then
        continue
    fi

    ((N_ASSIGNED += 1))
    ROW="${TARGET_ROWS[$INDEX]}"
    TARGET_CLADE="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $1}')"
    TAXID="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $7}')"
    SPECIES="$(printf '%s\n' "$ROW" | awk -F '\t' '{print $8}')"
    SPECIES_LABEL="$(printf '%s\n' "$SPECIES" | awk '{ if (NF >= 2) print $1, $2; else print $0 }')"

    SPECIES_SAFE="${TARGET_CLADE}_$(printf '%s\n' "$SPECIES_LABEL" | tr ' ' '_')"
    SPECIES_DIR="${OUTPUT_DIR}/${SPECIES_SAFE}"
    DONE_FILE="${DONE_DIR}/${SPECIES_SAFE}.done"
    RESULT_FILE="${SPECIES_DIR}/${TARGET_CLADE}_${PIPELINE_QUERY_NAME}.txt"

    if [[ -f "$DONE_FILE" ]]; then
        echo "Worker ${WORKER_ID} | ${SPECIES_SAFE} already completed; skipping."
        ((N_SKIPPED += 1))
        continue
    fi

    mkdir -p "$SPECIES_DIR"
    echo "Worker ${WORKER_ID} | taxid: ${TAXID} | species: ${SPECIES_SAFE}"
    echo "Running remote BLASTP: ${PIPELINE_QUERY_NAME} x ${SPECIES_SAFE}"

    if blastp \
        -query "$QUERY" \
        -db refseq_protein \
        -remote \
        -entrez_query "txid${TAXID}[Organism]" \
        -out "$RESULT_FILE" \
        -outfmt "6 qseqid saccver pident length qcovs evalue bitscore stitle" \
        -max_target_seqs 10
    then
        touch "$DONE_FILE"
        ((N_COMPLETED += 1))
        echo "Completed: ${PIPELINE_QUERY_NAME} x ${SPECIES_SAFE}"
    else
        BLAST_STATUS=$?
        ((N_FAILED += 1))
        echo "Warning: BLASTP failed with status ${BLAST_STATUS}: ${SPECIES_SAFE}" >&2
        echo "Continuing with this worker's remaining targets." >&2
    fi
done

echo "Worker ${WORKER_ID} summary: assigned=${N_ASSIGNED}, completed=${N_COMPLETED}, skipped=${N_SKIPPED}, failed=${N_FAILED}"

if [[ "$N_ASSIGNED" -ne 4 ]]; then
    echo "Error: expected 4 assigned targets for worker ${WORKER_ID}, got ${N_ASSIGNED}." >&2
    exit 1
fi
if [[ "$N_FAILED" -gt 0 ]]; then
    exit 1
fi
