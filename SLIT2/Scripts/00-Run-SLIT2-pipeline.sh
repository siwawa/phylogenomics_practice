#!/usr/bin/env bash
# Wrapper for the SLIT homolog search, filtering, alignment, and tree workflow.

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS_DIR="${BASE_DIR}/Scripts"

SPECIES_FILE="${SCRIPTS_DIR}/species.txt"
BLAST_SCRIPT="${SCRIPTS_DIR}/01-Blastp-against-species.sh"
EXPLORE_SCRIPT="${SCRIPTS_DIR}/02-Explore-blastp-result.R"
DOWNLOAD_SCRIPT="${SCRIPTS_DIR}/03-Download-fasta.sh"
FILTER_SCRIPT="${SCRIPTS_DIR}/04-Filter-isomers.py"
COUNT_HOMOLOGS_SCRIPT="${SCRIPTS_DIR}/05-Count-homologs-per-species.R"
ALIGN_MAFFT_SCRIPT="${SCRIPTS_DIR}/06-Align-MAFFT.sh"
ALIGN_PRANK_SCRIPT="${SCRIPTS_DIR}/06-Align-PRANK.sh"
ALIGN_SCRIPT="$ALIGN_PRANK_SCRIPT"
IQTREE_SCRIPT="${SCRIPTS_DIR}/07-Iq-tree.sh"
COMPARE_TREE_SCRIPT="${SCRIPTS_DIR}/08-Compare-tree.R"

BLAST_RESULTS_DIR="${BASE_DIR}/Blast/Results"
BLAST_DONE_DIR="${BASE_DIR}/Blast/Done"
HIGH_HITS_TABLE="${SCRIPTS_DIR}/Blast-high-scoring-hits.txt"
HIGH_HITS_FASTA="${SCRIPTS_DIR}/Blast-high-scoring-hits.fasta"
FILTERED_FASTA="${SCRIPTS_DIR}/SLIT-homologs.fasta"
HOMOLOG_COUNTS_PNG="${SCRIPTS_DIR}/SLIT_Homolog_Counts.png"
ALIGNMENT_FASTA="${BASE_DIR}/Alignments/SLIT_aligned.fasta"
ALIGNMENT_ALIGNER_FILE="${BASE_DIR}/Alignments/SLIT_aligned.aligner"
TREE_DIR="${BASE_DIR}/Tree"
TREE_PREFIX="${TREE_DIR}/SLIT"
TREE_FILE="${TREE_PREFIX}.treefile"
RAW_TREE_DIR="${TREE_DIR}/Raw-alignment"
RAW_TREE_FILE="${RAW_TREE_DIR}/SLIT.treefile"
TREE_PDF="${RAW_TREE_DIR}/Tree.pdf"

DRY_RUN=0
OFFLINE=0
FORCE_LOCAL=0
REFRESH_NCBI=0
ALIGNER="prank"
ALIGN_RAN=0
IQTREE_RAN=0
FROM_STEP="blast"
TO_STEP="compare_tree"

STEPS=(blast explore download filter count_homologs align iqtree compare_tree)

usage() {
    cat <<'USAGE'
Usage: Scripts/00-Run-SLIT2-pipeline.sh [options]

Runs the SLIT2 pipeline in order:
  blast -> explore -> download -> filter -> count_homologs -> align -> iqtree -> compare_tree

Default behavior is resume-safe:
  - Existing BLAST done markers/results are reused.
  - Existing downloaded NCBI FASTA is reused.
  - Existing filtered homolog FASTA is reused.
  - Existing homolog-count, alignment, tree, and comparison files are reused.
  - PRANK is used for alignment by default; use --aligner mafft to run MAFFT.

Options:
  --offline
      Never run NCBI-dependent stages. Fail if cached outputs are missing.

  --refresh-ncbi
      Allow rerunning NCBI-dependent download/filter stages even if their
      output files already exist. BLAST still respects per-species .done files.

  --force-local
      Rerun local stages even if outputs exist: explore plots/table, homolog
      counts, the selected aligner, IQ-TREE, and tree comparison. NCBI-heavy
      outputs are still reused unless --refresh-ncbi is also given.

  --aligner mafft|prank
      Select the aligner for the align step. Default: prank.

  --from STEP
      Start at STEP. Valid steps: blast, explore, download, filter,
      count_homologs, align, iqtree, compare_tree.

  --to STEP
      Stop after STEP.

  --dry-run
      Print what would run without executing commands.

  -h, --help
      Show this help.

Examples:
  Scripts/00-Run-SLIT2-pipeline.sh
  Scripts/00-Run-SLIT2-pipeline.sh --aligner mafft
  Scripts/00-Run-SLIT2-pipeline.sh --offline
  Scripts/00-Run-SLIT2-pipeline.sh --from align --to compare_tree --force-local
USAGE
}

log() {
    printf '[SLIT2] %s\n' "$*"
}

die() {
    printf '[SLIT2] ERROR: %s\n' "$*" >&2
    exit 1
}

is_step() {
    local candidate="$1"
    local step
    for step in "${STEPS[@]}"; do
        [[ "$candidate" == "$step" ]] && return 0
    done
    return 1
}

step_index() {
    local candidate="$1"
    local i
    for i in "${!STEPS[@]}"; do
        if [[ "${STEPS[$i]}" == "$candidate" ]]; then
            printf '%s\n' "$i"
            return 0
        fi
    done
    return 1
}

step_enabled() {
    local step="$1"
    local current from to
    current="$(step_index "$step")"
    from="$(step_index "$FROM_STEP")"
    to="$(step_index "$TO_STEP")"
    [[ "$current" -ge "$from" && "$current" -le "$to" ]]
}

run_cmd() {
    log "RUN: $*"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        "$@"
    fi
}

require_file() {
    local path="$1"
    local label="$2"
    [[ -s "$path" ]] || die "Missing ${label}: ${path}"
}

species_total() {
    awk 'NR > 1 && NF { n++ } END { print n + 0 }' "$SPECIES_FILE"
}

species_line_count() {
    awk 'END { print NR + 0 }' "$SPECIES_FILE"
}

done_total() {
    find "$BLAST_DONE_DIR" -type f -name '*.done' 2>/dev/null | wc -l | tr -d ' '
}

blast_result_count() {
    find "$BLAST_RESULTS_DIR" -type f -name '*SLIT1-2-3.txt' -size +0c 2>/dev/null | wc -l | tr -d ' '
}

expected_done_markers() {
    awk -v done_dir="$BLAST_DONE_DIR" '
        NR > 1 && NF {
            species = $9 " " $10
            gsub(/ /, "_", species)
            print done_dir "/" $1 "_" species ".done"
        }
    ' "$SPECIES_FILE"
}

missing_done_markers() {
    local done_file
    while IFS= read -r done_file; do
        [[ -f "$done_file" ]] || printf '%s\n' "$done_file"
    done < <(expected_done_markers)
}

blast_complete() {
    local total
    total="$(species_total)"
    [[ "$total" -gt 0 && -z "$(missing_done_markers)" ]]
}

ensure_dirs() {
    mkdir -p \
        "${BASE_DIR}/Blast/logs" \
        "$BLAST_RESULTS_DIR" \
        "$BLAST_DONE_DIR" \
        "${BASE_DIR}/logs" \
        "${BASE_DIR}/Alignments" \
        "$TREE_DIR" \
        "$RAW_TREE_DIR"
}

run_blast() {
    local total done results missing n_lines
    total="$(species_total)"
    done="$(done_total)"
    results="$(blast_result_count)"
    missing="$(missing_done_markers | wc -l | tr -d ' ')"

    if blast_complete; then
        log "SKIP blast: found expected done markers for ${total} species (${done} total .done files, ${results} result tables)."
        return
    fi

    if [[ "$OFFLINE" -eq 1 ]]; then
        die "BLAST cache is incomplete (${missing}/${total} expected done markers missing, ${done} total .done files, ${results} result tables) and --offline was requested."
    fi

    n_lines="$(species_line_count)"
    [[ "$n_lines" -ge 2 ]] || die "No species rows found in ${SPECIES_FILE}."

    log "Submitting BLAST SLURM array 2-${n_lines}%3. Existing .done files will be skipped by ${BLAST_SCRIPT}."
    run_cmd sbatch --array="2-${n_lines}%3" --wait "$BLAST_SCRIPT"

    blast_complete || {
    log "Missing expected BLAST done markers after submission:"
    missing_done_markers | head -n 20 >&2
    die "BLAST did not complete cleanly. Check ${BASE_DIR}/Blast/logs."
    } 
}

run_explore() {
    if [[ -s "$HIGH_HITS_TABLE" && "$FORCE_LOCAL" -eq 0 ]]; then
        log "SKIP explore: ${HIGH_HITS_TABLE} already exists."
        return
    fi

    [[ "$(blast_result_count)" -gt 0 ]] || die "No BLAST result tables found under ${BLAST_RESULTS_DIR}."

    if [[ "$DRY_RUN" -eq 1 ]]; then
        run_cmd Rscript "$EXPLORE_SCRIPT"
        return
    fi

    local had_table status
    had_table=0
    [[ -s "$HIGH_HITS_TABLE" ]] && had_table=1

    log "RUN: Rscript ${EXPLORE_SCRIPT}"
    set +e
    Rscript "$EXPLORE_SCRIPT"
    status=$?
    set -e

    if [[ "$status" -ne 0 ]]; then
        if [[ "$had_table" -eq 0 && -s "$HIGH_HITS_TABLE" && ! -s "$FILTERED_FASTA" ]]; then
            log "NOTE: ${EXPLORE_SCRIPT} wrote ${HIGH_HITS_TABLE} but exited later because ${FILTERED_FASTA} is not available yet."
            return
        fi
        die "${EXPLORE_SCRIPT} failed with exit status ${status}."
    fi

    require_file "$HIGH_HITS_TABLE" "high-scoring BLAST table"
}

run_download() {
    require_file "$HIGH_HITS_TABLE" "high-scoring BLAST table"

    if [[ -s "$HIGH_HITS_FASTA" && "$REFRESH_NCBI" -eq 0 ]]; then
        log "SKIP download: ${HIGH_HITS_FASTA} already exists."
        return
    fi

    if [[ "$OFFLINE" -eq 1 ]]; then
        die "Downloaded FASTA is missing or refresh was requested, but --offline forbids NCBI EFetch."
    fi

    run_cmd bash "$DOWNLOAD_SCRIPT"
    require_file "$HIGH_HITS_FASTA" "downloaded high-scoring hit FASTA"
}

run_filter() {
    require_file "$HIGH_HITS_TABLE" "high-scoring BLAST table"
    require_file "$HIGH_HITS_FASTA" "downloaded high-scoring hit FASTA"

    if [[ -s "$FILTERED_FASTA" && "$REFRESH_NCBI" -eq 0 ]]; then
        log "SKIP filter: ${FILTERED_FASTA} already exists."
        return
    fi

    if [[ "$OFFLINE" -eq 1 ]]; then
        die "Filtered FASTA is missing or refresh was requested, but --offline forbids NCBI ELink."
    fi

    run_cmd python3 "$FILTER_SCRIPT"
    require_file "$FILTERED_FASTA" "filtered SLIT homolog FASTA"
}

run_count_homologs() {
    require_file "$FILTERED_FASTA" "filtered SLIT homolog FASTA"

    if [[ -s "$HOMOLOG_COUNTS_PNG" && "$FORCE_LOCAL" -eq 0 ]]; then
        log "SKIP count_homologs: ${HOMOLOG_COUNTS_PNG} already exists."
        return
    fi

    run_cmd Rscript "$COUNT_HOMOLOGS_SCRIPT"
    require_file "$HOMOLOG_COUNTS_PNG" "homolog-count plot"
}

select_aligner() {
    case "$1" in
        mafft)
            ALIGNER="mafft"
            ALIGN_SCRIPT="$ALIGN_MAFFT_SCRIPT"
            ;;
        prank)
            ALIGNER="prank"
            ALIGN_SCRIPT="$ALIGN_PRANK_SCRIPT"
            ;;
        *)
            die "Unknown aligner: $1. Expected 'mafft' or 'prank'."
            ;;
    esac
}

alignment_matches_requested_aligner() {
    [[ -s "$ALIGNMENT_FASTA" && -s "$ALIGNMENT_ALIGNER_FILE" ]] || return 1
    [[ "$(tr -d '[:space:]' < "$ALIGNMENT_ALIGNER_FILE")" == "$ALIGNER" ]]
}

run_align() {
    require_file "$FILTERED_FASTA" "filtered SLIT homolog FASTA"

    if [[ "$FORCE_LOCAL" -eq 0 ]] && alignment_matches_requested_aligner; then
        log "SKIP align: ${ALIGNMENT_FASTA} already exists for ${ALIGNER}."
        return
    fi

    if [[ -s "$ALIGNMENT_FASTA" && "$FORCE_LOCAL" -eq 0 ]]; then
        log "Existing alignment is absent an aligner marker or was built with another aligner; rerunning with ${ALIGNER}."
    fi

    run_cmd sbatch --wait "$ALIGN_SCRIPT"
    [[ "$DRY_RUN" -eq 1 ]] && return

    require_file "$ALIGNMENT_FASTA" "${ALIGNER} alignment"
    printf '%s\n' "$ALIGNER" > "$ALIGNMENT_ALIGNER_FILE"
    ALIGN_RAN=1
}

sync_tree_outputs() {
    if [[ -s "$TREE_FILE" ]]; then
        if [[ "$FORCE_LOCAL" -eq 1 || "$IQTREE_RAN" -eq 1 || ! -s "$RAW_TREE_FILE" ]]; then
            run_cmd cp -p "${TREE_PREFIX}".* "$RAW_TREE_DIR"/
        fi
    fi
}

run_iqtree() {
    require_file "$ALIGNMENT_FASTA" "SLIT alignment"

    if [[ "$ALIGN_RAN" -eq 0 && -s "$RAW_TREE_FILE" && "$FORCE_LOCAL" -eq 0 ]]; then
        log "SKIP iqtree: ${RAW_TREE_FILE} already exists."
        return
    fi

    if [[ "$ALIGN_RAN" -eq 0 && -s "$TREE_FILE" && "$FORCE_LOCAL" -eq 0 ]]; then
        log "SKIP iqtree: ${TREE_FILE} already exists; syncing files into ${RAW_TREE_DIR} for plotting."
        sync_tree_outputs
        require_file "$RAW_TREE_FILE" "tree file for plotting"
        return
    fi

    run_cmd sbatch --wait "$IQTREE_SCRIPT"
    [[ "$DRY_RUN" -eq 1 ]] && return

    IQTREE_RAN=1
    require_file "$TREE_FILE" "IQ-TREE treefile"
    sync_tree_outputs
    require_file "$RAW_TREE_FILE" "tree file for plotting"
}

run_compare_tree() {
    if [[ "$IQTREE_RAN" -eq 0 && -s "$TREE_PDF" && "$FORCE_LOCAL" -eq 0 ]]; then
        log "SKIP compare_tree: ${TREE_PDF} already exists."
        return
    fi

    if [[ ! -s "$RAW_TREE_FILE" && -s "$TREE_FILE" ]]; then
        sync_tree_outputs
    fi

    require_file "$RAW_TREE_FILE" "tree file for plotting"
    run_cmd Rscript "$COMPARE_TREE_SCRIPT"
    require_file "$TREE_PDF" "tree-comparison PDF"
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --offline)
            OFFLINE=1
            shift
            ;;
        --refresh-ncbi)
            REFRESH_NCBI=1
            shift
            ;;
        --force-local)
            FORCE_LOCAL=1
            shift
            ;;
        --aligner)
            [[ "$#" -ge 2 ]] || die "--aligner requires 'mafft' or 'prank'."
            select_aligner "$2"
            shift 2
            ;;
        --from)
            [[ "$#" -ge 2 ]] || die "--from requires a step name."
            FROM_STEP="$2"
            is_step "$FROM_STEP" || die "Unknown --from step: ${FROM_STEP}"
            shift 2
            ;;
        --to)
            [[ "$#" -ge 2 ]] || die "--to requires a step name."
            TO_STEP="$2"
            is_step "$TO_STEP" || die "Unknown --to step: ${TO_STEP}"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

if [[ "$(step_index "$FROM_STEP")" -gt "$(step_index "$TO_STEP")" ]]; then
    die "--from step must come before or equal --to step."
fi

cd "$BASE_DIR"
ensure_dirs

step_enabled blast && run_blast
step_enabled explore && run_explore
step_enabled download && run_download
step_enabled filter && run_filter
step_enabled count_homologs && run_count_homologs
step_enabled align && run_align
step_enabled iqtree && run_iqtree
step_enabled compare_tree && run_compare_tree

log "Pipeline wrapper finished."
