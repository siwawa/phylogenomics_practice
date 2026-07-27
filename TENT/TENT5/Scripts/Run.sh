#!/bin/bash
#SBATCH --job-name=TENT5_tree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --chdir=/rna/liha/phylogenomics_practice/TENT/TENT5
#SBATCH --output=/rna/liha/phylogenomics_practice/TENT/TENT5/logs/run_profile_tree_%j.log
#SBATCH --error=/rna/liha/phylogenomics_practice/TENT/TENT5/logs/run_profile_tree_%j.log
set -euo pipefail

BASE_DIR="/rna/liha/phylogenomics_practice/TENT/TENT5"
SCRIPT_DIR="${BASE_DIR}/Scripts"
ALIGNMENT_FASTA="${BASE_DIR}/Alignments/profile-alignment.fasta"
OUTPUT_PREFIX="TENT5_profile_QMAMMAL_IG4"
IQTREE_MODEL="Q.MAMMAL+I+G4"
TREE_PDF="${BASE_DIR}/Tree/Raw-alignment/${OUTPUT_PREFIX}.pdf"

cd "$BASE_DIR"

if [[ ! -s "$ALIGNMENT_FASTA" ]]; then
    echo "Error: profile alignment not found or empty: $ALIGNMENT_FASTA" >&2
    exit 1
fi

if ! awk '
    function check_sequence() {
        if (!seen_sequence) return
        if (sequence_length == 0) {
            print "Error: empty sequence for " sequence_name > "/dev/stderr"
            failed = 1
        } else if (alignment_length == 0) {
            alignment_length = sequence_length
        } else if (sequence_length != alignment_length) {
            print "Error: unequal sequence length for " sequence_name \
                ": " sequence_length " versus " alignment_length > "/dev/stderr"
            failed = 1
        }
        sequence_count++
    }
    /^>/ {
        check_sequence()
        sequence_name = substr($0, 2)
        sequence_length = 0
        seen_sequence = 1
        next
    }
    {
        gsub(/[[:space:]]/, "")
        sequence_length += length($0)
    }
    END {
        check_sequence()
        if (sequence_count < 2) {
            print "Error: alignment contains fewer than two sequences." > "/dev/stderr"
            failed = 1
        }
        if (!failed) {
            print "Validated alignment: " sequence_count \
                " sequences x " alignment_length " columns."
        }
        exit failed
    }
' "$ALIGNMENT_FASTA"; then
    echo "Error: invalid aligned FASTA: $ALIGNMENT_FASTA" >&2
    exit 1
fi

if ! awk '
    /^>/ {
        name = substr($0, 2)
        if (seen[name]++) {
            print "Error: duplicate FASTA header: " name > "/dev/stderr"
            duplicate = 1
        }
    }
    END {
        exit duplicate
    }
' "$ALIGNMENT_FASTA"; then
    exit 1
fi

if ! awk '
    /^>/ && $0 ~ /(Amphimedon|Oscarella|Sycon)/ {
        found = 1
    }
    END {
        exit(found ? 0 : 1)
    }
' "$ALIGNMENT_FASTA"; then
    echo "Error: no poriferan outgroup sequence found in the alignment." >&2
    exit 1
fi

export PIPELINE_BASE_DIR="$BASE_DIR"
export PIPELINE_PROJECT_LABEL="TENT5-profile-tree"
export PIPELINE_OUTPUT_PREFIX="$OUTPUT_PREFIX"
export PIPELINE_ALIGNMENT_FASTA="$ALIGNMENT_FASTA"
export PIPELINE_IQTREE_MODEL="$IQTREE_MODEL"
export PIPELINE_OUTGROUP_PATTERN="Amphimedon|Oscarella|Sycon"
export PIPELINE_TREE_TITLE="TENT5 profile-alignment phylogeny (Q.MAMMAL+I+G4)"
export PIPELINE_TREE_PDF="$TREE_PDF"

bash "${SCRIPT_DIR}/00-Wrapper.sh" \
    --from iqtree \
    --to compare_tree \
    --force-local
