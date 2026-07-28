#!/usr/bin/env bash
#SBATCH --job-name=TENT2_RBH
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --nodelist=guanine
#SBATCH --partition=intel
#SBATCH --output=logs/TENT2_reciprocal_%j.log
#SBATCH --error=logs/TENT2_reciprocal_%j.log

set -euo pipefail

BASE_DIR="/rna/liha/phylogenomics_practice/TENT/TENT2"
RESULT_DIR="${BASE_DIR}/BLASTp-results"
OUTPUT_DIR="${RESULT_DIR}/Reciprocal"
FORWARD_HITS="${RESULT_DIR}/TENT2_BLASTp_clades.tsv"
REFSEQ_DB="/bank/ncbi/gnathostome/refseq_protein/refseq_protein"
HUMAN_DB="/bank/ncbi/gnathostome/HumanDB/human_prot"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

CANDIDATE_ACCESSIONS="${OUTPUT_DIR}/candidate_accessions.txt"
CANDIDATE_RECORDS="${OUTPUT_DIR}/candidate_accession_hash_sequence.tsv"
CANDIDATE_MAP="${OUTPUT_DIR}/candidate_to_representative.tsv"
CANDIDATE_FASTA="${OUTPUT_DIR}/candidate_representative_sequences.fasta"
HUMAN_TENT2="${OUTPUT_DIR}/human_TENT2_accessions.tsv"
RECIPROCAL_ALL="${OUTPUT_DIR}/TENT2_reciprocal_all.tsv"
ACCEPTED_REPRESENTATIVES="${OUTPUT_DIR}/accepted_representatives.txt"
ACCEPTED_ACCESSIONS="${OUTPUT_DIR}/accepted_accessions.txt"
RECIPROCAL_HITS="${OUTPUT_DIR}/TENT2_reciprocal_hits.tsv"
RECIPROCAL_FASTA="${OUTPUT_DIR}/TENT2_reciprocal_hit_sequences.fasta"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate NCBI-download

mkdir -p "$OUTPUT_DIR"

# Collect unique forward-hit accessions.
awk -F '\t' 'NR > 1 {print $2}' "$FORWARD_HITS" |
    sort -u > "$CANDIDATE_ACCESSIONS"

# Retrieve each requested accession and its exact sequence hash.
blastdbcmd \
    -db "$REFSEQ_DB" \
    -entry_batch "$CANDIDATE_ACCESSIONS" \
    -target_only \
    -outfmt $'%a\t%h\t%s' \
    > "$CANDIDATE_RECORDS"

# Collapse identical amino-acid sequences and retain the accession mapping.
: > "$CANDIDATE_FASTA"
: > "$CANDIDATE_MAP"
awk -F '\t' -v OFS='\t' \
    -v fasta="$CANDIDATE_FASTA" -v map="$CANDIDATE_MAP" '
    BEGIN {
        print "candidate_accession", "representative_accession" > map
    }
    {
        sequence_key = $2 SUBSEP $3
        if (!(sequence_key in representative)) {
            representative[sequence_key] = $1
            print ">" $1 > fasta
            print $3 > fasta
        }
        print $1, representative[sequence_key] > map
    }
' "$CANDIDATE_RECORDS"

# Identify every human TENT2/GLD2 accession in the human database.
blastdbcmd \
    -db "$HUMAN_DB" \
    -entry all \
    -outfmt $'%a\t%t' |
    awk -F '\t' -v OFS='\t' '
        BEGIN { IGNORECASE = 1 }
        $2 ~ /^poly\(A\) RNA polymerase GLD2([[:space:]]|$)/ ||
        $2 ~ /^terminal nucleotidyltransferase 2([[:space:]]|$)/ {
            print $1, $2
        }
    ' > "$HUMAN_TENT2"

# Retain ten human targets so tied best hits can be evaluated.
blastp \
    -query "$CANDIDATE_FASTA" \
    -db "$HUMAN_DB" \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -num_threads "$THREADS" \
    -outfmt "6 qseqid saccver evalue bitscore pident length qcovs stitle" \
    -out "$RECIPROCAL_ALL"

# Accept a representative when any maximum-bitscore hit is a human TENT2.
awk -F '\t' '
    FILENAME == ARGV[1] {
        tent2[$1] = 1
        next
    }
    {
        query = $1
        score = $4 + 0
        if (!(query in best_score) || score > best_score[query]) {
            best_score[query] = score
            accepted[query] = ($2 in tent2)
        } else if (score == best_score[query] && ($2 in tent2)) {
            accepted[query] = 1
        }
    }
    END {
        for (query in accepted) {
            if (accepted[query]) {
                print query
            }
        }
    }
' "$HUMAN_TENT2" "$RECIPROCAL_ALL" |
    sort > "$ACCEPTED_REPRESENTATIVES"

# Map accepted representatives back to every original accession.
awk -F '\t' '
    FILENAME == ARGV[1] {
        accepted[$1] = 1
        next
    }
    FNR > 1 && ($2 in accepted) {
        print $1
    }
' "$ACCEPTED_REPRESENTATIVES" "$CANDIDATE_MAP" |
    sort -u > "$ACCEPTED_ACCESSIONS"

# Preserve every organism, taxonomy, and clade row for accepted accessions.
awk -F '\t' '
    FILENAME == ARGV[1] {
        accepted[$1] = 1
        next
    }
    FNR == 1 || ($2 in accepted)
' "$ACCEPTED_ACCESSIONS" "$FORWARD_HITS" > "$RECIPROCAL_HITS"

# Save one sequence for each accepted exact-sequence representative.
awk '
    FILENAME == ARGV[1] {
        accepted[$1] = 1
        next
    }
    /^>/ {
        accession = substr($0, 2)
        emit = (accession in accepted)
    }
    emit
' "$ACCEPTED_REPRESENTATIVES" "$CANDIDATE_FASTA" > "$RECIPROCAL_FASTA"

echo "Candidate accessions: $(wc -l < "$CANDIDATE_ACCESSIONS")"
echo "Unique candidate sequences: $(grep -c '^>' "$CANDIDATE_FASTA")"
echo "Accepted accessions: $(wc -l < "$ACCEPTED_ACCESSIONS")"
echo "Reciprocal-hit table: $RECIPROCAL_HITS"
echo "Reciprocal-hit FASTA: $RECIPROCAL_FASTA"
