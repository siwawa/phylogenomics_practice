#!/usr/bin/env bash
#SBATCH --job-name=REFSEQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --nodelist=pichu
#SBATCH --output=logs/REFSEQ_%j.log
#SBATCH --error=logs/REFSEQ_%j.log


set -euo pipefail

SCRIPT_DIR=${CHECK_REFSEQ_WORK_DIR:-/rna/liha/phylogenomics_practice/TENT/TENT2/Scripts}
DB_DIR=/bank/ncbi/gnathostome/refseq_protein
DB="$DB_DIR/refseq_protein"
NAMES_DMP="$DB_DIR/taxonomy/names.dmp"
NODES_DMP="$DB_DIR/taxonomy/nodes.dmp"
BLASTDBCMD=/rna/liha/miniconda3/envs/NCBI-download/bin/blastdbcmd

SPECIES_LIST="$SCRIPT_DIR/species_list.tsv"
FINAL_SUMMARY="$SCRIPT_DIR/species_info.tsv"

for required_file in "$SPECIES_LIST" "$NAMES_DMP" "$NODES_DMP"; do
    if [[ ! -r "$required_file" ]]; then
        echo "ERROR: Required file is not readable: $required_file" >&2
        exit 1
    fi
done

if [[ ! -x "$BLASTDBCMD" ]]; then
    echo "ERROR: blastdbcmd is not executable: $BLASTDBCMD" >&2
    exit 1
fi

if [[ $(awk -F '\t' 'NR == 1 {print NF}' "$SPECIES_LIST") -lt 2 ]]; then
    echo "ERROR: Expected a tab-separated input with Taxonomy and Species_name columns: $SPECIES_LIST" >&2
    exit 1
fi

export BLASTDB="$DB_DIR"
export LC_ALL=C

TMP_DIR=$(mktemp -d "$SCRIPT_DIR/.curated_species.XXXXXX")
BLAST_PID=""

cleanup() {
    if [[ -n "$BLAST_PID" ]]; then
        kill "$BLAST_PID" 2>/dev/null || true
        wait "$BLAST_PID" 2>/dev/null || true
        BLAST_PID=""
    fi
    rm -rf -- "$TMP_DIR"
}

handle_signal() {
    cleanup
    trap - EXIT
    exit 130
}

trap cleanup EXIT
trap handle_signal INT TERM HUP

RESOLVED_ROOTS="$TMP_DIR/resolved_roots.tsv"
ROOT_TAXIDS="$TMP_DIR/root_taxids.txt"
RAW_RECORDS="$TMP_DIR/raw_taxid_accession.tsv"
TAXID_COUNTS="$TMP_DIR/taxid_prefix_counts.tsv"
SELECTED="$TMP_DIR/selected_taxids.tsv"
SELECTED_NAMES="$TMP_DIR/selected_taxid_names.tsv"
TMP_FINAL="$TMP_DIR/summary_with_header.tsv"

echo "[1/5] Resolving requested species against the local NCBI taxonomy..."

# Output columns: input_order, taxonomy_label, requested_species, root_taxid.
# Every NCBI name class is accepted so that current names and synonyms resolve.
awk '
BEGIN { OFS="\t" }

FILENAME == ARGV[1] {
    if (FNR == 1)
        next

    split($0, field, "\t")
    sub(/\r$/, "", field[2])

    if (field[2] != "") {
        input_order[++n] = n
        taxonomy[n] = field[1]
        requested[n] = field[2]
        wanted[field[2]] = 1
    }
    next
}

{
    split($0, field, "\t\\|\t")
    taxid = field[1]
    name = field[2]

    if (name in wanted) {
        key = name SUBSEP taxid
        if (!seen[key]++) {
            hits[name]++
            taxids[name] = taxids[name] ? taxids[name] "," taxid : taxid
        }
    }
}

END {
    # This legacy binomial is absent as an exact entry in the current
    # names.dmp. NCBI places the corresponding species under taxid 1890364.
    override["Protostelium aurantium"] = "1890364"

    for (i = 1; i <= n; i++) {
        name = requested[i]

        if (override[name] != "") {
            root = override[name]
            print "NOTICE: taxonomy override: " name " -> taxid " root > "/dev/stderr"
        } else if (hits[name] == 1) {
            root = taxids[name]
        } else {
            root = ""
            if (hits[name] == 0)
                print "WARNING: unresolved species: " name > "/dev/stderr"
            else
                print "WARNING: ambiguous species: " name " -> " taxids[name] > "/dev/stderr"
        }

        print i, taxonomy[i], name, root
    }
}
' "$SPECIES_LIST" "$NAMES_DMP" > "$RESOLVED_ROOTS"

awk -F '\t' '$4 ~ /^[0-9]+$/ {print $4}' "$RESOLVED_ROOTS" |
    sort -nu > "$ROOT_TAXIDS"

requested=$(wc -l < "$RESOLVED_ROOTS")
resolved=$(wc -l < "$ROOT_TAXIDS")
echo "      Resolved $resolved unique root taxids from $requested input rows."

if (( resolved == 0 )); then
    echo "ERROR: No input species could be resolved to an NCBI taxid." >&2
    exit 1
fi

echo "[2/5] Querying each species taxid and all of its descendants..."

# Taxid expansion is intentional. It includes strains, subspecies, cultivar
# groups, and other descendants. Each returned taxid is counted separately.
"$BLASTDBCMD" \
    -db "$DB" \
    -taxidlist "$ROOT_TAXIDS" \
    -target_only \
    -outfmt $'%T\t%a' \
    -out "$RAW_RECORDS" &
BLAST_PID=$!
wait "$BLAST_PID"
BLAST_PID=""

echo "[3/5] Counting protein records and accession prefixes per returned taxid..."

# Output columns: taxid, total, NP_, XP_, WP_, YP_, OTHER.
awk -F '\t' '
BEGIN { OFS="\t" }

NF >= 2 && $1 ~ /^[0-9]+$/ {
    taxid = $1
    prefix = "OTHER"

    if (match($2, /^[[:alpha:]]+_/))
        prefix = substr($2, RSTART, RLENGTH)

    total[taxid]++
    count[taxid SUBSEP prefix]++
}

END {
    for (taxid in total) {
        np = count[taxid SUBSEP "NP_"] + 0
        xp = count[taxid SUBSEP "XP_"] + 0
        wp = count[taxid SUBSEP "WP_"] + 0
        yp = count[taxid SUBSEP "YP_"] + 0
        other = total[taxid] - np - xp - wp - yp

        print taxid, total[taxid], np, xp, wp, yp, other
    }
}
' "$RAW_RECORDS" | sort -k1,1n > "$TAXID_COUNTS"

echo "[4/5] Selecting the most protein-rich descendant taxid per requested species..."

# Map every returned taxid upward through nodes.dmp. For each requested root,
# retain the descendant with the largest protein count. Ties prefer the root
# taxid itself, then the numerically smaller taxid.
# Output columns: root_taxid, selected_taxid, total, NP_, XP_, WP_, YP_, OTHER.
awk '
BEGIN { FS=OFS="\t" }

FILENAME == ARGV[1] {
    if ($4 ~ /^[0-9]+$/)
        root[$4] = 1
    next
}

FILENAME == ARGV[2] {
    split($0, field, "\t\\|\t")
    parent[field[1]] = field[2]
    next
}

{
    candidate = $1
    ancestor = candidate
    matched_root = ""

    while (ancestor != "" && ancestor != "1") {
        if (ancestor in root) {
            matched_root = ancestor
            break
        }
        ancestor = parent[ancestor]
    }

    if (matched_root == "")
        next

    candidate_total = $2 + 0
    replace = 0

    if (!(matched_root in best_taxid) || candidate_total > best_total[matched_root]) {
        replace = 1
    } else if (candidate_total == best_total[matched_root]) {
        if (candidate == matched_root && best_taxid[matched_root] != matched_root)
            replace = 1
        else if (candidate != matched_root && best_taxid[matched_root] != matched_root &&
                 candidate + 0 < best_taxid[matched_root] + 0)
            replace = 1
    }

    if (replace) {
        best_taxid[matched_root] = candidate
        best_total[matched_root] = candidate_total
        best_np[matched_root] = $3
        best_xp[matched_root] = $4
        best_wp[matched_root] = $5
        best_yp[matched_root] = $6
        best_other[matched_root] = $7
    }
}

END {
    for (root_taxid in best_taxid) {
        print root_taxid, best_taxid[root_taxid], best_total[root_taxid],
              best_np[root_taxid], best_xp[root_taxid], best_wp[root_taxid],
              best_yp[root_taxid], best_other[root_taxid]
    }
}
' "$RESOLVED_ROOTS" "$NODES_DMP" "$TAXID_COUNTS" > "$SELECTED"

# Resolve the current scientific name of every selected descendant taxid.
awk '
BEGIN { FS=OFS="\t" }

FILENAME == ARGV[1] {
    wanted[$2] = 1
    next
}

{
    split($0, field, "\t\\|\t")
    if ((field[1] in wanted) && field[4] ~ /^scientific name/)
        print field[1], field[2]
}
' "$SELECTED" "$NAMES_DMP" > "$SELECTED_NAMES"

echo "[5/5] Writing the final summary..."

printf 'Taxonomy\trequested_species\tselected_taxon\ttaxid\tin_refseq_protein\ttotal\tNP_\tXP_\tWP_\tYP_\tOTHER\tXP_ratio\n' \
    > "$TMP_FINAL"

awk '
BEGIN { FS=OFS="\t" }

FILENAME == ARGV[1] {
    root = $1
    selected_taxid[root] = $2
    total[root] = $3
    np[root] = $4
    xp[root] = $5
    wp[root] = $6
    yp[root] = $7
    other[root] = $8
    next
}

FILENAME == ARGV[2] {
    scientific_name[$1] = $2
    next
}

{
    taxonomy = $2
    requested = $3
    root = $4

    if (root == "" || root !~ /^[0-9]+$/) {
        print taxonomy, requested, "NA", "NA", "NO", 0, 0, 0, 0, 0, 0, "NA"
        next
    }

    if (!(root in selected_taxid)) {
        print taxonomy, requested, "NA", "NA", "NO", 0, 0, 0, 0, 0, 0, "NA"
        next
    }

    taxid = selected_taxid[root]
    selected_name = scientific_name[taxid]
    if (selected_name == "")
        selected_name = "(name unavailable)"

    ratio = (total[root] > 0) ? sprintf("%.4f", xp[root] / total[root]) : "NA"

    print taxonomy, requested, selected_name, taxid, "YES", total[root],
          np[root] + 0, xp[root] + 0, wp[root] + 0, yp[root] + 0,
          other[root] + 0, ratio
}
' "$SELECTED" "$SELECTED_NAMES" "$RESOLVED_ROOTS" >> "$TMP_FINAL"

# Replace the previous result only after every processing step succeeds.
mv -f -- "$TMP_FINAL" "$FINAL_SUMMARY"

echo "Done: $FINAL_SUMMARY"
echo ""
echo "Selection rule:"
echo "  For each requested species, select the species taxid or descendant taxid"
echo "  with the largest number of RefSeq protein records."
echo "  The taxid column and all protein statistics describe only that selected taxid."
echo "  All intermediate files are removed automatically."
