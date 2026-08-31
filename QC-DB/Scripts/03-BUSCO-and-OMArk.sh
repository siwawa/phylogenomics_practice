#!/bin/bash
# 03-BUSCO-and-OMArk.sh
#
# Runs the full proteome QC pipeline for one taxid/assembly:
#   1. Download protein.faa + genomic.gff + genomic.fna
#   2. agat: keep longest isoform per gene
#   3. gffread: translate filtered GFF into a protein fasta
#   4. Basic protein-set statistics
#   5. BUSCO (protein mode, eukaryota_odb10, offline)
#   6. OMArk (LUCA.h5)
#   7. Copy final results to a permanent location, then wipe the scratch dir
#
# Usage:
#   bash 03-BUSCO-and-OMArk.sh <TAXID> <ACCESSION> <TAXONOMY>
#
# Scratch work happens under /bank/ncbi/gnathostome/Tmp/<TAXID> and is
# deleted at the end of every run. Only final results are kept, under
# RESULTS_ROOT/<TAXONOMY>_<TAXID> (so results are readable and sortable by
# clade, not just by taxid — parse back out with e.g. `rev | cut -d_ -f1 | rev`).

set -euo pipefail

TAXID="${1:?Provide taxid as the first argument}"
ACCESSION="${2:?Provide assembly accession (e.g. GCF_...) as the second argument}"
TAXONOMY="${3:?Provide the Taxonomy group string as the third argument}"

TMPROOT=/bank/ncbi/gnathostome/Tmp
WORKDIR=${TMPROOT}/${TAXID}
RESULTS_ROOT=/rna/liha/phylogenomics_practice/QC-DB/Results
RESULTS_DIR=${RESULTS_ROOT}/${TAXONOMY}_${TAXID}
DONE_DIR=${RESULTS_ROOT}/.done
DONE_FLAG=${DONE_DIR}/${TAXONOMY}_${TAXID}.done

BUSCO_DB=/bank/ncbi/gnathostome/BUSCO
OMARK_DB=/bank/ncbi/gnathostome/OMArk/LUCA.h5

# agat and gffread both live in the busco env
ENV_NCBI=NCBI-download
ENV_AGAT=busco
ENV_BUSCO=busco
ENV_OMARK=omark

if [[ -f "${DONE_FLAG}" ]]; then
    echo "[INFO] taxid=${TAXID} already done (${DONE_FLAG}), skipping." >&2
    exit 0
fi

mkdir -p "${WORKDIR}"/{01_raw,02_filtered,03_stats,04_busco,05_omark,06_validation}
mkdir -p "${RESULTS_DIR}"
mkdir -p "${DONE_DIR}"

echo "[INFO] === taxid=${TAXID} accession=${ACCESSION} ===" >&2

# ---------------------------------------------------------------------------
# 1. Download protein.faa + genomic.gff + genomic.fna
# ---------------------------------------------------------------------------
echo "[INFO] Step 1: downloading assembly data..." >&2

source activate "${ENV_NCBI}"
cd "${WORKDIR}/01_raw"

datasets download genome accession "${ACCESSION}" --include protein,gff3,genome
unzip -q -o ncbi_dataset.zip
mv ncbi_dataset/data/"${ACCESSION}"/protein.faa .
mv ncbi_dataset/data/"${ACCESSION}"/genomic.gff .
mv ncbi_dataset/data/"${ACCESSION}"/*.fna genomic.fna
rm -rf ncbi_dataset ncbi_dataset.zip

n_raw_proteins=$(grep -c '^>' protein.faa || echo 0)
echo "[INFO] Raw protein.faa contains ${n_raw_proteins} sequences" >&2

if [[ "${n_raw_proteins}" -eq 0 ]]; then
    echo "[ERROR] protein.faa is empty for ${ACCESSION} — aborting" >&2
    rm -rf "${WORKDIR}"
    exit 1
fi

conda deactivate

# ---------------------------------------------------------------------------
# 2. agat: keep longest isoform per gene
# ---------------------------------------------------------------------------
echo "[INFO] Step 2: agat longest-isoform filtering..." >&2

source activate "${ENV_AGAT}"

agat_sp_keep_longest_isoform.pl \
    -f "${WORKDIR}/01_raw/genomic.gff" \
    -o "${WORKDIR}/02_filtered/longest_isoform.gff" \
    > "${WORKDIR}/02_filtered/agat.log" 2>&1

n_genes=$(awk -F'\t' '$3=="gene"' "${WORKDIR}/02_filtered/longest_isoform.gff" | wc -l)
n_mrna=$(awk -F'\t' '$3=="mRNA"' "${WORKDIR}/02_filtered/longest_isoform.gff" | wc -l)
echo "[INFO] After filtering: ${n_genes} genes, ${n_mrna} mRNA features (should be ~1:1)" >&2

# ---------------------------------------------------------------------------
# 3. gffread: translate filtered GFF into protein fasta (same env as agat)
# ---------------------------------------------------------------------------
echo "[INFO] Step 3: gffread translation..." >&2

if ! gffread "${WORKDIR}/02_filtered/longest_isoform.gff" \
    -g "${WORKDIR}/01_raw/genomic.fna" \
    -y "${WORKDIR}/02_filtered/longest_isoform.faa" \
    -E > "${WORKDIR}/02_filtered/gffread.log" 2>&1; then
    echo "[WARNING] gffread reported coordinate warnings or skipped invalid features (see gffread.log)" >&2
fi

n_filtered_proteins=$(grep -c '^>' "${WORKDIR}/02_filtered/longest_isoform.faa" || echo 0)
echo "[INFO] Filtered protein set contains ${n_filtered_proteins} sequences" >&2

conda deactivate

# ---------------------------------------------------------------------------
# 4. Basic protein-set statistics
# ---------------------------------------------------------------------------
echo "[INFO] Step 4: basic stats..." >&2

source activate "${ENV_BUSCO}"   # seqkit lives here

FAA="${WORKDIR}/02_filtered/longest_isoform.faa"
STATS="${WORKDIR}/03_stats/basic_stats.tsv"

# Deduplicate sequence header IDs if any exist
if seqkit fx2tab "$FAA" | awk '{print $1}' | sort | uniq -d | grep -q .; then
    echo "[WARNING] Duplicate sequence header IDs detected in FASTA — deduplicating using seqkit rmdup" >&2
    seqkit rmdup -n "$FAA" -o "${WORKDIR}/02_filtered/longest_isoform_dedup.faa"
    mv "${WORKDIR}/02_filtered/longest_isoform_dedup.faa" "$FAA"
fi

seqkit fx2tab -l -n "$FAA" | awk -F'\t' '{print $NF}' > "${WORKDIR}/03_stats/lengths.txt"

n_proteins=$(wc -l < "${WORKDIR}/03_stats/lengths.txt")
median_len=$(sort -n "${WORKDIR}/03_stats/lengths.txt" \
    | awk '{a[NR]=$1} END{print (NR%2==1)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2}')
n50=$(sort -nr "${WORKDIR}/03_stats/lengths.txt" \
    | awk '{a[NR]=$1; total+=$1} END{half=total/2; run=0; for(i=1;i<=NR;i++){run+=a[i]; if(run>=half){print a[i]; break}}}')
under50=$(awk '$1<50' "${WORKDIR}/03_stats/lengths.txt" | wc -l)
with_x=$(seqkit fx2tab -s "$FAA" | awk -F'\t' '$2 ~ /X/' | wc -l)
with_internal_stop=$(seqkit fx2tab -s "$FAA" \
    | awk -F'\t' '{seq=$2; sub(/\*$/,"",seq); if(seq ~ /\*/) print}' | wc -l)

echo -e "taxid\taccession\tn_proteins\tmedian_length\tprotein_length_N50\tunder_50aa\twith_X\twith_internal_stop" > "$STATS"
echo -e "${TAXID}\t${ACCESSION}\t${n_proteins}\t${median_len}\t${n50}\t${under50}\t${with_x}\t${with_internal_stop}" >> "$STATS"

# keep the full length distribution (for later plotting) but store it small —
# gzip a plain one-length-per-line file rather than keeping longest_isoform.faa
gzip -c "${WORKDIR}/03_stats/lengths.txt" > "${WORKDIR}/03_stats/lengths.txt.gz"
rm -f "${WORKDIR}/03_stats/lengths.txt"

conda deactivate

# ---------------------------------------------------------------------------
# 5. BUSCO
# ---------------------------------------------------------------------------
echo "[INFO] Step 5: BUSCO..." >&2

source activate "${ENV_BUSCO}"

busco -i "${WORKDIR}/02_filtered/longest_isoform.faa" \
    -l eukaryota_odb10 \
    -m protein \
    -o "busco_${TAXID}" \
    --out_path "${WORKDIR}/04_busco" \
    --download_path "${BUSCO_DB}" \
    --offline \
    -c 8 -f

conda deactivate

# ---------------------------------------------------------------------------
# 6. OMArk
# ---------------------------------------------------------------------------
echo "[INFO] Step 6: OMArk..." >&2

source activate "${ENV_OMARK}"

omamer search --db "${OMARK_DB}" \
    --query "${WORKDIR}/02_filtered/longest_isoform.faa" \
    --out "${WORKDIR}/05_omark/${TAXID}.omamer"

if ! omark -f "${WORKDIR}/05_omark/${TAXID}.omamer" \
    -d "${OMARK_DB}" \
    -o "${WORKDIR}/05_omark/"; then
    echo "[WARNING] OMArk taxonomy lineage lookup failed for taxid ${TAXID} (unmapped synthetic clade ID). Preserving omamer mappings." >&2
fi

conda deactivate

# ---------------------------------------------------------------------------
# 7. Copy final results to a permanent location, then wipe the scratch dir
# ---------------------------------------------------------------------------
echo "[INFO] Step 7: saving results to ${RESULTS_DIR} and cleaning up Tmp..." >&2

cp -r "${WORKDIR}/03_stats" "${RESULTS_DIR}/"
cp -r "${WORKDIR}/04_busco" "${RESULTS_DIR}/"
cp -r "${WORKDIR}/05_omark" "${RESULTS_DIR}/"
cp "${WORKDIR}/02_filtered/agat.log" "${RESULTS_DIR}/" 2>/dev/null || true

# scratch cleanup: everything under Tmp/<taxid> is disposable once results
# are safely copied out (raw downloads, filtered fasta/gff, genome fasta, etc.)
rm -rf "${WORKDIR}"

# Flag file: reaching this line means every step above ran without error
# (set -euo pipefail would have aborted otherwise). Each taxon can take a
# long time, so a flat directory of these — rather than checking inside
# every per-taxon Results/ subfolder — lets the full array, or a human,
# tell at a glance which taxa are already done. The early-exit guard near
# the top of this script reads this same flag to make reruns resumable.
{
    echo "taxid=${TAXID}"
    echo "accession=${ACCESSION}"
    echo "finished_at=$(date -Iseconds)"
} > "${DONE_FLAG}"

echo "[DONE] taxid=${TAXID} finished. Results under ${RESULTS_DIR}, Tmp cleaned." >&2
