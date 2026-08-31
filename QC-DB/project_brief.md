# Project brief: RefSeq proteome QC pipeline

## Goal
Assess the proteome (protein-coding gene set) quality of RefSeq taxa using
BUSCO, OMArk, and basic sequence statistics. Each taxon is represented by
exactly one selected RefSeq assembly (not the pooled set of every protein
under that taxid in the refseq_protein BLAST db).

## Species list
- Master table: `/rna/liha/phylogenomics_practice/QC-DB/Scripts/species_info.tsv`
- ~70 taxa spanning a very wide taxonomic range (Amoebozoa, Archaeplastida,
  Asgard Archaea, Cryptista, Discoba, Haptophyta, Obazoa/Holomycota/Holozoa
  including Metazoa subgroups, SAR), tracked via a `Taxonomy` column
  (e.g. `Amoebozoa-Evosea`).
- **`Asgard_Archaea` is excluded from analysis entirely** — not used
  downstream, skip it in every step.
- Columns include: Taxonomy, requested_species, selected_taxon, taxid,
  in_refseq_protein, total, NP_/XP_/WP_/YP_/OTHER counts, XP_ratio.

## Key design decision: single assembly per taxid
`refseq_protein` (local BLAST db) can pool proteins from multiple assemblies
under one taxid. Since AGAT-based isoform filtering requires a GFF3 tied to
one specific genome/protein FASTA, the workflow evaluates **one selected
RefSeq assembly's proteome per species**. The assembly is chosen as **the one
with the most pooled proteins** (determined via NCBI edirect IPG lookup:
count how many pooled-protein IPG hits map to each assembly, pick the max).

## Pipeline steps (per taxid + selected accession)
1. Download `protein.faa`, `genomic.gff`, `genomic.fna` via `datasets
   download genome accession <ACCESSION> --include protein,gff3,genome`
2. AGAT (`agat_sp_keep_longest_isoform.pl`) filters the GFF3 to one isoform
   (longest) per gene — avoids inflating BUSCO's "duplicated" count from
   alternative splicing.
3. `gffread -y` translates the filtered GFF3 (using genomic.fna) into a
   filtered protein FASTA — used as input for BUSCO/OMArk/stats. This
   intermediate FASTA is NOT kept in final results (large, regenerable).
   Noted caveat: gffread's mechanical translation may not perfectly match
   NCBI's original protein.faa for edge cases like selenoproteins.
4. Basic protein-set statistics: protein count, median length, protein-length
   N50 (explicitly a protein-length stat, not a genome-assembly N50), % under
   50 aa, % containing X, % containing internal stop codons. The full
   per-protein length list is kept but gzipped (`lengths.txt.gz`), not the
   large intermediate FASTA, to preserve the length distribution cheaply.
5. BUSCO (`-m protein`), lineage **unified to `eukaryota_odb10`** across all
   taxa (not taxon-specific lineages) — deliberate choice for cross-taxon
   comparability, matching the OMArk approach below. Run with `--offline`
   against a pre-downloaded local copy (already downloaded and verified, see
   Status below).
6. OMArk, database **unified to `LUCA.h5`** (not clade-specific DBs), same
   cross-taxon comparability reasoning. Isoforms already filtered upstream by
   AGAT, so OMArk's `--isoform_file` option is not needed.
7. Copy only final results (basic_stats.tsv + lengths.txt.gz, BUSCO
   short_summary + output dir, OMArk output, agat.log) to a permanent
   results directory, then **delete the entire per-taxid scratch directory**.

## Infrastructure
- HPC with Slurm; base data dir `/bank/ncbi/gnathostome`.
- `/bank/ncbi/gnathostome/BUSCO/lineages/eukaryota_odb10/` — BUSCO lineage
  dataset, **downloaded and confirmed complete** (hmms/, dataset.cfg,
  refseq_db.faa.gz all present).
- `/bank/ncbi/gnathostome/OMArk/LUCA.h5` — OMArk database (~8-9GB),
  downloaded via sbatch on the `pichu` node — confirm it finished before
  relying on it.
- `/bank/ncbi/gnathostome/refseq_protein/refseq_protein` — local
  refseq_protein BLAST db.
- `/bank/ncbi/gnathostome/Tmp/<taxid>/` — **scratch** working directory for
  the per-taxid pipeline run. Wiped (`rm -rf`) at the end of every run.
- `/rna/liha/phylogenomics_practice/QC-DB/Results/<Taxonomy>_<taxid>/` —
  **permanent** results directory, one per taxon. Named
  `<Taxonomy>_<taxid>` (not just taxid) for readability; parse back out with
  the trailing `_` since Taxonomy itself only contains hyphens.
- Conda envs: `busco` (BUSCO **and** agat/gffread — all three live here),
  `omark` (omamer/omark), `NCBI-download` (datasets, blastdbcmd, epost,
  efetch, seqkit).

## Scripts (in `/rna/liha/phylogenomics_practice/QC-DB/Scripts/`)
- `01-Check_pooled_counts.sh` (+ `.sbatch` wrapper) — adds
  n_assemblies_exact_taxid, n_proteins_pooled, n_assemblies_pooled columns
  to species_info.tsv → `species_info_pooled_counts.tsv`. **Already run once**
  — see Known issue below for gaps needing a retry pass.
- `02-Select-assemblies.sh species_info.tsv [output_tsv]` — for each taxid
  (excluding Asgard_Archaea), picks the assembly accession with the most
  pooled proteins. Writes `selected_assemblies.tsv` (default output path,
  next to the script) with columns: taxid, taxonomy, selected_accession,
  n_proteins_in_selected_assembly. Incremental: reuses any row already in
  the output file by taxid instead of redoing its IPG lookup, so rerunning
  after adding species to `species_info.tsv` only queries the new taxids
  (and retries any taxid that failed last time, since a failed lookup never
  gets a row written). Writes atomically (temp file + `mv`).
- `03-BUSCO-and-OMArk.sh <TAXID> <ACCESSION> <TAXONOMY>` — full per-taxid
  pipeline (download → agat → gffread → stats → BUSCO → OMArk → copy
  results → wipe scratch dir). Skips immediately if
  `Results/.done/<TAXONOMY>_<TAXID>.done` already exists, so resubmitting
  the array over a `selected_assemblies.tsv` that includes already-done
  taxa is safe and cheap.
- `04-Submit-array.sbatch` — Slurm array wrapper. Reads
  `selected_assemblies.tsv` line-by-line by `SLURM_ARRAY_TASK_ID`, calls
  `03-BUSCO-and-OMArk.sh` with taxid/accession/taxonomy. **Must be submitted
  with `--array=1-N%5` where N = actual data-row count of
  selected_assemblies.tsv**:
  ```bash
  N=$(($(wc -l < selected_assemblies.tsv) - 1))
  sbatch --array=1-${N}%5 04-Submit-array.sbatch
  ```
- `merge_summary.py` — NOT yet written; needs to parse BUSCO short_summary,
  OMArk output, and basic_stats.tsv (per Results/<Taxonomy>_<taxid>/ folder)
  into one combined row per taxon.

## Known issue: intermittent NCBI eutils 502 errors
While running `01-Check_pooled_counts.sh`, some `epost | efetch -format ipg`
calls hit `HTTP/1.0 502 Bad Gateway` from `eutils.ncbi.nlm.nih.gov` (transient
server-side, not a request problem — edirect auto-retries once but can still
fail twice in a row). The script does NOT crash on this — failed chunks are
silently skipped (`|| true`) and the affected taxid's `n_assemblies_pooled`
is left as `NA`, with a `[WARNING] IPG lookup failed for taxid $taxid` log
line. Two follow-ups needed:
1. **Immediate**: identify taxa with `n_assemblies_pooled == NA` (and
   `in_refseq_protein == YES`) in the completed
   `species_info_pooled_counts.tsv`, rebuild a subset input file for just
   those taxa, rerun `01-Check_pooled_counts.sh` on the subset, and merge
   the retry results back into the full table (replacing the NA rows).
2. **Structural**: add per-chunk retry logic (e.g. up to 3 attempts with a
   backoff sleep) to the `epost | efetch` calls in `01-Check_pooled_counts.sh`,
   `02-Select-assemblies.sh`, and (if it ends up needing IPG lookups)
   `03-BUSCO-and-OMArk.sh`, so this doesn't recur when those scripts run
   across all ~70 taxa. Consider logging permanently-failed taxids to a
   `needs_review.txt` file rather than silently leaving NA.

## Status / next steps
- BUSCO (`eukaryota_odb10`) and OMArk (`LUCA.h5`) databases: downloaded,
  BUSCO confirmed complete by directory listing; confirm OMArk's sbatch job
  finished too before relying on it.
- `01-Check_pooled_counts.sh` has been run once already but has NA gaps from
  the 502 issue above — resolve those first (see Known issue).
- Still need to: add chunk-retry logic to the scripts, rerun
  `02-Select-assemblies.sh` to produce `selected_assemblies.tsv`, run a
  smoke test of `03-BUSCO-and-OMArk.sh` on one small taxon end to end
  (verify env activation — agat/gffread both live in the `busco` env, output
  paths, gzip step, Tmp cleanup all work as expected), write
  `merge_summary.py`, then submit the full Slurm array (`--array=1-N%5`)
  across all taxa.
- Handoff note: prior work (through the pooled-counts run and initial
  drafts of scripts 01-04) was done interactively in Claude chat; this repo
  state is now being handed off to Claude Code for the remaining
  implementation and execution work.
- If running long/unattended tasks on the HPC: prefer working inside
  `tmux`/`screen` on the login node so work continues even if the local
  client disconnects; Slurm-submitted jobs (sbatch) already persist
  regardless of client connection.
