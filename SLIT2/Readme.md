# SLIT project 

This project runs a protein-family homolog search and tree-building workflow for SLIT. 
The results and analysis for this project is deeply descripted in: 
https://liha0-0.notion.site/2026-06-04-SLIT2-Mir-218-project-3735d58f3b8c8074a58cf628024585fa?pvs=74 




# Homolog Search and Phylogenomics Pipeline

The current data are for the SLIT family, but the scripts are parameterized so
the same pipeline can be reused for another protein family by changing the query
FASTA and a few environment variables.

The main entry point is:

```bash
Scripts/00-Run-SLIT2-pipeline.sh
```

The script name is historical. The pipeline itself is no longer tied to the
absolute path of this SLIT2 directory.

## What It Does

The wrapper runs these stages in order:

```text
blast -> explore -> download -> filter -> count_homologs -> align -> iqtree -> compare_tree
```

Stage summary:

```text
blast
  Runs remote BLASTP against NCBI refseq_protein for each species in
  Scripts/species.txt. Results are written under Blast/Results, with one
  .done marker per species under Blast/Done.

explore
  Reads the BLAST result tables, removes hits whose subject species do not
  match the target species, writes Scripts/Blast-high-scoring-hits.txt, and
  creates diagnostic plots in Scripts/.

download
  Downloads protein FASTA records for the high-scoring hit accessions.

filter
  Uses NCBI protein-to-gene links to collapse isoforms by Gene ID, keeping the
  longest isoform per gene, then writes the homolog FASTA.

count_homologs
  Counts retained homologs per species and writes a count plot.

align
  Aligns the homolog FASTA with PRANK by default, or MAFFT if requested.

iqtree
  Builds a maximum-likelihood tree from the alignment with IQ-TREE.

compare_tree
  Roots the tree using the configured outgroup pattern, drops NoGeneID tips
  with a warning listing the dropped tips, shortens tip labels, and writes a
  tree PDF.
```

## Required Inputs

```text
Scripts/species.txt
  Tab-delimited species table. The pipeline uses Target_Clade, TaxID, and
  Organism_Name.

Blast/Query/<PIPELINE_QUERY_NAME>.fasta
  Protein query FASTA used for BLASTP. For the current SLIT run, the default is:
  Blast/Query/SLIT1-2-3.fasta
```

The query FASTA is also used by `02-Explore-blastp-result.R` to calculate the
mean query protein length. That mean length is only used as the dashed vertical
line in `Scripts/E-value_vs_Aligned_Length.png`; it does not affect BLAST,
filtering, alignment, or tree inference.

## Basic Usage

Run from the project root:

```bash
cd /rna/liha/phylogenomics_practice/SLIT2
```

Preview what would run:

```bash
Scripts/00-Run-SLIT2-pipeline.sh --dry-run
```

Run the full pipeline with resume-safe defaults:

```bash
Scripts/00-Run-SLIT2-pipeline.sh
```

Run without NCBI network-dependent stages, using existing cached outputs:

```bash
Scripts/00-Run-SLIT2-pipeline.sh --offline
```

Rerun only local analysis stages:

```bash
Scripts/00-Run-SLIT2-pipeline.sh --from explore --to compare_tree --force-local
```

Rerun only alignment, tree inference, and tree plotting:

```bash
Scripts/00-Run-SLIT2-pipeline.sh --from align --to compare_tree --force-local
```

Use MAFFT instead of PRANK:

```bash
Scripts/00-Run-SLIT2-pipeline.sh --aligner mafft
```

## Resume Behavior

The wrapper skips stages when the expected output already exists:

```text
BLAST
  Reuses existing .done markers in Blast/Done.

download/filter
  Reuses existing NCBI FASTA and filtered homolog FASTA unless --refresh-ncbi
  is given.

local stages
  Reuses plots, alignments, trees, and tree PDFs unless --force-local is given.
```

Use `--refresh-ncbi` when you want to redownload NCBI-dependent outputs. BLAST
still respects per-species `.done` markers.

## Configuration

Set these variables before running the wrapper when reusing the pipeline for
another family:

```bash
export PIPELINE_PROJECT_LABEL=ROBO
export PIPELINE_OUTPUT_PREFIX=ROBO
export PIPELINE_QUERY_NAME=ROBO1-2-3-4
export PIPELINE_OUTGROUP_PATTERN=Branchiostoma
```

Mean query length is calculated automatically from:

```text
Blast/Query/${PIPELINE_QUERY_NAME}.fasta
```

Useful optional variables:

```text
PIPELINE_INPUT_FASTA
  Override the FASTA used by the align step. This is useful after manually
  removing divergent sequences.

PIPELINE_REFERENCE_SPECIES
  Species highlighted with red points in the explore plot. Default: Homo sapiens.

PIPELINE_MAIL_TO and PIPELINE_MAIL_SCRIPT
  Optional notification settings used by BLAST/PRANK scripts if both are set.
```

## Curated Late-Stage Rerun

`Scripts/Run.sh` is a convenience Slurm script for rerunning the late stages
from a manually curated FASTA:

```bash
sbatch Scripts/Run.sh
```

By default it sets:

```bash
PIPELINE_INPUT_FASTA=Scripts/SLIT-homologs-removed-divergent.fasta
```

Then it runs:

```bash
Scripts/00-Run-SLIT2-pipeline.sh --from align --to compare_tree --aligner prank --force-local
```

For another family, either set `PIPELINE_INPUT_FASTA` yourself or create a file
named:

```text
Scripts/<PIPELINE_OUTPUT_PREFIX>-homologs-removed-divergent.fasta
```

## Main Outputs

```text
Scripts/Blast-high-scoring-hits.txt
Scripts/Blast-high-scoring-hits.fasta
Scripts/<PIPELINE_OUTPUT_PREFIX>-homologs.fasta
Scripts/<PIPELINE_OUTPUT_PREFIX>_Hit_Counts.png
Scripts/<PIPELINE_OUTPUT_PREFIX>_Homolog_Counts.png
Scripts/E-value_Distribution.png
Scripts/E-value_vs_Aligned_Length.png

Alignments/<PIPELINE_OUTPUT_PREFIX>_aligned.fasta
Alignments/<PIPELINE_OUTPUT_PREFIX>_aligned.aligner

Tree/<PIPELINE_OUTPUT_PREFIX>.treefile
Tree/Raw-alignment/<PIPELINE_OUTPUT_PREFIX>.treefile
Tree/Raw-alignment/Tree.pdf
```

## Environment Notes

The wrapper uses:

```bash
conda run -n r_env Rscript
```

for R scripts.

The Slurm helper scripts expect these conda environments/tools to exist:

```text
NCBI-download
  blastp and remote NCBI access

alignments
  mafft and prank

iqtree
  iqtree3

r_env
  R packages: data.table, ggplot2, dplyr, Biostrings, ape
```

The Python filter script also needs:

```text
pandas
requests
```
