base_dir <- "/rna/liha/phylogenomics_practice/TENT/TENT2"
reciprocal_dir <- file.path(base_dir, "BLAST-reciprocal-hits")
result_dir <- file.path(base_dir, "BLASTp-results")
forward_hits <- read.delim(
  file.path(result_dir, "TENT2_BLASTp_results_isoform_resolved.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

hit_files <- list.files(
  reciprocal_dir,
  pattern = "_reciprocal_blast_hit\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

reciprocal_hits <- do.call(
  rbind,
  lapply(hit_files, read.delim, sep = "\t", check.names = FALSE)
)

# Exclude stale reciprocal files left by earlier forward-hit sets.
reciprocal_hits <- reciprocal_hits[
  reciprocal_hits$qseqid %in% forward_hits$saccver[!is.na(forward_hits$saccver)],
]

write.table(
  reciprocal_hits,
  file.path(result_dir, "TENT2_reciprocal_blast_result.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
