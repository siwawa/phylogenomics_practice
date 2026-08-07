base_dir <- "/rna/liha/phylogenomics_practice/TENT/TENT2"
species_info <- read.delim(
  file.path(base_dir, "Scripts", "species_info.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

blast_columns <- c(
  "qseqid", "saccver", "staxids", "sscinames", "pident",
  "length", "qcovs", "evalue", "bitscore", "stitle"
)

all_hits <- vector("list", nrow(species_info))

for (i in seq_len(nrow(species_info))) {
  species_id <- gsub(" ", "_", species_info$requested_species[i])
  hit_file <- file.path(
    base_dir, "BLAST-hits", species_id,
    paste0(species_id, "_blast_hit.txt")
  )

  if (file.size(hit_file) == 0) {
    hits <- as.data.frame(matrix(NA, nrow = 1, ncol = 10))
    names(hits) <- blast_columns
  } else {
    hits <- read.delim(
      hit_file,
      header = FALSE,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE,
      col.names = blast_columns
    )
  }

  # Species with no BLASTp hits receive one row with NA BLASTp values.
  all_hits[[i]] <- cbind(hits, species_info[rep(i, nrow(hits)), ])
}

result_dir <- file.path(base_dir, "BLASTp-results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  do.call(rbind, all_hits),
  file.path(result_dir, "TENT2_BLASTp_results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)
