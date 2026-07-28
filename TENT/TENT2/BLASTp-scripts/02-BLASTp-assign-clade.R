result_dir <- "/rna/liha/phylogenomics_practice/TENT/TENT2/BLASTp-results"
taxonomy_dir <- "/bank/ncbi/gnathostome/refseq_protein/taxonomy"

blast_columns <- c(
  "qseqid", "saccver", "staxids", "sscinames", "pident",
  "length", "qcovs", "evalue", "bitscore", "stitle"
)

blast <- read.delim(
  file.path(result_dir, "TENT2.txt"),
  header = FALSE,
  col.names = blast_columns,
  quote = "",
  comment.char = ""
)

taxid_sets <- strsplit(blast$staxids, ";", fixed = TRUE)
species_sets <- strsplit(blast$sscinames, ";", fixed = TRUE)
row_index <- rep(seq_len(nrow(blast)), lengths(taxid_sets))

blast <- blast[row_index, ]
blast$staxids <- unlist(taxid_sets)
blast$sscinames <- unlist(species_sets)
rownames(blast) <- NULL

nodes_file <- file.path(taxonomy_dir, "nodes.dmp")
nodes_connection <- pipe(paste("cut -f1,3", shQuote(nodes_file)))
nodes <- read.delim(
  nodes_connection,
  header = FALSE,
  col.names = c("taxid", "parent_taxid")
)
parents <- integer(max(nodes$taxid))
parents[nodes$taxid] <- nodes$parent_taxid

names_file <- file.path(taxonomy_dir, "names.dmp")
names_connection <- pipe(paste(
  "awk -F '\\t' '$7 == \"scientific name\" {print $1 \"\\t\" $3}'",
  shQuote(names_file)
))
taxon_names <- read.delim(
  names_connection,
  header = FALSE,
  col.names = c("taxid", "name"),
  quote = "",
  comment.char = ""
)
scientific_names <- character(length(parents))
scientific_names[taxon_names$taxid] <- taxon_names$name

clades <- c(
  Is_Bacteria = 2,
  Is_Archaea = 2157,
  Is_Eukaryota = 2759,
  Is_Opistokonta = 33154,
  Is_Metazoa = 33208,
  Is_Deuterostomia = 33511,
  Is_Gnathostomata = 7776
)

taxids <- unique(blast$staxids)

get_lineage <- function(taxid) {
  taxid <- as.integer(taxid)
  lineage <- taxid
  while (taxid != 1 && parents[taxid] != 0) {
    taxid <- parents[taxid]
    lineage <- c(lineage, taxid)
  }
  lineage
}

lineages <- setNames(lapply(taxids, get_lineage), taxids)

blast$Taxonomy <- vapply(
  lineages[blast$staxids],
  function(lineage) paste(
    scientific_names[rev(lineage[!lineage %in% c(1, 131567)])],
    collapse = ";"
  ),
  character(1)
)

membership <- t(vapply(
  lineages,
  function(lineage) clades %in% lineage,
  logical(length(clades))
))
rownames(membership) <- taxids
colnames(membership) <- names(clades)

for (clade in names(clades)) {
  blast[[clade]] <- membership[blast$staxids, clade]
}

write.table(
  blast,
  file.path(result_dir, "TENT2_BLASTp_clades.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


