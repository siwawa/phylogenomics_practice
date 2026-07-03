suppressWarnings(suppressPackageStartupMessages({
  library(ape)
}))

script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (!is.na(value) && nzchar(value)) value else default
}

scripts_dir <- script_dir()
base_dir <- normalizePath(file.path(scripts_dir, ".."), mustWork = TRUE)
output_prefix <- get_env("PIPELINE_OUTPUT_PREFIX", "SLIT")
outgroup_pattern <- get_env("PIPELINE_OUTGROUP_PATTERN", "Branchiostoma")
raw_tree_dir <- file.path(base_dir, "Tree", "Raw-alignment")
tree_file <- get_env("PIPELINE_TREE_FILE", file.path(raw_tree_dir, paste0(output_prefix, ".treefile")))
tree_pdf <- get_env("PIPELINE_TREE_PDF", file.path(raw_tree_dir, "Tree.pdf"))
tree_title <- get_env("PIPELINE_TREE_TITLE", paste(output_prefix, "phylogeny"))

if (!file.exists(tree_file) || file.info(tree_file)$size == 0) {
  stop("Missing or empty tree file: ", tree_file)
}

dir.create(dirname(tree_pdf), recursive = TRUE, showWarnings = FALSE)
raw_tree <- read.tree(tree_file)

outgroup_tips <- raw_tree$tip.label[grep(outgroup_pattern, raw_tree$tip.label, ignore.case = TRUE)]
if (length(outgroup_tips) > 0) {
  rooted_tree <- root(raw_tree, outgroup = outgroup_tips, resolve.root = TRUE, edgelabel = TRUE)
} else {
  warning("No outgroup tips matched PIPELINE_OUTGROUP_PATTERN=", outgroup_pattern, "; plotting unrooted tree.")
  rooted_tree <- raw_tree
}

no_gene_id_tips <- rooted_tree$tip.label[grep("NoGeneID", rooted_tree$tip.label, ignore.case = TRUE)]
if (length(no_gene_id_tips) > 0) {
  warning(
    "Dropping tips without Gene ID: ",
    paste(no_gene_id_tips, collapse = ", "),
    call. = FALSE
  )
  rooted_tree <- drop.tip(rooted_tree, no_gene_id_tips)
}

rooted_tree$tip.label <- gsub("_(H[0-9]+)_.*$", "_\\1", rooted_tree$tip.label)

pdf(tree_pdf, width = 15, height = 15)
on.exit(dev.off(), add = TRUE)

plot(rooted_tree, main = tree_title, cex = 0.6)

bs_values <- suppressWarnings(as.numeric(rooted_tree$node.label))
if (length(bs_values) > 0 && any(is.finite(bs_values))) {
  nodelabels(bs_values, cex = 0.5, frame = "none", col = "red", adj = c(1.2, -0.5))
}

add.scale.bar(x = 0.1, y = 5, length = 0.05, cex = 0.7, lwd = 2)
