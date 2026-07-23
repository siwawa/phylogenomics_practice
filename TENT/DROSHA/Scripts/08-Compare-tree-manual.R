suppressWarnings(suppressPackageStartupMessages({
  library(ape)
}))

tree_file <- "/rna/liha/phylogenomics_practice/TENT/TENT5/Tree/TENT5.treefile"
tree_pdf <- "/rna/liha/phylogenomics_practice/TENT/TENT5/Tree/Raw-alignment/Tree.pdf"
tree_title <- "TENT5 phylogeny"
outgroup_pattern <- "Cyclostomata"

raw_tree <- read.tree(tree_file)

outgroup_tips <- raw_tree$tip.label[grep(outgroup_pattern, raw_tree$tip.label, ignore.case = TRUE)]
if (length(outgroup_tips) > 0) {
  rooted_tree <- root(raw_tree, outgroup = outgroup_tips, resolve.root = TRUE, edgelabel = TRUE)
} else {
  warning("No outgroup tips matched outgroup_pattern=", outgroup_pattern, "; plotting unrooted tree.")
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

plot(rooted_tree, main = tree_title, cex = 0.6)

bs_values <- suppressWarnings(as.numeric(rooted_tree$node.label))
if (length(bs_values) > 0 && any(is.finite(bs_values))) {
  nodelabels(bs_values, cex = 0.5, frame = "none", col = "red", adj = c(1.2, -0.5))
}

add.scale.bar(x = 0.1, y = 5, length = 0.05, cex = 0.7, lwd = 2)
dev.off()