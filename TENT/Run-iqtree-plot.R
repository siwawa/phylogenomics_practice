suppressWarnings(suppressPackageStartupMessages({
  library(ape)
}))

alignment <- "/rna/liha/phylogenomics_practice/TENT/Query/Structurally-aligned-human-TENT.fasta"
prefix <- "/rna/liha/phylogenomics_practice/TENT/Query/Structurally-aligned-human-TENT"
tree_pdf <- paste0(prefix, ".pdf")
tree_title <- "Structurally aligned human TENT"

system2("conda", c("run", "-n", "iqtree", "iqtree3",
                   "-s", alignment,
                   "-B", "5000",
                   "--prefix", prefix,
                   "-T", "AUTO",
                   "--threads-max", "64"))

tree <- read.tree(paste0(prefix, ".treefile"))
tree <- unroot(tree)

# Clean tip labels
tree$tip.label <- gsub("Homo_sapiens", "", tree$tip.label, ignore.case = TRUE)
tree$tip.label <- gsub("model_0", "", tree$tip.label, ignore.case = TRUE)
tree$tip.label <- gsub("fold", "", tree$tip.label, ignore.case = TRUE)

# Remove repeated underscores/spaces and trim
tree$tip.label <- gsub("_+", "_", tree$tip.label)
tree$tip.label <- gsub(" +", " ", tree$tip.label)
tree$tip.label <- gsub("^[_ ]+|[_ ]+$", "", tree$tip.label)

pdf(tree_pdf, width = 15, height = 15)
plot(tree, type = "unrooted", main = tree_title, cex = 1.5)

bs_values <- suppressWarnings(as.numeric(tree$node.label))
if (length(bs_values) > 0 && any(is.finite(bs_values))) {
  nodelabels(bs_values, cex = 1, frame = "none", col = "red", adj = c(1.2, -0.5))
}

add.scale.bar(x = 0.1, y = 5, length = 1, cex = 0.7, lwd = 2)
dev.off()
