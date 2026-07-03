library(ape)
pdf("/rna/liha/phylogenomics_practice/SLIT2/Tree/03-Alignment-removed-species-PRANK/Tree.pdf", width = 15, height = 15)


raw_tree <- read.tree("/rna/liha/phylogenomics_practice/SLIT2/Tree/03-Alignment-removed-species-PRANK/SLIT.treefile")




# Root with Branchiostoma as outgroup; edgelabel=TRUE preserves bootstrap values on correct nodes
outgroup_tips <- raw_tree$tip.label[grep("Branchiostoma", raw_tree$tip.label, ignore.case = TRUE)] 
rooted_tree <- root(raw_tree, outgroup = outgroup_tips, resolve.root = TRUE, edgelabel = TRUE) 


# Remove species with "NoGeneID" in their names 
rooted_tree <- drop.tip(rooted_tree, rooted_tree$tip.label[grep("NoGeneID", rooted_tree$tip.label, ignore.case = TRUE)]) 

# Shorten tip labels to only include species names (remove after _H1, _H2, etc.) 
rooted_tree$tip.label <- gsub("_(H[0-9]+)_.*$", "_\\1", rooted_tree$tip.label)




plot(rooted_tree, main = "SLIT phylogeny", cex = 0.6)

bs_values <- as.numeric(rooted_tree$node.label)
nodelabels(bs_values, cex = 0.5, frame = "none", col = "red", adj = c(1.2, -0.5))

add.scale.bar(x = 0.1, y = 5, length = 0.05, cex = 0.7, lwd = 2)


dev.off()
