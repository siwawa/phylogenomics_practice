## This script compares the well-known species tree of mammals 
## with my NJ tree (geneious) and ML tree (IQ-TREE) made from COX2 sequences of the same species. 

# Libraries
library(tidyverse, quietly = TRUE) 
library(data.table, quietly = TRUE)  
library(ggplot2, quietly = TRUE)     
library(phangorn)   
library(ggtree)  





# UCE-1216 without filtering columns of low alignment quality 
uce_1216 <- read.tree("/rna/liha/phylogenomics_practice/UCE/Tree/uce-1216.treefile")
species_info <- fread("/rna/liha/phylogenomics_practice/UCE/harvest_UCE/00-species/species.tsv", header = TRUE, stringsAsFactors = FALSE) 

    # Organize species names for readability
    species_info$tip_format <- tolower(gsub(" ", "_", species_info$Organism_Name))

    clade_map <- setNames(species_info$Target_Clade, species_info$tip_format)
    new_labels <- sapply(uce_1216$tip.label, function(sp_name) {
    clade <- clade_map[sp_name]
    
    if (!is.na(clade)) {
        return(paste0(clade, "_", sp_name))
    } else {
        return(sp_name) 
    }
    })

    uce_1216$tip.label <- as.character(new_labels)

plot(uce_1216, cex = 0.6)
nodelabels(uce_1216$node.label, frame = "none", cex = 0.5, adj = c(1.2, -0.5)) 



# UCE-4587 with manually filtering columns of low alignment quality
# This results in about 200nt(...) of informative sites, which is likely not enough to resolve the phylogeny. 
uce_4587 <- read.tree("/rna/liha/phylogenomics_practice/UCE/Tree/uce-4587_aligned/uce-4587_aligned.treefile")
species_info <- fread("/rna/liha/phylogenomics_practice/UCE/harvest_UCE/00-species/species.tsv", header = TRUE, stringsAsFactors = FALSE) 

    # Organize species names for readability
    species_info$tip_format <- tolower(gsub(" ", "_", species_info$Organism_Name))

    clade_map <- setNames(species_info$Target_Clade, species_info$tip_format)
    new_labels <- sapply(uce_4587$tip.label, function(sp_name) {
    clade <- clade_map[sp_name]
    
    if (!is.na(clade)) {
        return(paste0(clade, "_", sp_name))
    } else {
        return(sp_name) 
    }
    })

    uce_4587$tip.label <- as.character(new_labels)

plot(uce_4587, cex = 0.6)
nodelabels(uce_4587$node.label, frame = "none", cex = 0.5, adj = c(1.2, -0.5)) 






