## This script compares the well-known species tree of mammals 
## with my NJ tree (geneious) and ML tree (IQ-TREE) made from COX2 sequences of the same species. 

# Libraries
library(tidyverse, quietly = TRUE) 
library(data.table, quietly = TRUE)  
library(ggplot2, quietly = TRUE)     
library(phangorn)   
library(ggtree)  





## Well known phylogeny of mammals from Alvarez-Carretero et al. 2022 
well_known_tree <- read.tree("/rna/liha/selection_project/Data/Species_Phylogeny/Final/mammalia_Alvarez-Carretero_2022_parsed.nwk") 
well_known_tree$tip.label <- gsub("_", " ", well_known_tree$tip.label)

# only retain species used for phylogenomics practice 
target_species <- readLines("/rna/liha/phylogenomics_practice/week3/species.txt")
target_species <- gsub("_", " ", target_species)
target_species <- trimws(target_species)

valid_species <- intersect(target_species, well_known_tree$tip.label)
missing_species <- setdiff(target_species, well_known_tree$tip.label)

well_known_tree <- keep.tip(well_known_tree, valid_species)

write.tree(well_known_tree, "/rna/liha/phylogenomics_practice/week4/well_known_tree.nwk")




## NJ Tree(geneious) made from COX2 sequences of the same species
cox2_nj_tree <- read.tree("/rna/liha/phylogenomics_practice/week3/cox2_primate_njtree.newick") 
cox2_nj_tree$tip.label <- gsub("_", " ", cox2_nj_tree$tip.label) 





## ML Tree (IQ-tree) made from COX2 sequences of the same species 
cox2_ml_tree <- read.tree("/rna/liha/phylogenomics_practice/week3/cox2_primate_aligned.fasta.treefile") 
cox2_ml_tree$tip.label <- gsub("_", " ", cox2_ml_tree$tip.label) 







## Plot the trees
par(mfrow = c(1, 3)) 

plot(well_known_tree, 
     main = "Known species tree of primates", 
     cex = 1.5, 
     edge.width = 2) 

plot(cox2_nj_tree, 
     main = "NJ tree of COX2 sequences", 
     cex = 1.5, 
     edge.width = 2)  

plot(cox2_ml_tree, 
     main = "ML tree of COX2 sequences", 
     cex = 1.5, 
     edge.width = 2 ) 

nodelabels(cox2_ml_tree$node.label, 
           cex = 1.5,  
           frame = "none", 
           adj = c(1.2, -0.2))

add.scale.bar(x = 1, y = 1.5,     
              length = 0.1,  
              cex = 1.5,  
              lwd = 2) 


par(mfrow = c(1, 1)) 