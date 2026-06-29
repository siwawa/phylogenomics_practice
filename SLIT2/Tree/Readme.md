Each tree is made up with different versions of input alignments:

01-Raw-alignment is from the original set of species, with baseline MAFFT. 
The alignment was not manually trimmed. 

02-Raw-alignment-added-species-PRANK added more species from 01, which includes 
C.milii(Chimera sharks), X.lavis(Frog with WGD), Teleosts(Including Salmoniformes, which underwent WGD) 
Phylogeny-aware alignments(PRANK) was used as an alternative alignment method, in hope to help with gappy regions. 

03-Alignment-removed-species-PRANK removed some species from 02, 
as some of them were highly divergent/gappy and may noise alignments or bias tree-building. 
The removed species were Gadiformes(Cods) and a highly divergent species in Urochordata. 
 
04-Alignment-trimmed-columns-PRANK removed some columns from 03, 
gappy columns were manually trimmed to remove noise from the alignments. 
The aligned length was reduced from 2600 to 1500. 