psiblast \
    -query mef2_I1only.fasta \
    -db nr \
    -remote \
    -entrez_query "Protostomia[Organism] OR Echinodermata[Organism] OR Hemichordata[Organism] OR Cephalochordata[Organism] OR Urochordata[Organism] OR Cyclostomata[Organism] OR Chondrichthyes[Organism] OR Euteleostomi[Organism]" \
    -num_iterations 4 \
    -inclusion_ethresh 0.01 \
    -num_alignments 1000 \
    -out mef2_psiblast.txt \
    -outfmt 6

