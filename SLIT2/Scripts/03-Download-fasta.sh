ACC_LIST=$(awk 'NR>1 {print $2}' /rna/liha/phylogenomics_practice/SLIT2/Scripts/Blast-high-scoring-hits.txt | paste -sd, -)

curl -s -d "db=protein&id=${ACC_LIST}&rettype=fasta&retmode=text" \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
    > /rna/liha/phylogenomics_practice/SLIT2/Scripts/Blast-high-scoring-hits.fasta

