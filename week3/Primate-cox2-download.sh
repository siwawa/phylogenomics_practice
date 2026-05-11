#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate NCBI-download

primates=(
  "Carlito syrichta"
  "Cheirogaleus medius"
  "Chlorocebus aethiops"
  "Daubentonia madagascariensis"
  "Galago senegalensis"
  "Gorilla gorilla"
  "Homo sapiens"
  "Hylobates syndactylus"
  "Lagothrix lagotrica"
  "Lemur catta"
  "Macaca mulatta"
  "Nycticebus coucang"
  "Pan paniscus"
  "Papio anubis"
)

for species in "${primates[@]}"; do
  echo "Downloading COX2 for $species..."

  filename="${species// /_}_COX2.zip"
  datasets download gene symbol COX2 --taxon "$species" --filename "$filename"
  
  sleep 1
done

