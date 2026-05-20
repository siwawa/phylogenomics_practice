#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate iqtree
iqtree3 -s /rna/liha/phylogenomics_practice/week3/cox2_primate_aligned.fasta -B 1000