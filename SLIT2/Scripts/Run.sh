#!/bin/bash
#SBATCH --job-name=SLIT
#SBATCH --nodes=1                 
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=4      
#SBATCH --mem=4G                 
#SBATCH --output=/rna/liha/phylogenomics_practice/SLIT2/logs/__%j.log    
#SBATCH --error=/rna/liha/phylogenomics_practice/SLIT2/logs/__%j.log      

cd /rna/liha/phylogenomics_practice/SLIT2/Scripts

export SLIT_INPUT_FASTA="/rna/liha/phylogenomics_practice/SLIT2/Scripts/SLIT-homologs-removed-divergent.fasta"

bash 00-Run-SLIT2-pipeline.sh \
  --from align \
  --to compare_tree \
  --aligner prank \
  --force-local