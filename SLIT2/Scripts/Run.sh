#!/bin/bash
#SBATCH --job-name=SLIT
#SBATCH --nodes=1                 
#SBATCH --ntasks=1               
#SBATCH --cpus-per-task=4      
#SBATCH --mem=4G                 
#SBATCH --output=/rna/liha/phylogenomics_practice/SLIT2/logs/__%j.log    
#SBATCH --error=/rna/liha/phylogenomics_practice/SLIT2/logs/__%j.log      

cd /rna/liha/phylogenomics_practice/SLIT2/Scripts

bash 00-Run-SLIT2-pipeline.sh \
  --from blast \
  --to compare_tree \
  --aligner prank \