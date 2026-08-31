import os
import glob
import pandas as pd
import urllib.request
import gzip

results_dir = '/rna/liha/phylogenomics_practice/QC-DB/Results'
summary_file = '/rna/liha/phylogenomics_practice/QC-DB/assembly_summary_refseq.txt'
output_file = '/rna/liha/phylogenomics_practice/QC-DB/Results/genome_quality_stats.tsv'

# 1. Get all unique accessions from the basic_stats files
accessions = set()
for tsv in glob.glob(os.path.join(results_dir, '*', '03_stats', 'basic_stats.tsv')):
    try:
        df = pd.read_csv(tsv, sep='\t')
        if 'accession' in df.columns:
            for acc in df['accession']:
                accessions.add(acc)
    except Exception as e:
        print(f"Error reading {tsv}: {e}")

print(f"Found {len(accessions)} accessions in Results.")

# 2. Get FTP paths from assembly_summary
ftp_paths = {}
with open(summary_file, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) > 19:
            acc = parts[0]
            ftp = parts[19]
            if acc in accessions:
                ftp_paths[acc] = ftp

print(f"Found {len(ftp_paths)} FTP paths in summary.")

# 3. Fetch stats
stats_data = []

for acc, ftp in ftp_paths.items():
    if ftp == "na": continue
    # ftp looks like: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
    base_name = ftp.rstrip('/').split('/')[-1]
    stats_url = f"{ftp}/{base_name}_assembly_stats.txt"
    
    try:
        req = urllib.request.Request(stats_url)
        with urllib.request.urlopen(req, timeout=10) as response:
            lines = response.read().decode('utf-8').splitlines()
            
            scaffold_n50 = None
            contig_n50 = None
            coverage = None
            
            for line in lines:
                if line.startswith('#'):
                    # Some files might have coverage in comments like:
                    # # Genome coverage: 100x
                    if "Genome coverage:" in line:
                        coverage = line.split(":")[-1].strip()
                elif line.startswith('all\tall\tall\tall\tscaffold-N50'):
                    scaffold_n50 = line.split('\t')[-1]
                elif line.startswith('all\tall\tall\tall\tcontig-N50'):
                    contig_n50 = line.split('\t')[-1]
            
            stats_data.append({
                'accession': acc,
                'scaffold_N50': scaffold_n50,
                'contig_N50': contig_n50,
                'genome_coverage': coverage
            })
            print(f"Fetched {acc}: scaffold-N50={scaffold_n50}, contig-N50={contig_n50}, cov={coverage}")
    except Exception as e:
        print(f"Error fetching {stats_url}: {e}")

df_out = pd.DataFrame(stats_data)
df_out.to_csv(output_file, sep='\t', index=False)
print(f"Saved stats to {output_file}")
