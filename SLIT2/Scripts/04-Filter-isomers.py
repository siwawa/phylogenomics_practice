import pandas as pd
import requests
import xml.etree.ElementTree as ET
import time

metadata_file = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/Blast-high-scoring-hits.txt"
fasta_file = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/Blast-high-scoring-hits.fasta"
output_fasta = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/SLIT-homologs.fasta"

print("0. Loading accession-to-organism metadata...")
# Read the metadata table from the tab-delimited file.
try:
    meta_df = pd.read_csv(metadata_file, sep='\t')
    # Map each accession (saccver) to its organism label.
    acc_to_organism = dict(zip(meta_df['saccver'], meta_df['organism']))
    print(f"-> Loaded {len(acc_to_organism)} unique metadata mappings.\n")
except Exception as e:
    print(f"Failed to read the metadata file: {e}")
    exit(1)


print("1. Extracting sequences, lengths, and accessions from FASTA...")
sequences = {}
current_id = ""

with open(fasta_file, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # Keep only the accession from headers such as '>XP_019635574.1 ...'.
            current_id = line.split()[0][1:]
            sequences[current_id] = ""
        else:
            sequences[current_id] += line


# Build a dataframe and add organism labels from the metadata mapping.
data = []
for acc, seq in sequences.items():
    # Use a fallback label when metadata is missing.
    organism_name = acc_to_organism.get(acc, "Unknown_Clade_Species")
    data.append({
        "saccver": acc, 
        "length": len(seq), 
        "sequence": seq,
        "organism": organism_name
    })

df = pd.DataFrame(data)
print(f"-> Found {len(df)} protein sequences.\n")


print("2. Looking up Gene IDs for each protein through NCBI; this may take 1-2 minutes...")
gene_ids = []
for index, acc in enumerate(df['saccver']):
    base_acc = acc.split('.')[0] 
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=gene&id={base_acc}"
    
    try:
        response = requests.get(url)
        root = ET.fromstring(response.text)
        
        link = root.find(".//LinkSetDb/Link/Id")
        if link is not None:
            gene_id = link.text
        else:
            print(f"   -> Info: {acc} has no registered Gene ID; treating it as an independent sequence.")
            gene_id = f"NoGeneID_{acc}" 
            
    except Exception as e:
        gene_id = f"Error_{acc}"
        
    gene_ids.append(gene_id)
    
    if (index + 1) % 10 == 0:
        print(f"   [{index + 1}/{len(df)}] mappings completed...")
    time.sleep(0.4) 

df['GeneID'] = gene_ids


print("\n3. Removing isoforms and assigning paralog labels (H1, H2...)...")
# 1) Group by GeneID and keep the longest sequence as the isoform representative.
filtered_df = df.sort_values('length', ascending=False).groupby('GeneID').head(1)

# 2) Group by organism and assign H labels in descending length order.
filtered_df = filtered_df.sort_values(by=['organism', 'length'], ascending=[True, False])
filtered_df['H_num'] = filtered_df.groupby('organism').cumcount() + 1
filtered_df['H_label'] = 'H' + filtered_df['H_num'].astype(str)

print(f"-> Filtered {len(df)} sequences down to {len(filtered_df)}.\n")


print("4. Writing the final FASTA file...")
with open(output_fasta, "w") as f_out:
    for _, row in filtered_df.iterrows():
        # New FASTA header format: >Target_Clade_Binomial_name_H1_geneID_original_acc_XP_...
        header = f">{row['organism']}_{row['H_label']}_{row['GeneID']}_original_acc_{row['saccver']}"
        f_out.write(header + "\n")
        
        seq = row['sequence']
        for i in range(0, len(seq), 80):
            f_out.write(seq[i:i+80] + "\n")

print(f"Done. Output FASTA: {output_fasta}")