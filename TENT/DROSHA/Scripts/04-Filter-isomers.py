from pathlib import Path
import os
import time
import xml.etree.ElementTree as ET

import pandas as pd
import requests

script_dir = Path(__file__).resolve().parent
output_prefix = os.environ.get("PIPELINE_OUTPUT_PREFIX", "SLIT")
metadata_file = Path(os.environ.get("PIPELINE_HIGH_HITS_TABLE", script_dir / "Blast-high-scoring-hits.txt"))
fasta_file = Path(os.environ.get("PIPELINE_HIGH_HITS_FASTA", script_dir / "Blast-high-scoring-hits.fasta"))
output_fasta = Path(os.environ.get("PIPELINE_FILTERED_FASTA", script_dir / f"{output_prefix}-homologs.fasta"))

print("0. Loading accession-to-organism metadata...")
try:
    meta_df = pd.read_csv(metadata_file, sep="\t")
    acc_to_organism = dict(zip(meta_df["saccver"], meta_df["organism"]))
    print(f"-> Loaded {len(acc_to_organism)} unique metadata mappings.\n")
except Exception as e:
    print(f"Failed to read the metadata file {metadata_file}: {e}")
    raise SystemExit(1)

print("1. Extracting sequences, lengths, and accessions from FASTA...")
if not fasta_file.is_file() or fasta_file.stat().st_size == 0:
    print(f"Missing or empty FASTA file: {fasta_file}")
    raise SystemExit(1)

sequences = {}
current_id = ""
with fasta_file.open() as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            current_id = line.split()[0][1:]
            sequences[current_id] = ""
        elif current_id:
            sequences[current_id] += line

data = []
for acc, seq in sequences.items():
    organism_name = acc_to_organism.get(acc, "Unknown_Clade_Species")
    data.append({
        "saccver": acc,
        "length": len(seq),
        "sequence": seq,
        "organism": organism_name,
    })

df = pd.DataFrame(data)
print(f"-> Found {len(df)} protein sequences.\n")
if df.empty:
    print(f"No sequences were parsed from {fasta_file}")
    raise SystemExit(1)

print("2. Looking up Gene IDs for each protein through NCBI; this may take 1-2 minutes...")
gene_ids = []
for index, acc in enumerate(df["saccver"]):
    base_acc = acc.split(".")[0]
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=gene&id={base_acc}"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        root = ET.fromstring(response.text)
        link = root.find(".//LinkSetDb/Link/Id")
        if link is not None:
            gene_id = link.text
        else:
            print(f"   -> Info: {acc} has no registered Gene ID; treating it as an independent sequence.")
            gene_id = f"NoGeneID_{acc}"
    except Exception:
        gene_id = f"Error_{acc}"

    gene_ids.append(gene_id)

    if (index + 1) % 10 == 0:
        print(f"   [{index + 1}/{len(df)}] mappings completed...")
    time.sleep(0.4)

df["GeneID"] = gene_ids

print("\n3. Removing isoforms and assigning paralog labels (H1, H2...)...")
filtered_df = df.sort_values("length", ascending=False).groupby("GeneID").head(1)
filtered_df = filtered_df.sort_values(by=["organism", "length"], ascending=[True, False])
filtered_df["H_num"] = filtered_df.groupby("organism").cumcount() + 1
filtered_df["H_label"] = "H" + filtered_df["H_num"].astype(str)

print(f"-> Filtered {len(df)} sequences down to {len(filtered_df)}.\n")

print("4. Writing the final FASTA file...")
with output_fasta.open("w") as f_out:
    for _, row in filtered_df.iterrows():
        header = f">{row['organism']}_{row['H_label']}_{row['GeneID']}_original_acc_{row['saccver']}"
        f_out.write(header + "\n")
        seq = row["sequence"]
        for i in range(0, len(seq), 80):
            f_out.write(seq[i:i + 80] + "\n")

print(f"Done. Output FASTA: {output_fasta}")
