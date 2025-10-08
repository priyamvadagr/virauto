#!/usr/bin/env python3
"""
======================================================================
Script: merge_netmhcpan_predictions.py
Author: Priyamvada Guha Roy
Description:
    This script merges NetMHCpan prediction output files with their 
    corresponding paired viral–human FASTA chunks to produce a combined 
    table of MHC binding predictions for non-identical viral–human peptide 
    pairs. The script is modular, allowing the user to specify the MHC 
    type (e.g., type1), class (e.g., A, B, or C), and peptide length 
    (k-mer size) via command-line arguments.

Enhancement:
    If the output file already exists, previously processed files are 
    detected via `source_file_human`, and only new NetMHCpan results 
    not already included in the output are processed and appended.

Workflow:
    1. Load all NetMHCpan .xls prediction files from the specified directory.
    2. Match each .xls file with its corresponding FASTA chunk (by chunk number).
    3. Read FASTA headers to extract peptide metadata.
    4. Merge FASTA metadata with NetMHCpan predictions (1:1 order).
    5. Split into viral and human subsets, merge on pair ID and allele.
    6. Remove duplicate and identical peptide pairs.
    7. Append only new results if output file exists.

CLI Arguments:
    --type   : MHC binding type (e.g., type1, type2)
    --class  : MHC class letter (e.g., A, B, C)
    --kmer   : Peptide length (e.g., 8, 9, 10)

Example Usage:
    python merge_netmhcpan_predictions.py --type type1 --class B --kmer 9
======================================================================
"""

# ====================================================
# IMPORTS & CLI ARGUMENTS
# ====================================================
import pandas as pd
import glob
import os
import re
import argparse

# ----------------------------------------------------
# Parse command-line arguments
# ----------------------------------------------------
parser = argparse.ArgumentParser(description="Merge NetMHCpan predictions with paired FASTA chunks.")
parser.add_argument("--type", required=True, help="MHC type (e.g., type1, type2)")
parser.add_argument("--class", required=True, dest="class_", help="MHC class letter (e.g., A, B, C)")
parser.add_argument("--kmer", required=True, type=int, help="Peptide length (e.g., 8, 9, 10)")
args = parser.parse_args()

type_ = args.type
class_ = args.class_
k_mer = args.kmer

# ====================================================
# USER-CONFIGURED PATHS (auto-generated)
# ====================================================
base_results = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan"
base_data = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan"

netmhc_dir = f"{base_results}/{k_mer}_mers/{type_}/{type_}_{class_}_chunks"
fasta_dir = f"{base_data}/paired_k_mers/{k_mer}_mers/chunks"
out_dir = f"{base_results}/{k_mer}_mers/{type_}/{type_}_processed"
out_file = os.path.join(out_dir, f"{type_}_{class_}_all_predictions_processed.txt")

# ====================================================
# Helper: extract chunk number from filename
# ====================================================
def extract_chunk_number(filename):
    """Extract chunk number from filename like 'HLA-C_04_157_matched_pairs_chunk_1049.xls'"""
    match = re.search(r'chunk_(\d+)', filename)
    return match.group(1) if match else None

# ====================================================
# Helper: read FASTA headers in order
# ====================================================
def read_fasta_headers(fasta_file):
    """Read FASTA headers into a DataFrame preserving order."""
    records = []
    with open(fasta_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                full_id = line[1:]
                parts = full_id.split("_")
                if len(parts) >= 4:
                    typ = parts[0]                    # VIRAL or HUMAN
                    pair_id = "_".join(parts[1:3])    # e.g. 118008_46
                    uniprot = parts[-1]               # e.g. Q9DSS3
                else:
                    typ, pair_id, uniprot = None, None, None
                records.append({
                    "full_id": full_id,
                    "type": typ,
                    "pair_id": pair_id,
                    "uniprot": uniprot
                })
    return pd.DataFrame(records)

# ====================================================
# Load previously processed files if output exists
# ====================================================
processed_files = set()
if os.path.exists(out_file):
    print(f"🔄 Existing output found: {out_file}")
    try:
        prev_df = pd.read_csv(out_file, sep="\t", usecols=["source_file_human"])
        processed_files = set(prev_df["source_file_human"].unique())
        print(f"Found {len(processed_files)} previously processed files.")
    except Exception as e:
        print(f"⚠️ Could not read existing output properly: {e}")
else:
    print("No previous output found — starting fresh.")

# ====================================================
# Process each .xls file
# ====================================================
xls_files = sorted(glob.glob(f"{netmhc_dir}/*.xls"))
print(f"\nFound {len(xls_files)} NetMHCpan .xls files total")

all_processed = []

for xls_file in xls_files:
    xls_basename = os.path.basename(xls_file)

    # Skip files already processed
    if xls_basename in processed_files:
        print(f"⏩ Skipping already processed file: {xls_basename}")
        continue

    print(f"\n{'='*60}")
    print(f"Processing: {xls_basename}")
    print(f"{'='*60}")
    
    # Extract chunk number
    chunk_num = extract_chunk_number(xls_file)
    if not chunk_num:
        print(f"Warning: Could not extract chunk number from {xls_basename}")
        continue
    
    fasta_path = os.path.join(fasta_dir, f"matched_pairs_chunk_{chunk_num}.fasta")
    if not os.path.exists(fasta_path):
        print(f"Warning: FASTA file not found: {fasta_path}")
        continue
    
    # Read FASTA metadata
    fasta_info = read_fasta_headers(fasta_path)
    print(f"  Loaded {len(fasta_info):,} FASTA headers")
    
    # Extract allele
    with open(xls_file) as fh:
        first_line = fh.readline().strip()
        allele = first_line.split()[0] if first_line else "Unknown"
    print(f"  Allele: {allele}")
    
    # Read NetMHCpan predictions
    df = pd.read_csv(xls_file, sep="\t", skiprows=1)
    print(f"  Loaded {len(df):,} predictions")
    
    if len(df) != len(fasta_info):
        print(f"  ⚠️ Row count mismatch — skipping {xls_basename}")
        continue
    
    # Attach metadata
    df["Allele"] = allele
    df["Peptide_ID_full"] = fasta_info["full_id"].values
    df["chunk"] = chunk_num
    df["source_file"] = xls_basename
    df["type"] = fasta_info["type"].values
    df["pair_id"] = fasta_info["pair_id"].values
    
    # Split by type
    viral = df[df["type"] == "VIRAL"].copy()
    human = df[df["type"] == "HUMAN"].copy()
    
    merged = pd.merge(
        viral,
        human,
        on=["pair_id", "Allele"],
        suffixes=("_viral", "_human"),
        how="inner"
    )
    print(f"  Merged {len(merged):,} viral–human pairs")
    
    # Filter duplicates and identicals
    merged_diff = merged.drop_duplicates()
    merged_diff = merged_diff[merged_diff["Peptide_viral"] != merged_diff["Peptide_human"]].copy()
    
    all_processed.append(merged_diff)

# ====================================================
# Combine and save
# ====================================================
if not all_processed:
    print("No new files to process.")
else:
    new_df = pd.concat(all_processed, ignore_index=True)
    print(f"✅ Total new non-identical pairs: {len(new_df):,}")

    os.makedirs(out_dir, exist_ok=True)
    if os.path.exists(out_file):
        print(f"Appending new data to existing file → {out_file}")
        with open(out_file, "a") as f:
            new_df.to_csv(f, sep="\t", index=False, header=False)
    else:
        print(f"Creating new output file → {out_file}")
        new_df.to_csv(out_file, sep="\t", index=False)

print("\n✅ Script complete.")
