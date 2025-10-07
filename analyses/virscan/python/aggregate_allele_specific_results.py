import pandas as pd
import glob
import os
import re

# ====================================================
# Input/output paths
# ====================================================
netmhc_dir = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_B_chunks"
fasta_dir = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks"
out_dir = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_concat"

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
    records = []
    with open(fasta_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                full_id = line[1:]  # remove ">"
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
# Process each .xls file
# ====================================================
xls_files = sorted(glob.glob(f"{netmhc_dir}/*.xls"))
print(f"Found {len(xls_files)} NetMHCpan .xls files")

all_processed = []

for xls_file in xls_files:
    xls_basename = os.path.basename(xls_file)
    print(f"\n{'='*60}")
    print(f"Processing: {xls_basename}")
    print(f"{'='*60}")
    
    # Extract chunk number from filename
    chunk_num = extract_chunk_number(xls_file)
    if not chunk_num:
        print(f"Warning: Could not extract chunk number from {xls_basename}")
        continue
    
    print(f"Chunk number: {chunk_num}")
    
    # Find corresponding FASTA file
    fasta_path = os.path.join(fasta_dir, f"matched_pairs_chunk_{chunk_num}.fasta")
    
    if not os.path.exists(fasta_path):
        print(f"Warning: FASTA file not found: {fasta_path}")
        continue
    
    # Read FASTA metadata (ordered)
    print(f"Reading FASTA: matched_pairs_chunk_{chunk_num}.fasta")
    fasta_info = read_fasta_headers(fasta_path)
    print(f"  Loaded {len(fasta_info):,} FASTA headers")
    
    # Extract allele from first line of .xls file
    with open(xls_file) as fh:
        first_line = fh.readline().strip()
        allele = first_line.split()[0] if first_line else "Unknown"
    
    print(f"  Allele: {allele}")
    
    # Read NetMHCpan predictions
    df = pd.read_csv(xls_file, sep="\t", skiprows=1)
    print(f"  Loaded {len(df):,} predictions")
    
    # Check if row counts match
    if len(df) != len(fasta_info):
        print(f"  WARNING: Row count mismatch!")
        print(f"    NetMHCpan predictions: {len(df):,}")
        print(f"    FASTA entries: {len(fasta_info):,}")
        continue
    
    # Attach FASTA metadata (1:1 ordered match) - only full peptide ID
    df["Allele"] = allele
    df["Peptide_ID_full"] = fasta_info["full_id"].values
    df["chunk"] = chunk_num
    df["source_file"] = xls_basename
    
    # Parse type and pair_id from full_id for merging purposes
    df["type"] = fasta_info["type"].values
    df["pair_id"] = fasta_info["pair_id"].values
    
    # Split by type
    viral = df[df["type"] == "VIRAL"].copy()
    human = df[df["type"] == "HUMAN"].copy()
    
    print(f"  Viral peptides: {len(viral):,}")
    print(f"  Human peptides: {len(human):,}")
    
    # Merge viral–human pairs
    merged = pd.merge(
        viral,
        human,
        on=["pair_id", "Allele"],
        suffixes=("_viral", "_human"),
        how="inner"
    )
    print(f"  Merged {len(merged):,} viral–human pairs")
    
    # Drop rows that are fully identical (all columns match)
    merged_diff = merged.drop_duplicates()
    print(f"  After removing duplicates: {len(merged_diff):,}")
    
    # Additionally filter out pairs with identical peptide sequences
    merged_diff = merged_diff[merged_diff["Peptide_viral"] != merged_diff["Peptide_human"]].copy()
    print(f"  Non-identical peptide pairs: {len(merged_diff):,}")
    
    # Drop the temporary type and pair_id columns from viral/human sides (keep the main ones)
    cols_to_drop = [col for col in merged_diff.columns if col in ['type_viral', 'type_human', 'pair_id']]
    merged_diff = merged_diff.drop(columns=cols_to_drop, errors='ignore')
    
    all_processed.append(merged_diff)

# ====================================================
# Combine all files and save final output
# ====================================================
print(f"\n{'='*60}")
print("Combining all files...")
print(f"{'='*60}")

if not all_processed:
    print("ERROR: No files were successfully processed!")
else:
    final_df = pd.concat(all_processed, ignore_index=True)
    print(f"Total non-identical pairs: {len(final_df):,}")
    
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "type1_B_all_predictions_processed.txt")
    final_df.to_csv(out_file, sep="\t", index=False)
    print(f"Saved processed table → {out_file}")
    
