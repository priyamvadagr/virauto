#!/usr/bin/env python3
"""
======================================================================
Script: aggregate_allele_specific_results.py
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
    7. Batched disk writes (~10–20 merges at a time) to existing file

CLI Arguments:
    --type   : MHC binding type (e.g., type1, type2)
    --class  : MHC class letter (e.g., A, B, C)
    --kmer   : Peptide length (e.g., 8, 9, 10)

Example Usage:
    python merge_netmhcpan_predictions.py --type type1 --class B --kmer 9
======================================================================
"""

# ===================================================
# Import Required Libraries
# ===================================================
import pandas as pd
import os, re, glob, gzip, argparse
import concurrent.futures
from tqdm import tqdm


# ===================================================
# Parse Command-Line Arguments
# ===================================================
parser = argparse.ArgumentParser(
    description="Merge NetMHCpan predictions with paired viral–human FASTA metadata."
)
parser.add_argument(
    "--type", dest="type_", required=True,
    help="MHC type (e.g., type1 or type2)"
)
parser.add_argument(
    "--class", dest="class_", required=False, default="all",
    help="MHC class (e.g., A, B, C, all)"
)
parser.add_argument(
    "--chunk", dest="chunk", required=False, default=None,
    help="Specific chunk number to process"
)
parser.add_argument(
    "--kmer", dest="k_mer", type=int, required=True,
    help="Peptide length (e.g., 8, 9, 10)"
)
parser.add_argument(
    "--results_dir", default="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan",
    help="Base directory containing NetMHCpan results"
)
parser.add_argument(
    "--data_dir", default="/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan",
    help="Base directory containing FASTA chunks"
)
parser.add_argument(
    "--batch_size", type=int, default=15,
    help="Number of files to batch before writing to disk"
)
parser.add_argument(
    "--workers", type=int, default=8,
    help="Number of parallel worker processes"
)

args = parser.parse_args()
type_ = args.type_
class_ = args.class_
k_mer = args.k_mer
BATCH_SIZE = args.batch_size
MAX_WORKERS = args.workers
base_results = args.results_dir
base_data = args.data_dir
chunk = args.chunk


# ===================================================
# Configuration and Directory Setup
# ===================================================
if chunk:
    netmhc_dir = f"{base_results}/{k_mer}_mers/{type_}/{type_}_{class_}_chunks/{chunk}"
else:
    netmhc_dir = f"{base_results}/{k_mer}_mers/{type_}/{type_}_{class_}_chunks"
    
fasta_dir = f"{base_data}/paired_k_mers/{k_mer}_mers/chunks"
out_dir = f"{base_results}/{k_mer}_mers/{type_}/{type_}_processed"
os.makedirs(out_dir, exist_ok=True)
out_file = os.path.join(out_dir, f"{type_}_{class_}_{chunk}_all_predictions_processed.txt.gz")

print(f"⚙️  Configuration:")
print(f"   MHC Type:     {type_}")
print(f"   MHC Class:    {class_}")
print(f"   k-mer Length: {k_mer}")
print(f"   Output Dir:   {out_dir}")
print(f"   Batch Size:   {BATCH_SIZE}")
print(f"   Workers:      {MAX_WORKERS}")


# ===================================================
# Resume from Previously Processed Files
# ===================================================
processed_files = set()
if os.path.exists(out_file):
    prev = pd.read_csv(
        out_file,
        sep="\t",
        usecols=["source_file_human_human"],
        compression="gzip",
        low_memory=False
    )
    processed_files = set(prev["source_file_human_human"].unique())
    print(f"Found {len(processed_files)} processed files in {out_file}")


# ===================================================
# Define Regex Patterns
# ===================================================
chunk_re = re.compile(r'chunk_(\d+)')


# ===================================================
# Extract Chunk Number from Filename
# ===================================================
def extract_chunk_number(path):
    m = chunk_re.search(path)
    return m.group(1) if m else None


# ===================================================
# Read FASTA File (Extract type, pair_id, and uniprot)
# ===================================================
def read_fasta_fast(path):
    recs = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                full_id = line[1:].strip()  # e.g., "VIRAL_45607_22_P50822"
                parts = full_id.split('_')
                
                if len(parts) >= 4:
                    typ = parts[0]                      # "VIRAL" or "HUMAN"
                    pair_id = f"{parts[1]}_{parts[2]}"  # "45607_22"
                    uniprot = parts[3]                  # "P50822"
                    recs.append((typ, pair_id, uniprot))
                    
    return pd.DataFrame(recs, columns=["type", "pair_id", "uniprot"])


# ===================================================
# Annotate and Classify Viral-Human Peptide Pairs
# ===================================================
def process_one(xls_file):
    """Process a single NetMHCpan .xls file"""
    base = os.path.basename(xls_file)
    chunk = extract_chunk_number(base)
    fasta = os.path.join(fasta_dir, f"matched_pairs_chunk_{chunk}.fasta")
    if not os.path.exists(fasta):
        return None

    fasta_df = read_fasta_fast(fasta)
    with open(xls_file) as fh:
        first = fh.readline().strip()
        allele = first.split()[0] if first else "Unknown"

    df = pd.read_csv(xls_file, sep="\t", skiprows=1, engine="c", dtype=str)
    if len(df) != len(fasta_df):
        return None

    # Annotate with FASTA metadata
    df["Allele"] = allele
    df["chunk"] = chunk
    df["source_file_human"] = base
    df["type"] = fasta_df["type"].values
    df["pair_id"] = fasta_df["pair_id"].values
    df["uniprot"] = fasta_df["uniprot"].values

    # Separate viral and human peptides, then merge on pair_id + Allele
    viral = df[df["type"] == "VIRAL"]
    human = df[df["type"] == "HUMAN"]
    merged = pd.merge(
        viral,
        human,
        on=["pair_id", "Allele"],
        suffixes=("_viral", "_human")
    )
    
    # Remove identical peptide pairs
    merged = merged[merged["Peptide_viral"] != merged["Peptide_human"]]
    
    # Drop unnecessary columns - keep only one pair_id and the two uniprot columns
    # Note: pair_id has no suffix since it's a merge key
    if not merged.empty:
        # Drop the type columns (we don't need them after merge)
        cols_to_drop = [col for col in merged.columns if col in ['type_viral', 'type_human']]
        merged = merged.drop(columns=cols_to_drop)
    
    return merged if not merged.empty else None


# ===================================================
# Identify Unprocessed Files
# ===================================================
xls_files = sorted(glob.glob(f"{netmhc_dir}/*.xls"))
unprocessed = [f for f in xls_files if os.path.basename(f) not in processed_files]
print(f"Total XLS: {len(xls_files)} | To process: {len(unprocessed)}")

if not unprocessed:
    print("All files already processed.")
    exit(0)


# ===================================================
# Parallel Execution with Batching and Gzip Writing
# ===================================================
batch = []
total_written = 0

with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as exe:
    with tqdm(total=len(unprocessed), desc="Processing NetMHCpan chunks", dynamic_ncols=True) as pbar:
        for result in exe.map(process_one, unprocessed, chunksize=100):
            if result is not None:
                batch.append(result)

                # Write every BATCH_SIZE results
                if len(batch) >= BATCH_SIZE:
                    df_batch = pd.concat(batch, ignore_index=True)
                    header = not os.path.exists(out_file)
                    df_batch.to_csv(
                        out_file,
                        sep="\t",
                        index=False,
                        mode="a",
                        header=header,
                        compression="gzip"
                    )
                    total_written += len(df_batch)
                    batch.clear()  # Free memory after write
            pbar.update(1)


# ===================================================
# Write Remaining Batches (If Any)
# ===================================================
if batch:
    df_batch = pd.concat(batch, ignore_index=True)
    header = not os.path.exists(out_file)
    df_batch.to_csv(
        out_file,
        sep="\t",
        index=False,
        mode="a",
        header=header,
        compression="gzip"
    )
    total_written += len(df_batch)


# ===================================================
# Summary and Completion Message
# ===================================================
print(f"All chunks processed and appended ({total_written:,} rows total).")
print(f"Compressed output saved to: {out_file}")