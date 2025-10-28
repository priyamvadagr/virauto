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
    - Uses manifest file to track processed files (fast resume)
    - Splits output into multiple files to avoid large file issues
    - If the output file already exists, previously processed files are 
      detected via manifest, and only new NetMHCpan results 
      not already included in the output are processed and appended.

Workflow:
    1. Load all NetMHCpan .xls prediction files from the specified directory.
    2. Match each .xls file with its corresponding FASTA chunk (by chunk number).
    3. Read FASTA headers to extract peptide metadata.
    4. Merge FASTA metadata with NetMHCpan predictions (1:1 order).
    5. Split into viral and human subsets, merge on pair ID and allele.
    6. Remove duplicate and identical peptide pairs.
    7. Batched disk writes (~10–20 merges at a time) to split output files

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
import json
from datetime import datetime


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
parser.add_argument(
    "--max_rows_per_file", type=int, default=1000000,
    help="Maximum rows per output file before creating a new one"
)
parser.add_argument(
    "--force_reprocess", action="store_true",
    help="Ignore manifest and reprocess all files"
)

args = parser.parse_args()
type_ = args.type_
class_ = args.class_
k_mer = args.k_mer
BATCH_SIZE = args.batch_size
MAX_WORKERS = args.workers
MAX_ROWS_PER_FILE = args.max_rows_per_file
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

# Manifest file for tracking progress
manifest_file = os.path.join(out_dir, f"{type_}_{class_}_{chunk}_manifest.json")

print(f"⚙️  Configuration:")
print(f"   MHC Type:     {type_}")
print(f"   MHC Class:    {class_}")
print(f"   k-mer Length: {k_mer}")
print(f"   Output Dir:   {out_dir}")
print(f"   Batch Size:   {BATCH_SIZE}")
print(f"   Workers:      {MAX_WORKERS}")
print(f"   Max Rows/File: {MAX_ROWS_PER_FILE}")
print(f"   Manifest:     {manifest_file}")


# ===================================================
# Manifest Management Functions
# ===================================================
def load_manifest(manifest_file):
    """Load the manifest file that tracks processed files and output parts."""
    if os.path.exists(manifest_file) and not args.force_reprocess:
        with open(manifest_file, 'r') as f:
            manifest_data = json.load(f)
        return manifest_data
    else:
        return {
            "processed_files": {},  # {filename: {"part": part_num, "timestamp": timestamp}}
            "output_parts": {},      # {part_num: {"row_count": count, "files": [list of files]}}
            "last_part": 0,
            "total_rows": 0
        }


def save_manifest(manifest_file, manifest_data):
    """Save the manifest file."""
    with open(manifest_file, 'w') as f:
        json.dump(manifest_data, f, indent=2)


def get_next_output_file(out_dir, type_, class_, chunk, manifest_data):
    """Determine the next output file to write to based on manifest."""
    last_part = manifest_data["last_part"]
    
    # Check if we need a new part file
    if last_part == 0:
        # No files exist yet
        next_part = 1
        create_new = True
    else:
        # Check row count of last part
        last_part_info = manifest_data["output_parts"].get(str(last_part), {})
        row_count = last_part_info.get("row_count", 0)
        
        if row_count >= MAX_ROWS_PER_FILE:
            # Current file is full, create new one
            next_part = last_part + 1
            create_new = True
        else:
            # Continue with current file
            next_part = last_part
            create_new = False
    
    chunk_suffix = f"_{chunk}" if chunk else ""
    out_file = os.path.join(out_dir, f"{type_}_{class_}{chunk_suffix}_predictions_part{next_part:04d}.txt.gz")
    
    return out_file, next_part, create_new


# ===================================================
# Resume from Previously Processed Files
# ===================================================
# Load manifest to get previously processed files
manifest_data = load_manifest(manifest_file)
processed_files = set(manifest_data["processed_files"].keys())
print(f"📊 Manifest: {len(processed_files)} files already processed")


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
# Parallel Execution with Batching and Split File Writing
# ===================================================
batch = []
batch_files = []  # Track which files are in the current batch
total_written = 0
timestamp = datetime.now().isoformat()

with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as exe:
    with tqdm(total=len(unprocessed), desc="Processing NetMHCpan chunks", dynamic_ncols=True) as pbar:
        for xls_file, result in zip(unprocessed, exe.map(process_one, unprocessed, chunksize=100)):
            if result is not None:
                batch.append(result)
                batch_files.append(os.path.basename(xls_file))

                # Write every BATCH_SIZE results
                if len(batch) >= BATCH_SIZE:
                    df_batch = pd.concat(batch, ignore_index=True)
                    
                    # Get the appropriate output file (may create new part if needed)
                    out_file, part_num, create_new = get_next_output_file(
                        out_dir, type_, class_, chunk, manifest_data
                    )
                    
                    # Write to file
                    header = create_new or not os.path.exists(out_file)
                    df_batch.to_csv(
                        out_file,
                        sep="\t",
                        index=False,
                        mode="a",
                        header=header,
                        compression="gzip"
                    )
                    
                    # Update manifest
                    for file_name in batch_files:
                        manifest_data["processed_files"][file_name] = {
                            "part": part_num,
                            "timestamp": timestamp
                        }
                    
                    # Update output part info
                    part_key = str(part_num)
                    if part_key not in manifest_data["output_parts"]:
                        manifest_data["output_parts"][part_key] = {
                            "row_count": 0,
                            "files": []
                        }
                    
                    manifest_data["output_parts"][part_key]["row_count"] += len(df_batch)
                    manifest_data["output_parts"][part_key]["files"].extend(batch_files)
                    manifest_data["last_part"] = part_num
                    manifest_data["total_rows"] += len(df_batch)
                    
                    total_written += len(df_batch)
                    
                    # Save manifest periodically
                    save_manifest(manifest_file, manifest_data)
                    
                    print(f"\n💾 Batch written to part{part_num:04d} ({len(df_batch)} rows)")
                    
                    batch.clear()  # Free memory after write
                    batch_files.clear()
            pbar.update(1)


# ===================================================
# Write Remaining Batches (If Any)
# ===================================================
if batch:
    df_batch = pd.concat(batch, ignore_index=True)
    
    # Get the appropriate output file
    out_file, part_num, create_new = get_next_output_file(
        out_dir, type_, class_, chunk, manifest_data
    )
    
    header = create_new or not os.path.exists(out_file)
    df_batch.to_csv(
        out_file,
        sep="\t",
        index=False,
        mode="a",
        header=header,
        compression="gzip"
    )
    
    # Update manifest
    for file_name in batch_files:
        manifest_data["processed_files"][file_name] = {
            "part": part_num,
            "timestamp": timestamp
        }
    
    # Update output part info
    part_key = str(part_num)
    if part_key not in manifest_data["output_parts"]:
        manifest_data["output_parts"][part_key] = {
            "row_count": 0,
            "files": []
        }
    
    manifest_data["output_parts"][part_key]["row_count"] += len(df_batch)
    manifest_data["output_parts"][part_key]["files"].extend(batch_files)
    manifest_data["last_part"] = part_num
    manifest_data["total_rows"] += len(df_batch)
    
    total_written += len(df_batch)
    
    # Save final manifest
    save_manifest(manifest_file, manifest_data)


# ===================================================
# Summary and Completion Message
# ===================================================
print(f"\n{'='*60}")
print(f"✅ Processing Complete!")
print(f"{'='*60}")
print(f"   Total rows written: {total_written:,}")
print(f"   Files processed: {len(manifest_data['processed_files'])}")
print(f"   Output parts created: {len(manifest_data['output_parts'])}")
print(f"   Manifest saved: {manifest_file}")

# List output files
output_files = sorted(glob.glob(os.path.join(out_dir, f"{type_}_{class_}*_predictions_part*.txt.gz")))
if output_files:
    print(f"\n📁 Output files:")
    for of in output_files:
        size_mb = os.path.getsize(of) / (1024 * 1024)
        part_match = re.search(r'part(\d+)', of)
        if part_match:
            part_num = part_match.group(1)
            part_info = manifest_data["output_parts"].get(part_num, {})
            row_count = part_info.get("row_count", "unknown")
            print(f"   {os.path.basename(of)}: {size_mb:.1f} MB, {int(row_count):,} rows")