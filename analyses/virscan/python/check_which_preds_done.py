import os
import glob
import shutil
import math

# === CONFIG ===
base_dir = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks"
out_dir = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers"
processed_files_dir = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/Type1_NR/All_types"
files_per_dir = 50

# === Step 1: Collect all peptide FASTA files ===
all_fastas = sorted(glob.glob(os.path.join(base_dir, "*.fasta")))
print(f"Found {len(all_fastas)} FASTA files total.")

# === Step 2: Identify already processed peptides (based on .xls outputs) ===
processed_files = [os.path.basename(f) for f in glob.glob(os.path.join(processed_files_dir, "*.xls"))] 
processed_peptides = [str(f).split("_matched_pairs_chunk_")[1].replace(".xls", "") for f in processed_files]
processed_peptides = list(set(processed_peptides))  # unique

print(f"Detected {len(processed_peptides)} processed peptide chunks.")

# === Step 3: Map FASTA files by chunk id (numeric suffix) ===
fasta_map = {}
for fasta in all_fastas:
    name = os.path.basename(fasta)
    try:
        chunk_id = name.split("matched_pairs_chunk_")[1].replace(".fasta", "")
        fasta_map[chunk_id] = fasta
    except IndexError:
        print(f"Warning: Skipping file '{name}' as it doesn't match the expected format")

# === Step 4: Determine which are processed vs unprocessed ===
processed_fastas = [fasta_map[c] for c in processed_peptides if c in fasta_map]
unprocessed_fastas = [f for c, f in fasta_map.items() if c not in processed_peptides]

print(f"{len(processed_fastas)} processed FASTAs will go into chunk_1.")
print(f"{len(unprocessed_fastas)} remaining FASTAs will be distributed into new chunks.")

# === Step 5: Create output directories ===
parent_out = os.path.join(out_dir, "split_chunks")
os.makedirs(parent_out, exist_ok=True)

# Chunk 1 (processed ones)
chunk_1_dir = os.path.join(parent_out, "chunk_1")
os.makedirs(chunk_1_dir, exist_ok=True)

for f in processed_fastas[:files_per_dir]:  # keep only 50
    shutil.copy2(f, chunk_1_dir)

# === Step 6: Split the remaining into groups of 50 ===
remaining = unprocessed_fastas
num_extra_chunks = math.ceil(len(remaining) / files_per_dir)

for i in range(num_extra_chunks):
    chunk_dir = os.path.join(parent_out, f"chunk_{i + 2}")
    os.makedirs(chunk_dir, exist_ok=True)

    subset = remaining[i * files_per_dir : (i + 1) * files_per_dir]
    for f in subset:
        shutil.copy2(f, chunk_dir)

    print(f"Created {chunk_dir} with {len(subset)} FASTA files.")

print("\n✅ Split complete.")
print(f"Parent output directory: {parent_out}")
