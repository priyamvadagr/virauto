#!/usr/bin/env python3
"""
Split HLA allele file into A, B, C, and Others groups,
and write each group into its own directory with chunks of 30 alleles each.

Output filenames are zero-padded: allele_chunk_000.txt, allele_chunk_001.txt, etc.
"""

import os

# =============================================================================
# Configuration
# =============================================================================
input_file = "/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Class_I_NR/HLA_all_NR_alleles_with_names.txt"
base_output_dir = "/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Class_I_NR"
chunk_size = 30

# Define groups
groups = {
    "HLA-A": [],
    "HLA-B": [],
    "HLA-C": [],
    "HLA-others": []
}

# =============================================================================
# Read and categorize alleles
# =============================================================================
with open(input_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith("HLA-A"):
            groups["HLA-A"].append(line)
        elif line.startswith("HLA-B"):
            groups["HLA-B"].append(line)
        elif line.startswith("HLA-C"):
            groups["HLA-C"].append(line)
        else:
            groups["HLA-others"].append(line)

# =============================================================================
# Write each group into its own directory, chunked by 30 alleles
# =============================================================================
for group_name, lines in groups.items():
    if not lines:
        print(f"⚠️ No alleles found for {group_name}, skipping.")
        continue
    
    group_dir = os.path.join(base_output_dir, f"{group_name}_chunks")
    os.makedirs(group_dir, exist_ok=True)
    
    num_chunks = (len(lines) - 1) // chunk_size + 1
    
    for i in range(num_chunks):
        start = i * chunk_size
        end = start + chunk_size
        chunk = lines[start:end]
        
        # Zero-padded filename
        output_file = os.path.join(group_dir, f"allele_chunk_{i:03d}.txt")
        with open(output_file, "w") as out:
            out.write("\n".join(chunk) + "\n")
    
    print(f"✅ {group_name}: {len(lines)} alleles split into {num_chunks} chunks in {group_dir}/")

print("\n🎉 All allele groups processed successfully!")
