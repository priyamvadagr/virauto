#!/usr/bin/env python3
"""
===============================================================================
Script: generate_paired_viral_human_kmers_and_fasta_chunks.py
Description:
    Generate paired viral–human k-mers from BLAST hits, excluding exact viral–human sequence matches.
    Write the resulting pairs into chunked FASTA files (alternating VIRAL and HUMAN sequences).

Workflow:
    1. Load filtered BLAST results.
    2. Load viral and human protein FASTA files into dictionaries.
    3. For each BLAST hit:
        - Extract aligned human subsequence.
        - Generate k-mers for viral and human regions.
        - Pair positionally, skipping identical sequences.
    4. Write alternating VIRAL–HUMAN FASTA records in chunks (default 100).

Inputs:
    - similarity_filtered_blast_out.csv
    - VirScan_peptides_api.fasta
    - uniprot_human_all.fasta

Outputs:
    - matched_pairs_chunk_*.fasta (e.g., matched_pairs_chunk_1.fasta)

Parameters:
    --kmer_size <int>    (default: 9)
    --chunk_size <int>   (default: 100)

Dependencies:
    - pandas
    - biopython
===============================================================================
"""

import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ====================================================
# Command-line arguments
# ====================================================
parser = argparse.ArgumentParser(
    description="Generate paired viral–human k-mers (excluding identicals) and FASTA chunks."
)
parser.add_argument("--kmer", "-k", type=int, default=9, help="Length of k-mers (default: 9)")
parser.add_argument("--chunk_size", "-c", type=int, default=100, help="FASTA chunk size (default: 100)")
args = parser.parse_args()

k = args.kmer
chunk_size = args.chunk_size

# ====================================================
# File paths
# ====================================================
blast_hits_file = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/similarity_filtered_blast_out.csv"
viral_fasta = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/VirScan_peptides_api.fasta"
human_fasta = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta"

out_dir = f"/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/{k}_mers"
os.makedirs(out_dir, exist_ok=True)
chunk_dir = os.path.join(out_dir, "chunks")
os.makedirs(chunk_dir, exist_ok=True)

# ====================================================
# Load BLAST hits
# ====================================================
blast_hits = pd.read_csv(blast_hits_file)
print(f"[INFO] Loaded {len(blast_hits):,} BLAST alignments")

# ====================================================
# Load viral and human FASTAs into dictionaries
# ====================================================
viral_dict = {rec.id.split("|")[0]: str(rec.seq) for rec in SeqIO.parse(viral_fasta, "fasta")}
print(f"[INFO] Loaded {len(viral_dict):,} viral peptide sequences")

human_dict = {}
for rec in SeqIO.parse(human_fasta, "fasta"):
    acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
    human_dict[acc] = str(rec.seq)
print(f"[INFO] Loaded {len(human_dict):,} human protein sequences")

# ====================================================
# Generate paired k-mers and FASTA records
# ====================================================
records = []
dropped_pairs = 0

for _, row in blast_hits.iterrows():
    pep_id = str(row["qseqid"]).split("|")[0]
    viral_protein = str(row["qseqid"]).split("|")[1]
    human_id = str(row["sseqid"]).split("|")[1] if "|" in str(row["sseqid"]) else str(row["sseqid"])

    vseq = viral_dict.get(pep_id)
    hseq = human_dict.get(human_id)

    if not vseq or not hseq:
        continue

    start, end = int(row["sstart"]), int(row["send"])
    if start <= end:
        h_aligned = hseq[start - 1:end]
    else:
        continue  # Skip reverse alignments

    v_kmers = [vseq[i:i + k] for i in range(len(vseq) - k + 1)]
    h_kmers = [h_aligned[i:i + k] for i in range(len(h_aligned) - k + 1)]

    min_len = min(len(v_kmers), len(h_kmers))
    for i in range(min_len):
        vk = v_kmers[i].upper()
        hk = h_kmers[i].upper()
        if vk == hk:
            dropped_pairs += 1
            continue  # Skip identical viral-human k-mers

        viral_record = SeqRecord(Seq(vk), id=f"VIRAL_{pep_id}_{i+1}_{viral_protein}", description="")
        human_record = SeqRecord(Seq(hk), id=f"HUMAN_{pep_id}_{i+1}_{human_id}", description="")

        records.append(viral_record)
        records.append(human_record)

print(f"[INFO] Generated {len(records)//2:,} non-identical viral–human pairs "
      f"(dropped {dropped_pairs:,} identical matches)")

# ====================================================
# Write chunked FASTA files
# ====================================================
for i in range(0, len(records), chunk_size):
    chunk = records[i:i + chunk_size]
    chunk_num = i // chunk_size + 1
    output_file = os.path.join(chunk_dir, f"matched_pairs_chunk_{chunk_num}.fasta")
    SeqIO.write(chunk, output_file, "fasta")
    print(f"   [Chunk {chunk_num}] Wrote {len(chunk)} sequences → {output_file}")

print(f"\n✅ Done: {len(records)} total sequences "
      f"({len(records)//2} viral–human pairs) written to {len(os.listdir(chunk_dir))} chunk files.")
