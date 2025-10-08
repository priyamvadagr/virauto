#!/usr/bin/env python3
"""
======================================================================
Script: generate_paired_viral_human_kmers.py
Author: Priyamvada Guha Roy
Description:
    This script pairs viral peptide sequences with their corresponding 
    aligned human protein regions based on BLAST results, then generates 
    positionally matched k-mer pairs (default: 9-mers) for downstream 
    similarity or mimicry analyses.

Workflow:
    1. Load filtered BLAST alignment results between viral peptides 
       (VirScan) and human proteins.
    2. Load viral and human protein FASTA files into dictionaries.
    3. For each BLAST hit:
        - Retrieve the viral peptide and aligned human subsequence.
        - Extract the aligned region (based on sstart–send coordinates).
        - Slice both viral and human sequences into k-mers (e.g., 9-mers).
        - Pair viral and human k-mers positionally up to the shortest length.
    4. Save all paired k-mers with metadata to a CSV file.

Inputs:
    - similarity_filtered_blast_out.csv: Filtered BLAST alignments.
    - VirScan_peptides_api.fasta: Viral peptide sequences.
    - uniprot_human_all.fasta: Human protein sequences (UniProt reference).

Outputs:
    - paired_viral_human_<k>mers.csv: Table of paired viral–human k-mers with columns:
        [pep_id, viral_protein, human_protein, offset, viral_kmer, human_kmer]

Parameters:
    - k: k-mer length (default = 9; can be set via --kmer_size)

Dependencies:
    - pandas
    - biopython

Usage:
    python generate_paired_viral_human_kmers.py --kmer_size 9
======================================================================
"""

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# ====================================================
# Command-line arguments
# ====================================================
parser = argparse.ArgumentParser(description="Generate paired viral–human k-mers from BLAST hits.")
parser.add_argument(
    "--kmer_size", "-k",
    type=int,
    default=9,
    help="Length of k-mers to generate (default: 9)"
)
args = parser.parse_args()
k = args.kmer_size

# ====================================================
# File paths
# ====================================================
blast_hits_file = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/similarity_filtered_blast_out.csv"
viral_fasta = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/VirScan_peptides_api.fasta"
human_fasta = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta"
out_dir = f"/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/{k}_mers"

# ====================================================
# Load BLAST hits
# ====================================================
blast_hits = pd.read_csv(blast_hits_file)

# ====================================================
# Load viral and human FASTAs into dictionaries
# ====================================================
viral_dict = {
    rec.id.split("|")[0]: str(rec.seq)
    for rec in SeqIO.parse(viral_fasta, "fasta")
}

human_dict = {}
for rec in SeqIO.parse(human_fasta, "fasta"):
    acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
    human_dict[acc] = str(rec.seq)

# ====================================================
# Generate paired k-mers
# ====================================================
records = []

for _, row in blast_hits.iterrows():
    pep_id = str(row["qseqid"]).split("|")[0]   # VirScan peptide ID
    viral_protein = str(row["qseqid"]).split("|")[1]
    human_id = str(row["sseqid"]).split("|")[1] if "|" in str(row["sseqid"]) else str(row["sseqid"])

    vseq = viral_dict.get(pep_id)
    hseq = human_dict.get(human_id)

    if not vseq or not hseq:
        continue

    # Extract aligned region from human
    start, end = int(row["sstart"]), int(row["send"])
    if start <= end:
        h_aligned = hseq[start-1:end]  # 1-based inclusive indexing
    else:
        continue  # Skip reverse alignments (rare for proteins)

    # Chop into k-mers
    v_kmers = [vseq[i:i+k] for i in range(len(vseq) - k + 1)]
    h_kmers = [h_aligned[i:i+k] for i in range(len(h_aligned) - k + 1)]

    # Pair positionally (truncate to shortest)
    min_len = min(len(v_kmers), len(h_kmers))
    for i in range(min_len):
        records.append({
            "pep_id": pep_id,
            "viral_protein": viral_protein,
            "human_protein": human_id,
            "offset": i + 1,
            "viral_kmer": v_kmers[i],
            "human_kmer": h_kmers[i]
        })

# ====================================================
# Save results
# ====================================================
os.makedirs(out_dir, exist_ok=True)
out_csv = os.path.join(out_dir, f"paired_viral_human_{k}mers.csv")
pd.DataFrame(records).to_csv(out_csv, index=False)
print(f"✅ Paired {len(records)} k-mers written to {out_csv}")
