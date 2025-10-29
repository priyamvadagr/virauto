#!/usr/bin/env python3
"""
======================================================================
Script: filter_blast_hits_extract_human_sequences.py
Description:
    This script filters BLAST alignment results between VirScan peptides 
    and the human proteome to identify human protein segments with 
    75–95% sequence similarity to viral peptides and remove uniprot ids
    being dropped from new release. It then extracts the 
    corresponding subsequences from the UniProt human proteome FASTA file 
    and writes them to a new FASTA file for downstream analyses.

Workflow:
    1. Load BLAST output (tab-separated) and assign column names.
    2. Compute coverage-adjusted similarity = pident × (alignment_length / query_length).
    3. Filter hits with 75–95% normalized similarity.
    4. Save filtered hits to CSV.
    5. Parse the human proteome FASTA into a dictionary {UniProt ID → SeqRecord}.
    6. Extract aligned subsequences (start–end) for matching hits.
    7. Write extracted sequences to a new FASTA file.

Inputs:
    - blast_file: VirScan peptide vs human proteome BLAST output (.tsv)
    - human_fasta: UniProt human proteome FASTA file

Outputs:
    - out_df: Filtered BLAST hits (CSV)
    - out_fasta: Extracted human protein subsequences (FASTA)

Dependencies:
    pandas, biopython

Usage:
    python filter_blast_hits_extract_human_sequences.py
======================================================================
"""

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# ====================================================
# Config
# ====================================================
blast_file = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/all_peptides_blast_out.tsv"
out_df = '/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/similarity_filtered_blast_out.csv'
human_fasta = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta"
out_fasta = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/filetered_human_protein_seqs.fasta"
uniprot_filter = '/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/proteins_to_remove_from_UniProtKB.txt'

# ====================================================
# Load blast results and filter to 75-95% similarity
# ====================================================
cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
              "qlen","slen","qstart","qend","sstart","send","evalue","bitscore"]

df = pd.read_csv(blast_file, sep="\t", names=cols)

uniprot_df = pd.read_csv(uniprot_filter , sep = '\t', names = ['hu_prot_id'])

# Compute coverage-normalized similarity
df["similarity"] = df["pident"] * (df["length"] / df["qlen"])

hits = df[(df["similarity"] >= 75) & (df["similarity"] <= 95)]
print(hits.shape)

hits['hu_prot_id'] = hits['sseqid'].str.split("|").str[1]

hits = hits[~hits['hu_prot_id'].isin(uniprot_df['hu_prot_id'])]
print(hits.shape)


# Save to file
hits.to_csv(out_df, index=False)
print(f"✅ Filtered BLAST results saved to {out_df}")
# ====================================================
# Load human proteome into dictionary {UniProt ID → SeqRecord}
# ====================================================
seq_dict = {}
for rec in SeqIO.parse(human_fasta, "fasta"):
    # UniProt FASTA headers are like >sp|P12345|PROT_HUMAN
    acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
    seq_dict[acc] = rec

# ====================================================
# Extract only the aligned subsequence
# ====================================================
records = []
for _, row in hits.iterrows():
    pep_id = row["qseqid"].split("|")[0]   # VirScan peptide ID
    sseqid = row["sseqid"].split("|")[1] if "|" in row["sseqid"] else row["sseqid"]

    if sseqid not in seq_dict:
        print(f"⚠️ {sseqid} not found in FASTA")
        continue

    full_seq = seq_dict[sseqid].seq
    start, end = int(row["sstart"]), int(row["send"])

    if start <= end:
        subseq = full_seq[start-1:end]   # 1-based inclusive indexing
    else:
        subseq = full_seq[end-1:start].reverse_complement()  # rarely for proteins

    header = f"{pep_id}|HUMAN_{sseqid}|{start}-{end}|sim={row['similarity']:.2f}"
    records.append(SeqRecord(Seq(str(subseq)), id=header, description=""))

# ====================================================
# Write FASTA
# ====================================================
SeqIO.write(records, out_fasta, "fasta")
print(f"✅ Extracted {len(records)} aligned subsequences to {out_fasta}")
