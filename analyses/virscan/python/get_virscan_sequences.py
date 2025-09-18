#!/usr/bin/env python3
"""
Extract peptide sequences from UniProt FASTA files using coordinates
and write them into a FASTA file for BLAST.
"""

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# =====================
# CONFIG
# =====================
virscan_csv = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/peptide_metadata.csv"
human_fasta = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_sprot.fasta"
virus_fasta = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_viruses_sprot.fasta"
out_fasta = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/VirScan_peptides.fasta"

# =====================
# Load VirScan metadata
# =====================
df = pd.read_csv(virscan_csv)

# =====================
# Load FASTA sequences
# =====================
def load_fasta_to_dict(fasta_file):
    """Return dict {UniProt_acc: SeqRecord}"""
    records = {}
    for rec in SeqIO.parse(fasta_file, "fasta"):
        # UniProt FASTA header usually like:
        # >sp|P12345|PROT_HUMAN Some protein name
        acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
        records[acc] = rec
    return records

human_dict = load_fasta_to_dict(human_fasta)
virus_dict = load_fasta_to_dict(virus_fasta)

# Merge both dicts
seq_dict = {**human_dict, **virus_dict}

# =====================
# Extract peptides
# =====================
peptide_records = []

for _, row in df.iterrows():
    acc = str(row["UniProt_acc"])
    start = int(row["start"])
    end = int(row["end"])
    pep_id = str(row["pep_id"])

    if acc not in seq_dict:
        print(f"⚠️ UniProt accession {acc} not found in FASTAs")
        continue

    full_seq = seq_dict[acc].seq
    # UniProt coordinates are usually 1-based inclusive
    peptide_seq = full_seq[start-1:end]

    # Create FASTA entry
    header = f"{pep_id}|{acc}|{start}-{end}"
    rec = SeqRecord(Seq(str(peptide_seq)), id=header, description="")
    peptide_records.append(rec)

# =====================
# Write FASTA
# =====================
SeqIO.write(peptide_records, out_fasta, "fasta")
print(f"✅ Wrote {len(peptide_records)} peptides to {out_fasta}")
