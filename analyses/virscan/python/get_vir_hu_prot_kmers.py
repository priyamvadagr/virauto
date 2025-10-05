#!/usr/bin/env python3
import pandas as pd
from Bio.Seq import Seq

# Load BLAST hits with viral/human alignments (filtered)
blast_hits = pd.read_csv(
    "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/similarity_filtered_blast_out.csv"
)

# Load peptide FASTAs (dict by ID → sequence)
from Bio import SeqIO

viral_fasta = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/VirScan_peptides_api.fasta"
human_fasta = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta"

viral_dict = {}
for rec in SeqIO.parse(viral_fasta, "fasta"):
        acc = rec.id.split("|")[0]  # adjust if IDs contain extra parts
        viral_dict[acc] = str(rec.seq)
   
human_dict = {}
for rec in SeqIO.parse(human_fasta, "fasta"):
    # UniProt FASTA headers are like >sp|P12345|PROT_HUMAN
    acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
    human_dict[acc] = str(rec.seq)

# Parameters
k = 9
records = []

# Iterate over BLAST hits
for _, row in blast_hits.iterrows():
    pep_id = str(row["qseqid"]).split("|")[0]   # VirScan peptide ID
    vir_protein = str(row["qseqid"]).split("|")[1]
    human_id = str(row["sseqid"]).split("|")[1] if "|" in str(row["sseqid"]) else str(row["sseqid"])
    vseq = viral_dict.get(pep_id)
    hseq = human_dict.get(human_id)

    if not vseq or not hseq:
        continue

    # Extract aligned region from human
    start, end = int(row["sstart"]), int(row["send"])
    if start <= end:
        h_aligned = hseq[start-1:end]   # 1-based inclusive indexing
    else:
        continue  # reverse alignments rare in protein BLAST

    # Chop into k-mers
    v_kmers = [vseq[i:i+k] for i in range(len(vseq)-k+1)]
    h_kmers = [h_aligned[i:i+k] for i in range(len(h_aligned)-k+1)]

    # Pair positionally (truncate to shortest)
    min_len = min(len(v_kmers), len(h_kmers))
    for i in range(min_len):
        records.append({
            "pep_id": pep_id,
            "viral_protein": vir_proteinf,
            "human_protein": human_id,
            "offset": i+1,
            "viral_kmer": v_kmers[i],
            "human_kmer": h_kmers[i]
        })

# Save paired kmers
out_csv = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/paired_viral_human_9mers.csv"
pd.DataFrame(records).to_csv(out_csv, index=False)

print(f"✅ Paired k-mers written to {out_csv}")
