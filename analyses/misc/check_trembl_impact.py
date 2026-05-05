#!/usr/bin/env python3
"""
Check the impact of filtering to Swiss-Prot only entries
on the current pipeline results.
"""
import pandas as pd

BLAST_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla.csv.gz"
STRONG_FILE = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_i/parsed/iedb_mhci_strong_mimicry.csv.gz"
PROTEOME_FASTA = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta"

# ====================================================
# Step 1: Count Swiss-Prot vs TrEMBL in proteome FASTA
# ====================================================
print("=" * 60)
print("Proteome FASTA composition")
print("=" * 60)

sp_acs = set()
tr_acs = set()
with open(PROTEOME_FASTA) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        if line.startswith(">sp|"):
            ac = line.split("|")[1]
            sp_acs.add(ac)
        elif line.startswith(">tr|"):
            ac = line.split("|")[1]
            tr_acs.add(ac)

print(f"  Swiss-Prot (sp|): {len(sp_acs):,}")
print(f"  TrEMBL (tr|):     {len(tr_acs):,}")
print(f"  Total:            {len(sp_acs) + len(tr_acs):,}")
print(f"  TrEMBL fraction:  {len(tr_acs) / (len(sp_acs) + len(tr_acs)) * 100:.1f}%")

# ====================================================
# Step 2: Impact on BLAST pairs
# ====================================================
print(f"\n{'=' * 60}")
print("BLAST pairs impact")
print("=" * 60)

blast = pd.read_csv(BLAST_FILE, low_memory=False)
print(f"  Total pairs: {len(blast):,}")

blast["is_sp"] = blast["hu_prot_id"].isin(sp_acs)
blast["is_tr"] = blast["hu_prot_id"].isin(tr_acs)
blast["is_neither"] = ~blast["is_sp"] & ~blast["is_tr"]

n_sp = blast["is_sp"].sum()
n_tr = blast["is_tr"].sum()
n_neither = blast["is_neither"].sum()

print(f"  Swiss-Prot hits: {n_sp:,} ({n_sp/len(blast)*100:.1f}%)")
print(f"  TrEMBL hits:     {n_tr:,} ({n_tr/len(blast)*100:.1f}%)")
print(f"  Neither:         {n_neither:,}")

if n_neither > 0:
    print(f"    Sample 'neither' hu_prot_ids: {blast[blast['is_neither']]['hu_prot_id'].head(5).tolist()}")

# Unique epitopes lost
blast_sp = blast[blast["is_sp"]]
all_epitopes = set(blast["structure_id"].unique())
sp_epitopes = set(blast_sp["structure_id"].unique())
lost_epitopes = all_epitopes - sp_epitopes
print(f"\n  Unique epitopes (structure_id):")
print(f"    Total:           {len(all_epitopes):,}")
print(f"    In sp| hits:     {len(sp_epitopes):,}")
print(f"    Lost (tr| only): {len(lost_epitopes):,}")

# Unique human proteins
all_prots = set(blast["hu_prot_id"].unique())
sp_prots = set(blast_sp["hu_prot_id"].unique())
print(f"\n  Unique human proteins hit:")
print(f"    Total:     {len(all_prots):,}")
print(f"    sp| only:  {len(sp_prots):,}")
print(f"    Lost:      {len(all_prots - sp_prots):,}")

# ====================================================
# Step 3: Impact on strong mimicry
# ====================================================
print(f"\n{'=' * 60}")
print("Strong mimicry impact")
print("=" * 60)

strong = pd.read_csv(STRONG_FILE, low_memory=False)
print(f"  Total strong pairs: {len(strong):,}")

strong["is_sp"] = strong["hu_prot_id"].isin(sp_acs) if "hu_prot_id" in strong.columns else False

n_sp_strong = strong["is_sp"].sum()
n_tr_strong = len(strong) - n_sp_strong

print(f"  Swiss-Prot: {n_sp_strong:,} ({n_sp_strong/len(strong)*100:.1f}%)")
print(f"  TrEMBL:     {n_tr_strong:,} ({n_tr_strong/len(strong)*100:.1f}%)")

# Check if TrEMBL pairs are redundant with sp| pairs
# (same viral epitope + same mhc_allele, just different human protein)
if "viral_peptide" in strong.columns and "mhc_allele" in strong.columns:
    strong_sp = strong[strong["is_sp"]]
    strong_tr = strong[~strong["is_sp"]]

    # For TrEMBL pairs, check if the same (viral_peptide, mhc_allele)
    # also has a Swiss-Prot hit
    sp_keys = set(
        strong_sp[["viral_peptide", "mhc_allele"]]
        .apply(tuple, axis=1)
    )
    tr_keys = set(
        strong_tr[["viral_peptide", "mhc_allele"]]
        .apply(tuple, axis=1)
    )

    redundant = tr_keys & sp_keys
    unique_tr = tr_keys - sp_keys

    print(f"\n  TrEMBL (viral_peptide, mhc_allele) pairs:")
    print(f"    Also have sp| hit (redundant): {len(redundant):,}")
    print(f"    TrEMBL-only (would be lost):   {len(unique_tr):,}")

# ====================================================
# Step 4: Top TrEMBL-only proteins
# ====================================================
print(f"\n{'=' * 60}")
print("Top proteins that would be lost (TrEMBL-only)")
print("=" * 60)

if "hu_prot_name" in strong.columns:
    tr_strong = strong[~strong["is_sp"]]
    tr_protein_counts = tr_strong["hu_prot_name"].value_counts().head(15)
    if not tr_protein_counts.empty:
        print(f"  Top 15 TrEMBL proteins by pair count:")
        for prot, n in tr_protein_counts.items():
            print(f"    {prot}: {n:,} pairs")
    else:
        print("  No TrEMBL proteins in strong mimicry")

print(f"\n✅ Done.")