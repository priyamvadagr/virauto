#!/usr/bin/env python3
"""
filter_blast_swissprot.py

Filter BLAST mimicry pairs to Swiss-Prot (reviewed) entries only.
Removes TrEMBL (unreviewed) hits to eliminate isoform/fragment redundancy.

Input:
    - iedb_mhci_pairs_4digit_hla.csv.gz (BLAST output with all hits)
    - uniprot_human_all.fasta (to identify sp| vs tr| entries)

Output:
    - iedb_mhci_pairs_4digit_hla_swissprot.csv.gz

Usage:
    python filter_blast_swissprot.py
"""

import pandas as pd
import os

# ====================================================
# Config
# ====================================================
BASE_DIR = "/ix/djishnu/Priyamvada/virauto"
BLAST_FILE = os.path.join(BASE_DIR, "data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla.csv.gz")
PROTEOME_FASTA = os.path.join(BASE_DIR, "data/refs/uniprot/uniprot_human_all.fasta")
OUTPUT_FILE = os.path.join(BASE_DIR, "data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla_swissprot.csv.gz")

# ====================================================
# Step 1: Identify Swiss-Prot ACs from FASTA
# ====================================================
print("=" * 60)
print("Step 1: Identifying Swiss-Prot entries from FASTA")
print("=" * 60)

sp_acs = set()
with open(PROTEOME_FASTA) as fh:
    for line in fh:
        if line.startswith(">sp|"):
            ac = line.split("|")[1]
            sp_acs.add(ac)

print(f"  Swiss-Prot ACs: {len(sp_acs):,}")

# ====================================================
# Step 2: Filter BLAST pairs
# ====================================================
print(f"\n{'=' * 60}")
print("Step 2: Filtering BLAST pairs")
print("=" * 60)

blast = pd.read_csv(BLAST_FILE, low_memory=False)
print(f"  Input pairs: {len(blast):,}")
print(f"  Unique epitopes: {blast['structure_id'].nunique():,}")
print(f"  Unique human proteins: {blast['hu_prot_id'].nunique():,}")
print(f"  Unique HLA alleles: {blast['mhc_allele'].nunique()}")

# Filter
blast_sp = blast[blast["hu_prot_id"].isin(sp_acs)].copy()

print(f"\n  After Swiss-Prot filter:")
print(f"  Pairs: {len(blast_sp):,} ({len(blast_sp)/len(blast)*100:.1f}%)")
print(f"  Unique epitopes: {blast_sp['structure_id'].nunique():,}")
print(f"  Unique human proteins: {blast_sp['hu_prot_id'].nunique():,}")
print(f"  Unique HLA alleles: {blast_sp['mhc_allele'].nunique()}")

# Report losses
lost_epitopes = set(blast["structure_id"].unique()) - set(blast_sp["structure_id"].unique())
lost_prots = set(blast["hu_prot_id"].unique()) - set(blast_sp["hu_prot_id"].unique())
print(f"\n  Lost:")
print(f"    Epitopes (no sp| hit): {len(lost_epitopes):,}")
print(f"    Human proteins:        {len(lost_prots):,}")
print(f"    Pairs:                 {len(blast) - len(blast_sp):,}")

# ====================================================
# Step 3: Save
# ====================================================
print(f"\n{'=' * 60}")
print("Step 3: Saving filtered output")
print("=" * 60)

blast_sp.to_csv(OUTPUT_FILE, index=False, compression="gzip")
print(f"  Output: {OUTPUT_FILE}")
print(f"  Rows: {len(blast_sp):,}")

# ====================================================
# Summary of what to re-run
# ====================================================
print(f"\n{'=' * 60}")
print("Next steps")
print("=" * 60)
print(f"""
  1. Update INPUT_FILE in prepare_netmhcpan_iedb.py to:
     {OUTPUT_FILE}

  2. Re-run the pipeline:
     python prepare_netmhcpan_iedb.py
     sbatch submit_netmhcpan_iedb.sh
     python parse_netmhcpan_iedb_results.py

  3. Then re-run downstream analyses:
     python compute_mimetopes_per_protein.py
     python compute_protein_enrichment.py
     python compute_domain_enrichment.py
     python plot_mimicry_summary.py
""")

print("✅ Done.")