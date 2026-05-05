#!/usr/bin/env python3
"""
Quick check: do organism names from the IEDB FASTA headers
match source_organism in the strong mimicry file?
"""
import pandas as pd

FASTA_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/fasta/iedb_epitopes_mhc_i.fasta"
STRONG_FILE = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_i/parsed/iedb_mhci_strong_mimicry.csv.gz"

# ====================================================
# Parse FASTA headers
# ====================================================
print("Parsing FASTA headers...")
# Format: >structure_id|sequence|organism|HLA
fasta_records = {}  # structure_id → organism
with open(FASTA_FILE) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        parts = line[1:].strip().split("|")
        if len(parts) >= 3:
            sid = parts[0]
            org = parts[2].replace("_", " ")
            fasta_records[sid] = org

fasta_orgs = set(fasta_records.values())
fasta_sids = set(fasta_records.keys())
print(f"  Unique structure_ids: {len(fasta_sids):,}")
print(f"  Unique organisms: {len(fasta_orgs)}")

# ====================================================
# Load strong mimicry
# ====================================================
print("\nLoading strong mimicry...")
strong = pd.read_csv(STRONG_FILE, low_memory=False)
strong_orgs = set(strong["source_organism"].unique())
strong_sids = set(strong["structure_id"].astype(str).unique())
print(f"  Unique structure_ids: {len(strong_sids):,}")
print(f"  Unique organisms: {len(strong_orgs)}")

# ====================================================
# Compare organism names
# ====================================================
print(f"\n{'=' * 60}")
print("Organism name comparison")
print(f"{'=' * 60}")

in_both = fasta_orgs & strong_orgs
in_fasta_only = fasta_orgs - strong_orgs
in_strong_only = strong_orgs - fasta_orgs

print(f"\n  In both:        {len(in_both)}")
print(f"  FASTA only:     {len(in_fasta_only)}")
print(f"  Strong only:    {len(in_strong_only)}")

if in_strong_only:
    print(f"\n  Organisms in strong but NOT in FASTA ({len(in_strong_only)}):")
    for org in sorted(in_strong_only):
        n = strong[strong["source_organism"] == org]["pair_id"].nunique()
        print(f"    {org}  ({n} pairs)")

# ====================================================
# Compare structure_ids
# ====================================================
print(f"\n{'=' * 60}")
print("structure_id comparison")
print(f"{'=' * 60}")

sid_both = fasta_sids & strong_sids
sid_fasta_only = fasta_sids - strong_sids
sid_strong_only = strong_sids - fasta_sids

print(f"\n  In both:        {len(sid_both):,}")
print(f"  FASTA only:     {len(sid_fasta_only):,}")
print(f"  Strong only:    {len(sid_strong_only):,}")

if sid_strong_only:
    print(f"\n  Sample structure_ids in strong but NOT in FASTA:")
    for sid in sorted(sid_strong_only)[:10]:
        org = strong[strong["structure_id"].astype(str) == sid]["source_organism"].iloc[0]
        print(f"    {sid}  (organism: {org})")

# ====================================================
# For organisms only in strong, check what their
# structure_ids map to in the FASTA
# ====================================================
if in_strong_only:
    print(f"\n{'=' * 60}")
    print("Tracing mismatched organisms via structure_id")
    print(f"{'=' * 60}")
    
    for org in sorted(in_strong_only):
        sids = set(strong[strong["source_organism"] == org]["structure_id"].astype(str).unique())
        fasta_mapped_orgs = set()
        for sid in sids:
            if sid in fasta_records:
                fasta_mapped_orgs.add(fasta_records[sid])
        
        if fasta_mapped_orgs:
            print(f"\n  Strong: {org}")
            print(f"    FASTA maps to: {fasta_mapped_orgs}")
        else:
            print(f"\n  Strong: {org}")
            print(f"    NO structure_id match in FASTA (sids: {sorted(sids)[:3]})")

print(f"\n✅ Done.")