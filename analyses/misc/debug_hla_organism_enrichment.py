#!/usr/bin/env python3
"""
Debug: trace through Figures 1 and 2 logic step by step
"""
import pandas as pd

STRONG_FILE = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_i/parsed/iedb_mhci_strong_mimicry.csv.gz"
BLAST_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla.csv.gz"
PAIR_MAP_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/netmhcpan/mhc_i/pair_id_mapping.csv.gz"
FASTA_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/fasta/iedb_epitopes_mhc_i.fasta"

# ====================================================
# Load files
# ====================================================
print("=" * 60)
print("Loading files")
print("=" * 60)

strong = pd.read_csv(STRONG_FILE, low_memory=False)
blast = pd.read_csv(BLAST_FILE, low_memory=False)
pair_map = pd.read_csv(PAIR_MAP_FILE, low_memory=False)

print(f"  Strong: {len(strong):,} rows, columns: {strong.columns.tolist()}")
print(f"  BLAST:  {len(blast):,} rows, columns: {blast.columns.tolist()}")
print(f"  PairMap: {len(pair_map):,} rows, columns: {pair_map.columns.tolist()}")

# ====================================================
# FIGURE 1 DEBUG: HLA conversion rate
# ====================================================
print(f"\n{'=' * 60}")
print("FIGURE 1: HLA conversion rate")
print("=" * 60)

# Step 1: Does BLAST file have pair_id?
print(f"\n  Step 1: Does BLAST file have pair_id?")
print(f"    'pair_id' in blast columns: {'pair_id' in blast.columns}")

# Step 2: Reconstruct pair_key and merge pair_ids
print(f"\n  Step 2: Merge pair_ids into BLAST file")
blast["pair_key"] = (
    blast["structure_id"].astype(str) + "_" +
    blast["hu_prot_id"].astype(str) + "_" +
    blast["viral_sequence"].astype(str) + "_" +
    blast["human_mimic_sequence"].astype(str)
)
print(f"    BLAST unique pair_keys: {blast['pair_key'].nunique():,}")
print(f"    PairMap unique pair_keys: {pair_map['pair_key'].nunique():,}")

# Check overlap
blast_keys = set(blast["pair_key"].unique())
map_keys = set(pair_map["pair_key"].unique())
print(f"    Overlap: {len(blast_keys & map_keys):,}")
print(f"    BLAST only: {len(blast_keys - map_keys):,}")
print(f"    PairMap only: {len(map_keys - blast_keys):,}")

key_to_pair_id = pair_map.drop_duplicates("pair_key").set_index("pair_key")["pair_id"]
blast["pair_id"] = blast["pair_key"].map(key_to_pair_id)
n_mapped = blast["pair_id"].notna().sum()
n_unmapped = blast["pair_id"].isna().sum()
print(f"    Mapped: {n_mapped:,}, Unmapped: {n_unmapped:,}")

blast_with_pid = blast.dropna(subset=["pair_id"])

# Step 3: Count per HLA in BLAST (denominator)
print(f"\n  Step 3: Count unique pair_ids per HLA (denominator)")
pairs_per_hla_input = blast_with_pid.groupby("mhc_allele")["pair_id"].nunique()
print(f"    HLA alleles: {len(pairs_per_hla_input)}")
print(f"    Top 5:")
for allele, n in pairs_per_hla_input.nlargest(5).items():
    print(f"      {allele}: {n:,} unique pair_ids")

# Step 4: Count per HLA in strong (numerator)
print(f"\n  Step 4: Count unique pair_ids per HLA in strong (numerator)")
strong_per_hla = strong.groupby("mhc_allele")["pair_id"].nunique()
print(f"    HLA alleles: {len(strong_per_hla)}")
print(f"    Top 5:")
for allele, n in strong_per_hla.nlargest(5).items():
    print(f"      {allele}: {n:,} unique pair_ids")

# Step 5: Compute fractions
print(f"\n  Step 5: Compute fractions")
hla_conv = pairs_per_hla_input.rename("input_pairs").to_frame().join(
    strong_per_hla.rename("strong_pairs"), how="left"
).fillna(0)
hla_conv["strong_pairs"] = hla_conv["strong_pairs"].astype(int)
hla_conv["fraction"] = hla_conv["strong_pairs"] / hla_conv["input_pairs"]

print(f"    Top 10 by fraction:")
top = hla_conv[hla_conv["strong_pairs"] > 0].nlargest(10, "fraction")
for allele, row in top.iterrows():
    print(f"      {allele}: {int(row['strong_pairs'])}/{int(row['input_pairs'])} = {row['fraction']:.4f}")

all_one = (hla_conv[hla_conv["strong_pairs"] > 0]["fraction"] == 1.0).all()
print(f"\n    ALL fractions == 1.0? {all_one}")
if all_one:
    print(f"    ⚠️ BUG: every strong pair_id appears exactly once per HLA in both files")
    # Check if pair_ids are the same sets
    for allele in strong_per_hla.nlargest(3).index:
        strong_pids = set(strong[strong["mhc_allele"] == allele]["pair_id"].unique())
        blast_pids = set(blast_with_pid[blast_with_pid["mhc_allele"] == allele]["pair_id"].unique())
        print(f"\n      {allele}:")
        print(f"        Strong pair_ids: {len(strong_pids)}")
        print(f"        BLAST pair_ids:  {len(blast_pids)}")
        print(f"        Strong - BLAST:  {len(strong_pids - blast_pids)}")
        print(f"        BLAST - Strong:  {len(blast_pids - strong_pids)}")

# ====================================================
# FIGURE 2 DEBUG: Organism conversion rate
# ====================================================
print(f"\n{'=' * 60}")
print("FIGURE 2: Organism conversion rate")
print("=" * 60)

# Step 1: Parse FASTA
print(f"\n  Step 1: Parse FASTA")
fasta_orgs = {}  # sid → organism
with open(FASTA_FILE) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        parts = line[1:].strip().split("|")
        if len(parts) >= 3:
            sid = parts[0]
            org = parts[2].replace("_", " ")
            fasta_orgs[sid] = org

print(f"    Unique structure_ids: {len(fasta_orgs):,}")
print(f"    Sample sids: {list(fasta_orgs.keys())[:3]}")
print(f"    Sample orgs: {list(fasta_orgs.values())[:3]}")

# Step 2: Count per organism in FASTA (denominator)
from collections import Counter
org_counts_fasta = Counter(fasta_orgs.values())
print(f"\n  Step 2: FASTA epitopes per organism (denominator)")
print(f"    Unique organisms: {len(org_counts_fasta)}")
print(f"    Top 5:")
for org, n in Counter(org_counts_fasta).most_common(5):
    print(f"      {org}: {n}")

# Step 3: Count per organism in strong (numerator)
print(f"\n  Step 3: Strong mimicry structure_ids per organism (numerator)")
print(f"    structure_id dtype: {strong['structure_id'].dtype}")
print(f"    Sample values: {strong['structure_id'].head(3).tolist()}")

strong_epi_per_org = strong.groupby("source_organism")["structure_id"].nunique()
print(f"    Unique organisms: {len(strong_epi_per_org)}")
print(f"    Top 5:")
for org, n in strong_epi_per_org.nlargest(5).items():
    print(f"      {org}: {n}")

# Step 4: Check join alignment
print(f"\n  Step 4: Join alignment")
epitopes_per_org_fasta = pd.Series(org_counts_fasta, name="fasta_epitopes")

org_conv = epitopes_per_org_fasta.to_frame().join(
    strong_epi_per_org.rename("strong_epitopes"), how="left"
).fillna(0)
org_conv["strong_epitopes"] = org_conv["strong_epitopes"].astype(int)
org_conv["fraction"] = org_conv["strong_epitopes"] / org_conv["fasta_epitopes"]

print(f"    Top 10 by fraction:")
top_org = org_conv[org_conv["strong_epitopes"] > 0].nlargest(10, "fraction")
for org, row in top_org.iterrows():
    print(f"      {org}: {int(row['strong_epitopes'])}/{int(row['fasta_epitopes'])} = {row['fraction']:.4f}")

all_one_org = (org_conv[org_conv["strong_epitopes"] > 0]["fraction"] == 1.0).all()
print(f"\n    ALL fractions == 1.0? {all_one_org}")

if all_one_org:
    # Deeper check: are the structure_ids identical?
    print(f"\n    ⚠️ Checking structure_id overlap for top organism:")
    top_org_name = strong_epi_per_org.idxmax()
    strong_sids = set(strong[strong["source_organism"] == top_org_name]["structure_id"].astype(str).unique())
    fasta_sids = set(sid for sid, org in fasta_orgs.items() if org == top_org_name)
    print(f"      Organism: {top_org_name}")
    print(f"      Strong structure_ids: {len(strong_sids)}")
    print(f"      FASTA structure_ids:  {len(fasta_sids)}")
    print(f"      Strong - FASTA:       {len(strong_sids - fasta_sids)}")
    print(f"      FASTA - Strong:       {len(fasta_sids - strong_sids)}")
    print(f"      Sample strong sids:   {sorted(strong_sids)[:3]}")
    print(f"      Sample fasta sids:    {sorted(fasta_sids)[:3]}")

print(f"\n✅ Done.")