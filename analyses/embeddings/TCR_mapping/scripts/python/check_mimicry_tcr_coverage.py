#!/usr/bin/env python3
"""
======================================================================
Script: check_mimicry_tcr_coverage.py

Description:
    Count how many strong mimicry candidates have associated TCR
    sequences in the IEDB receptor_full_v3 export.

    Join chain:
        strong mimicry (pair_id)
            → pair_id_mapping (pair_id → structure_id)
            → iedb_mhci_with_tcr (structure_id → TCR data)

Output:
    Console summary + mimicry_with_tcr.csv.gz (candidates with TCRs)

Usage:
    python check_mimicry_tcr_coverage.py
======================================================================
"""
# %%
import os
import pandas as pd

# ====================================================================
# Config
# ====================================================================
MIMICRY_FILE  = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_i/swissprot/parsed/iedb_mhci_strong_mimicry.csv.gz"
PAIR_MAP_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/netmhcpan/mhc_i/swissprot/pair_id_mapping.csv.gz"
TCR_FILE      = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/mhc_i/iedb_epitopes_with_tcr_mhc_i.csv.gz"
OUTPUT_FILE   = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/embeddings/mhc_i/tcr_mapping/mimicry_with_tcr.csv.gz"

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)


# %% 
# ====================================================================
# Load
# ====================================================================
print("Loading strong mimicry candidates...")
mimicry = pd.read_csv(MIMICRY_FILE, low_memory=False)
print(f"  Records        : {len(mimicry):,}")
print(f"  Unique pair_id : {mimicry['pair_id'].nunique():,}")
print(f"  Unique viral peptides: {mimicry['viral_peptide'].nunique():,}")

print("\nLoading pair_id mapping...")
pair_map = pd.read_csv(PAIR_MAP_FILE, low_memory=False)
print(f"  Columns: {list(pair_map.columns)}")

print("\nLoading IEDB epitope-TCR data...")
tcr = pd.read_csv(TCR_FILE, low_memory=False)
print(f"  Records        : {len(tcr):,}")
print(f"  Unique structure_id: {tcr['structure_id'].nunique():,}")
print(f"  Columns: {list(tcr.columns)}")

# %%

# ====================================================================
# Step 1: Map pair_id → structure_id
# ====================================================================
print("\nMapping pair_id → structure_id...")

# pair_map has one structure_id per pair_id (possibly multiple rows
# for multi-position BLAST hits — take one per pair_id)
pair_to_struct = (
    pair_map[["pair_id", "structure_id"]]
    .drop_duplicates("pair_id")
)

mimicry = mimicry.merge(pair_to_struct, on="pair_id", how="left")
n_with_struct = mimicry["structure_id"].notna().sum()
print(f"  Mimicry records with structure_id: {n_with_struct:,} / {len(mimicry):,}")

# %%

# ====================================================================
# Step 2: Join to TCR data on structure_id
# ====================================================================
print("\nJoining to TCR data...")

# TCR columns to keep — adjust if your column names differ
tcr_cols = ["structure_id"]
for col in ["receptor_id", "receptor_group_iri", "t_cell_id",
            "alpha_cdr3", "beta_cdr3",
            "alpha_v_gene", "alpha_j_gene",
            "beta_v_gene",  "beta_j_gene",
            "alpha_full_seq", "beta_full_seq"]:
    if col in tcr.columns:
        tcr_cols.append(col)

tcr_slim = tcr[tcr_cols].drop_duplicates()

merged = mimicry.merge(tcr_slim, on="structure_id", how="left")

# %%

# ====================================================================
# Step 3: Coverage summary
# ====================================================================
print(f"\n{'=' * 55}")
print(f"  TCR coverage summary")
print(f"{'=' * 55}")

# At the pair level (one row per mimicry pair regardless of how many TCRs)
has_any_tcr    = merged["receptor_id"].notna() if "receptor_id" in merged.columns \
                 else merged["alpha_cdr3"].notna()

pairs_with_tcr = merged[has_any_tcr]["pair_id"].nunique()
pairs_total    = mimicry["pair_id"].nunique()
viral_with_tcr = merged[has_any_tcr]["viral_peptide"].nunique()
viral_total    = mimicry["viral_peptide"].nunique()

print(f"  Unique mimicry pairs total    : {pairs_total:,}")
print(f"  Pairs with ≥1 TCR             : {pairs_with_tcr:,}  "
      f"({100*pairs_with_tcr/pairs_total:.1f}%)")
print(f"  Pairs without TCR             : {pairs_total - pairs_with_tcr:,}  "
      f"({100*(pairs_total-pairs_with_tcr)/pairs_total:.1f}%)")
print(f"\n  Unique viral peptides total   : {viral_total:,}")
print(f"  Viral peptides with ≥1 TCR    : {viral_with_tcr:,}  "
      f"({100*viral_with_tcr/viral_total:.1f}%)")

# TCR chain completeness (for DecoderTCR — needs both chains)
if "alpha_cdr3" in merged.columns and "beta_cdr3" in merged.columns:
    has_both = (
        merged["alpha_cdr3"].notna() & merged["beta_cdr3"].notna()
    )
    pairs_both_chains = merged[has_both]["pair_id"].nunique()
    print(f"\n  Pairs with αβ CDR3 (both chains): {pairs_both_chains:,}  "
          f"({100*pairs_both_chains/pairs_total:.1f}%)")
    has_alpha_only = merged["alpha_cdr3"].notna() & merged["beta_cdr3"].isna()
    has_beta_only  = merged["alpha_cdr3"].isna()  & merged["beta_cdr3"].notna()
    print(f"  Pairs with α CDR3 only          : "
          f"{merged[has_alpha_only]['pair_id'].nunique():,}")
    print(f"  Pairs with β CDR3 only          : "
          f"{merged[has_beta_only]['pair_id'].nunique():,}")

# Unique receptor groups (deduplicated TCR clonotypes)
if "receptor_group_iri" in merged.columns:
    n_receptor_groups = merged.loc[has_any_tcr, "receptor_group_iri"].nunique()
    print(f"\n  Unique TCR clonotypes (receptor_group): {n_receptor_groups:,}")

# Per-allele breakdown
print(f"\n  TCR coverage by HLA allele (top 15):")
allele_coverage = (
    merged.groupby("mhc_allele")
    .apply(lambda g: pd.Series({
        "total_pairs"    : g["pair_id"].nunique(),
        "pairs_with_tcr" : g.loc[g["receptor_id"].notna()
                                  if "receptor_id" in g.columns
                                  else g["alpha_cdr3"].notna(),
                                  "pair_id"].nunique(),
    }))
    .assign(pct=lambda d: 100 * d["pairs_with_tcr"] / d["total_pairs"])
    .sort_values("total_pairs", ascending=False)
    .head(15)
)
print(allele_coverage.to_string())

# %%

# ====================================================================
# Save candidates with TCR
# ====================================================================
mimicry_with_tcr = merged[has_any_tcr].copy()
mimicry_with_tcr.to_csv(OUTPUT_FILE, index=False, compression="gzip")
print(f"\n  Saved {len(mimicry_with_tcr):,} records (mimicry pairs × TCRs) → {OUTPUT_FILE}")
print(f"  (One row per mimicry pair × TCR combination)")
print(f"\n✅ Done.")