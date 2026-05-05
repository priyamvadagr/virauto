#!/usr/bin/env python3
"""
======================================================================
Script: prepare_stitchr_input.py

Description:
    Prepare thimble (batch stitchr) input TSVs from strong mimicry
    candidates with paired TCR data.

    Produces two TSVs (alpha and beta chains run separately in thimble),
    plus a metadata table linking each TCR_name back to pair_id and
    receptor_group for downstream joining to DecoderTCR inputs.

Input:
    mimicry_with_tcr.csv.gz

Output:
    stitchr_input_alpha.tsv   → run with: thimble -i ... -r a -s HUMAN
    stitchr_input_beta.tsv    → run with: thimble -i ... -r b -s HUMAN
    stitchr_metadata.tsv      → TCR_name → pair_id, receptor_group, CDR3s

Usage:
    python prepare_stitchr_input.py

Then run thimble:
    thimble -i stitchr_input_alpha.tsv -o stitchr_output_alpha.tsv -r a -s HUMAN
    thimble -i stitchr_input_beta.tsv  -o stitchr_output_beta.tsv  -r b -s HUMAN
======================================================================
"""

import os
import re
import pandas as pd

# ====================================================================
# Config
# ====================================================================
MIMICRY_TCR_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/mimicry_candidates/mhc_i/tcr_mapping/mimicry_with_tcr.csv.gz"
OUT_DIR          = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/mimicry_candidates/mhc_i/stitchr"
os.makedirs(OUT_DIR, exist_ok=True)

ALPHA_TSV    = os.path.join(OUT_DIR, "stitchr_input_alpha.tsv")
BETA_TSV     = os.path.join(OUT_DIR, "stitchr_input_beta.tsv")
METADATA_TSV = os.path.join(OUT_DIR, "stitchr_metadata.tsv")

# ====================================================================
# Load
# ====================================================================
print("Loading mimicry TCR data...")
tcr = pd.read_csv(MIMICRY_TCR_FILE, low_memory=False)
print(f"  Total rows: {len(tcr):,}")

# ====================================================================
# Filter to rows with both CDR3 chains
# ====================================================================
both = tcr["alpha_cdr3"].notna() & tcr["beta_cdr3"].notna()
df = tcr[both].copy()
print(f"  Rows with both CDR3α and CDR3β: {len(df):,}")

# ====================================================================
# Normalize V/J gene names
# ====================================================================
EXPECTED_PREFIX = {
    "alpha_v_gene": "TRAV",
    "alpha_j_gene": "TRAJ",
    "beta_v_gene" : "TRBV",
    "beta_j_gene" : "TRBJ",
}

def normalize_vj_gene(name, expected_prefix):
    """
    Fix common IEDB formatting issues in V/J gene names:
      - Strip leading/trailing whitespace and non-breaking spaces (\xa0)
      - Fix wrong prefixes: TCRAV→TRAV, TCRBV→TRBV, TCRAJ→TRAJ, TCRBJ→TRBJ
      - Preserve allele suffixes (*01, *02) — stitchr uses these if present
      - Return None for genuinely unrecognizable names (e.g. TRDV in alpha col)
    """
    if pd.isna(name):
        return None
    name = str(name).strip()
    name = re.sub(r"^TCRA([VJ])", r"TRA\1", name)
    name = re.sub(r"^TCRB([VJ])", r"TRB\1", name)
    if not name.startswith(expected_prefix):
        return None
    return name

for col, pfx in EXPECTED_PREFIX.items():
    df[col + "_clean"] = df[col].apply(lambda x, p=pfx: normalize_vj_gene(x, p))

# ====================================================================
# Require all 4 V/J genes to be cleanly normalized
# ====================================================================
all_clean = (
    df["alpha_v_gene_clean"].notna() &
    df["alpha_j_gene_clean"].notna() &
    df["beta_v_gene_clean"].notna()  &
    df["beta_j_gene_clean"].notna()
)
df = df[all_clean].copy()
print(f"  Rows with all 4 V/J genes normalized: {len(df):,}")
print(f"  Unique pair_ids: {df['pair_id'].nunique():,}")

# ====================================================================
# Deduplicate on actual sequence content
# Unique unit = (alpha_cdr3, alpha_v, alpha_j, beta_cdr3, beta_v, beta_j)
# Many receptor_group IRIs can map to the same clonotype sequence
# ====================================================================
seq_cols = [
    "alpha_cdr3", "alpha_v_gene_clean", "alpha_j_gene_clean",
    "beta_cdr3",  "beta_v_gene_clean",  "beta_j_gene_clean",
]
df["tcr_seq_key"] = (
    df["alpha_v_gene_clean"] + "|" + df["alpha_cdr3"] + "|" + df["alpha_j_gene_clean"] + "||" +
    df["beta_v_gene_clean"]  + "|" + df["beta_cdr3"]  + "|" + df["beta_j_gene_clean"]
)

# Keep one row per (pair_id, unique TCR sequence)
df_dedup = df.drop_duplicates(subset=["pair_id", "tcr_seq_key"]).copy()
print(f"\n  After deduplication on sequence content:")
print(f"    Unique (pair_id, TCR sequence) jobs : {len(df_dedup):,}")
print(f"    Unique TCR sequences                : {df_dedup['tcr_seq_key'].nunique():,}")
print(f"    Unique pair_ids                     : {df_dedup['pair_id'].nunique():,}")

# ====================================================================
# Assign stable TCR_name
# Format: TCR_{zero-padded index}
# This is the join key across alpha TSV, beta TSV, metadata, and
# eventually DecoderTCR input
# ====================================================================
# One TCR_name per unique TCR sequence (not per pair — same TCR can
# appear across multiple pairs; we score it once per unique sequence,
# then fan out to all pairs in the metadata table)
unique_tcrs = (
    df_dedup[["tcr_seq_key"] + seq_cols]
    .drop_duplicates("tcr_seq_key")
    .reset_index(drop=True)
)
unique_tcrs["TCR_name"] = ["TCR_{:05d}".format(i) for i in range(len(unique_tcrs))]
print(f"\n  Unique TCR sequences to stitch: {len(unique_tcrs):,}")

# Map TCR_name back to df_dedup
df_dedup = df_dedup.merge(unique_tcrs[["tcr_seq_key", "TCR_name"]], on="tcr_seq_key", how="left")

# ====================================================================
# Build thimble input TSVs
# Thimble mandatory columns:
#   Alpha: TCR_name, TRAV, TRAJ, TRA_CDR3
#   Beta:  TCR_name, TRBV, TRBJ, TRB_CDR3
# ====================================================================
# Alpha — full template columns, optional ones left empty
alpha_input = unique_tcrs[["TCR_name", "alpha_v_gene_clean", "alpha_j_gene_clean", "alpha_cdr3"]].rename(columns={
    "alpha_v_gene_clean": "TRAV",
    "alpha_j_gene_clean": "TRAJ",
    "alpha_cdr3"        : "TRA_CDR3",
})
for col in ["TRBV", "TRBJ", "TRB_CDR3", "TRAC", "TRBC",
            "TRA_leader", "TRB_leader", "Linker", "Link_order",
            "TRA_5_prime_seq", "TRA_3_prime_seq",
            "TRB_5_prime_seq", "TRB_3_prime_seq"]:
    alpha_input[col] = ""

# Beta — same full template
beta_input = unique_tcrs[["TCR_name", "beta_v_gene_clean", "beta_j_gene_clean", "beta_cdr3"]].rename(columns={
    "beta_v_gene_clean": "TRBV",
    "beta_j_gene_clean": "TRBJ",
    "beta_cdr3"        : "TRB_CDR3",
})
for col in ["TRAV", "TRAJ", "TRA_CDR3", "TRAC", "TRBC",
            "TRA_leader", "TRB_leader", "Linker", "Link_order",
            "TRA_5_prime_seq", "TRA_3_prime_seq",
            "TRB_5_prime_seq", "TRB_3_prime_seq"]:
    beta_input[col] = ""

# Reorder to match exact template column order
template_cols = ["TCR_name", "TRAV", "TRAJ", "TRA_CDR3", "TRBV", "TRBJ", "TRB_CDR3",
                 "TRAC", "TRBC", "TRA_leader", "TRB_leader", "Linker", "Link_order",
                 "TRA_5_prime_seq", "TRA_3_prime_seq", "TRB_5_prime_seq", "TRB_3_prime_seq"]

alpha_input[template_cols].to_csv(ALPHA_TSV, sep="\t", index=False)
beta_input[template_cols].to_csv(BETA_TSV,   sep="\t", index=False)
print(f"\n  Alpha thimble input → {ALPHA_TSV}  ({len(alpha_input):,} rows)")
print(f"  Beta  thimble input → {BETA_TSV}  ({len(beta_input):,} rows)")

# ====================================================================
# Metadata table: TCR_name → pair_id + all context needed downstream
# One row per (TCR_name, pair_id) — a TCR can appear across multiple
# pairs, and a pair can have multiple TCRs
# ====================================================================
meta_cols = [
    "TCR_name", "pair_id",
    "viral_peptide", "human_peptide", "mhc_allele",
    "receptor_group_iri",
    "alpha_cdr3", "beta_cdr3",
    "alpha_v_gene_clean", "alpha_j_gene_clean",
    "beta_v_gene_clean",  "beta_j_gene_clean",
]
meta_cols = [c for c in meta_cols if c in df_dedup.columns]
metadata = df_dedup[meta_cols].drop_duplicates().reset_index(drop=True)
metadata.to_csv(METADATA_TSV, sep="\t", index=False)
print(f"  Metadata table       → {METADATA_TSV}  ({len(metadata):,} rows)")

# ====================================================================
# Summary + thimble run commands
# ====================================================================
print(f"\n{'=' * 60}")
print("  Next steps: run thimble on the cluster")
print(f"{'=' * 60}")
print(f"""
  # Initialize stitchr data (first time only)
  stitchrdl -s HUMAN

  # Run alpha chain reconstruction
  thimble \\
      -in {ALPHA_TSV} \\
      -o {os.path.join(OUT_DIR, 'stitchr_output_alpha.tsv')} \\
      -r a -s HUMAN

  # Run beta chain reconstruction
  thimble \\
      -in {BETA_TSV} \\
      -o {os.path.join(OUT_DIR, 'stitchr_output_beta.tsv')} \\
      -r b -s HUMAN
""")
print("✅ Done.")
