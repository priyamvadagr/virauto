#!/usr/bin/env python3
"""
===============================================================================
Script: summarize_pair_binding_classes.py
Author: Priyamvada Guha Roy
===============================================================================
Description:
    Summarize viral–human peptide pair binding classes across all HLA alleles.

    For each unique (viral, human) peptide pair:
      - Count number of alleles tested
      - Count how many alleles show Human_dominant, Viral_dominant, Equivalent_binding, or Non_binder
      - Compute fractions of each class
      - Generate both global (all alleles) and HLA-A–specific summaries

Inputs:
    Parquet dataset partitioned by allele
Outputs:
    pair_summary_all.csv / pair_summary_HLA_A.csv
===============================================================================
"""

import os
import pandas as pd
import pyarrow.parquet as pq
from tqdm import tqdm
import re

# =============================================================================
# 1. Configuration
# =============================================================================
BASE_PATH = "/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/virscan/9_mers/Class_I/HLA-A/"
DATASET_PATH = os.path.join(BASE_PATH, "data")
OUT_DIR = os.path.join(BASE_PATH, "summaries")
os.makedirs(OUT_DIR, exist_ok=True)

# =============================================================================
# 2. Helper function
# =============================================================================
def summarize_pairwise_binding(df: pd.DataFrame, label: str):
    """Summarize number of alleles per binding class per peptide pair."""
    print(f"[INFO] Summarizing {label} data ({len(df):,} rows, {df['allele'].nunique()} alleles)")

    # --- Step 1: Get unique mapping of pair_id → UniProt IDs ---
    pair_to_proteins = (
        df.groupby("pair_id")
          .agg({
              "uniprot_viral": lambda x: ";".join(sorted(set(map(str, x)))),
              "uniprot_human": lambda x: ";".join(sorted(set(map(str, x)))),
              "peptide_viral": "first",
              "peptide_human": "first"
          })
          .reset_index()
    )

    # --- Step 2: Count number of alleles per binding class ---
    summary = (
        df.groupby(["pair_id", "Binding_Class"])
          .agg(n_alleles=("allele", "nunique"))
          .reset_index()
    )

    # --- Step 3: Pivot binding class counts ---
    summary_pivot = (
        summary.pivot(index="pair_id", columns="Binding_Class", values="n_alleles")
        .fillna(0)
        .reset_index()
    )

    # --- Step 4: Merge UniProt / peptide metadata back ---
    summary_pivot = summary_pivot.merge(pair_to_proteins, on="pair_id", how="left")

    # --- Step 5: Compute totals and fractions ---
    class_cols = ["Human_dominant", "Viral_dominant", "Equivalent_binding"]
    summary_pivot["n_alleles_tested"] = summary_pivot[class_cols].sum(axis=1)

    for col in class_cols:
        summary_pivot[f"{col}_fraction"] = (
            summary_pivot[col] / summary_pivot["n_alleles_tested"]
        ).fillna(0)

    # --- Step 6: Sort by strongest human bias ---
    summary_pivot.sort_values("Human_dominant_fraction", ascending=False, inplace=True)

    # --- Step 7: Save output ---
    out_csv = os.path.join(OUT_DIR, f"pair_binding_summary_{label}.csv")
    summary_pivot.to_csv(out_csv, index=False)
    print(f"[INFO] Saved {label} summary → {out_csv}")
    return summary_pivot


# =============================================================================
# 3. Load Parquet dataset
# =============================================================================
print(f"[INFO] Reading Parquet dataset from {DATASET_PATH}")
dataset = pq.ParquetDataset(DATASET_PATH)
df = dataset.read_pandas().to_pandas()

# Keep only relevant columns
keep_cols = ["peptide_viral", "uniprot_viral", "uniprot_human", "pair_id", "peptide_human", "allele", "Binding_Class"]
df = df[keep_cols].dropna(subset=["Binding_Class"])

# Normalize allele names and extract HLA type
df["HLA_class"] = df["allele"].apply(
    lambda x: (m.group(1)
               if (m := re.search(r"HLA[-_: ]?([A-G])", str(x), flags=re.IGNORECASE))
               else "Unknown")
)

print(f"[INFO] Loaded {len(df):,} rows from {df['allele'].nunique()} alleles.")
print(f"[INFO] HLA class distribution:\n{df['HLA_class'].value_counts()}")

# =============================================================================
# 4. Global summary (all alleles)
# =============================================================================
summary_all = summarize_pairwise_binding(df, "all")

# =============================================================================
# 5. HLA-A specific summary
# =============================================================================
df_hla_a = df[df["HLA_class"] == "A"].copy()
if len(df_hla_a) > 0:
    summary_hla_a = summarize_pairwise_binding(df_hla_a, "HLA_A")
else:
    print("[WARN] No HLA-A alleles found in dataset — skipping HLA-A summary.")
