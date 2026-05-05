#!/usr/bin/env python3
"""
compute_mimetopes_per_protein.py

Count unique mimetope positions (Ni) per human protein per peptide length k,
then merge with precomputed Li (k-mer window counts) from the proteome.

Input:
    - iedb_mhci_strong_mimicry.csv.gz
    - pair_id_mapping.csv.gz (with hu_prot_id, sstart, send, human_mimic_sequence)
    - uniprot_human_kmer_counts_filtered.tsv (Li per protein, precomputed)

Output:
    - mimetopes_per_protein.tsv
      Columns: hu_prot_id, entry_name, description, length,
               Li_8..Li_k, Ni_8..Ni_k

Usage:
    python compute_mimetopes_per_protein.py \
        --strong-mimicry iedb_mhci_strong_mimicry.csv.gz \
        --pair-map pair_id_mapping.csv.gz \
        --kmer-counts uniprot_human_kmer_counts_filtered.tsv \
        --outfile mimetopes_per_protein.tsv
"""

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Compute Ni (mimetopes per protein per k) and merge with Li"
    )
    parser.add_argument("--strong-mimicry", required=True,
                        help="iedb_mhci_strong_mimicry.csv.gz")
    parser.add_argument("--pair-map", required=True,
                        help="pair_id_mapping.csv.gz")
    parser.add_argument("--kmer-counts", required=True,
                        help="uniprot_human_kmer_counts_filtered.tsv (Li per protein)")
    parser.add_argument("--outfile", required=True,
                        help="Output TSV with Li and Ni per protein")
    args = parser.parse_args()

    # --- Load strong mimicry ---
    print("Loading strong mimicry...")
    strong = pd.read_csv(args.strong_mimicry, low_memory=False)
    strong_pair_ids = set(strong["pair_id"].unique())
    print(f"  {len(strong):,} rows, {len(strong_pair_ids):,} unique pair_ids")

    # --- Load pair_id_mapping ---
    print("\nLoading pair ID mapping...")
    pair_map = pd.read_csv(args.pair_map, low_memory=False)
    print(f"  {len(pair_map):,} rows")
    print(f"  Columns: {pair_map.columns.tolist()}")

    # Filter to strong pairs
    mimetopes = pair_map[pair_map["pair_id"].isin(strong_pair_ids)].copy()
    print(f"  {len(mimetopes):,} rows matching strong pair_ids")

    # Compute peptide length
    mimetopes["pep_len"] = mimetopes["human_mimic_sequence"].str.len()

    # Deduplicate: unique (hu_prot_id, sstart, send) positions
    mimetopes = mimetopes.drop_duplicates(subset=["hu_prot_id", "sstart", "send"])
    print(f"  {len(mimetopes):,} unique mimetope positions after dedup")

    k_values = sorted(mimetopes["pep_len"].unique())
    print(f"  Peptide lengths: {k_values}")

    # --- Count Ni per (protein, k) ---
    print("\nCounting Ni per protein per k...")
    ni_counts = (
        mimetopes
        .groupby(["hu_prot_id", "pep_len"])
        .size()
        .reset_index(name="Ni")
    )

    # Pivot to wide format: one row per protein, Ni_8, Ni_9, etc.
    ni_wide = ni_counts.pivot(
        index="hu_prot_id", columns="pep_len", values="Ni"
    ).fillna(0).astype(int)
    ni_wide.columns = [f"Ni_{k}" for k in ni_wide.columns]
    ni_wide = ni_wide.reset_index()

    print(f"  {len(ni_wide):,} proteins with at least one mimetope")
    for col in ni_wide.columns:
        if col.startswith("Ni_"):
            n_pos = (ni_wide[col] > 0).sum()
            print(f"    {col}: {n_pos:,} proteins with Ni > 0")

    # --- Load Li (precomputed k-mer counts) ---
    print(f"\nLoading precomputed Li: {args.kmer_counts}")
    li_df = pd.read_csv(args.kmer_counts, sep="\t")
    print(f"  {len(li_df):,} proteins")
    print(f"  Columns: {li_df.columns.tolist()}")

    # Rename Li columns to standard format
    # Input has: 8mer_count, 9mer_count, etc.
    li_rename = {}
    for col in li_df.columns:
        if col.endswith("mer_count"):
            k = int(col.replace("mer_count", ""))
            li_rename[col] = f"Li_{k}"
    li_df = li_df.rename(columns=li_rename)
    print(f"  Renamed: {li_rename}")

    # --- Merge Li and Ni ---
    print("\nMerging Li and Ni...")
    merged = li_df.merge(ni_wide, on="hu_prot_id", how="left")

    # Fill missing Ni with 0 (proteins with no mimetopes)
    ni_cols = [c for c in merged.columns if c.startswith("Ni_")]
    for col in ni_cols:
        merged[col] = merged[col].fillna(0).astype(int)

    # Ensure we have matching Li and Ni columns
    li_cols = sorted([c for c in merged.columns if c.startswith("Li_")])
    ni_cols = sorted([c for c in merged.columns if c.startswith("Ni_")])
    print(f"  Li columns: {li_cols}")
    print(f"  Ni columns: {ni_cols}")

    n_with_mimetopes = (merged[ni_cols].sum(axis=1) > 0).sum()
    print(f"  Total proteins: {len(merged):,}")
    print(f"  Proteins with ≥1 mimetope: {n_with_mimetopes:,}")

    # Sanity check: Ni should never exceed Li for any protein
    for li_col, ni_col in zip(li_cols, ni_cols):
        k_li = int(li_col.split("_")[1])
        k_ni = int(ni_col.split("_")[1])
        if k_li == k_ni:
            violations = merged[merged[ni_col] > merged[li_col]]
            if len(violations) > 0:
                print(f"  ⚠️ k={k_li}: {len(violations)} proteins where Ni > Li!")

    # --- Save ---
    merged.to_csv(args.outfile, sep="\t", index=False)
    print(f"\n  Output: {args.outfile}")
    print(f"  {len(merged):,} rows")

    print(f"\n✅ Done.")


if __name__ == "__main__":
    main()