#!/usr/bin/env python3
"""
======================================================================
Script: merge_iedb_epitopes_with_tcr.py
Description:
    Merge IEDB T cell epitope data with TCR sequences from the IEDB 
    tcr_full_v3 bulk export. Saves clean, minimal-column gzipped CSVs.

Inputs:
    - IEDB epitope CSVs (from query_iedb_epitopes.py)
    - tcr_full_v3.csv (IEDB bulk export)

Outputs (gzipped):
    - iedb_epitopes_with_tcr_{mhc_class}.csv.gz
    - iedb_epitopes_no_tcr_{mhc_class}.csv.gz
    - tcr_full_v3_parsed.csv.gz

Dependencies:
    pandas

Usage:
    python merge_iedb_epitopes_with_tcr.py
======================================================================
"""

import pandas as pd
import os
import ast
from collections import OrderedDict

# ====================================================
# Config
# ====================================================
IEDB_DIR = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb"
TCR_FILE = os.path.join(IEDB_DIR, "all_receptor_data", "tcr_full_v3.csv")

MHC_CLASSES = {
    "mhc_i": "MHC Class I",
    "mhc_ii": "MHC Class II",
}


# ====================================================
# Step 1: Load and parse TCR bulk export
# ====================================================
def load_tcr_export(filepath):
    """Load tcr_full_v3.csv with its two-row header."""
    print("=" * 60)
    print("Loading TCR bulk export")
    print("=" * 60)

    # Build unique column names from two-row header
    headers = pd.read_csv(filepath, nrows=0, header=[0, 1])
    col_names = []
    for cat, field in headers.columns:
        clean_name = f"{cat.strip()} - {field.strip()}"
        if clean_name in col_names:
            i = 2
            while f"{clean_name}_{i}" in col_names:
                i += 1
            clean_name = f"{clean_name}_{i}"
        col_names.append(clean_name)

    df = pd.read_csv(filepath, skiprows=2, header=None, names=col_names,
                     low_memory=False, on_bad_lines="skip")
    print(f"  Raw TCR records: {len(df):,}")
    return df


def parse_tcr_export(tcr_raw):
    """
    Extract only the columns we need from the TCR export and rename 
    them to clean, short names.
    """
    print("  Parsing and trimming columns...")

    # Direct column name mapping — no pattern matching, just exact names
    rename = {
        "Receptor - Group IRI": "receptor_group_iri",
        "Receptor - IEDB Receptor ID": "receptor_id",
        "Receptor - Reference Name": "receptor_name",
        "Receptor - Type": "receptor_type",
        "Epitope - IEDB IRI": "epitope_iri",
        "Epitope - Name": "epitope_sequence",
        "Epitope - Source Organism": "epitope_source_organism",
        "Assay - MHC Allele Names": "tcr_mhc_allele",
        # Chain 1 (alpha)
        "Chain 1 - Type": "chain1_type",
        "Chain 1 - Curated V Gene": "alpha_v_gene",
        "Chain 1 - Curated J Gene": "alpha_j_gene",
        "Chain 1 - Protein Sequence": "alpha_full_seq",
        "Chain 1 - CDR3 Curated": "alpha_cdr3_curated",
        "Chain 1 - CDR3 Calculated": "alpha_cdr3_calc",
        "Chain 1 - CDR1 Curated": "alpha_cdr1_curated",
        "Chain 1 - CDR1 Calculated": "alpha_cdr1_calc",
        "Chain 1 - CDR2 Curated": "alpha_cdr2_curated",
        "Chain 1 - CDR2 Calculated": "alpha_cdr2_calc",
        # Chain 2 (beta)
        "Chain 2 - Type": "chain2_type",
        "Chain 2 - Curated V Gene": "beta_v_gene",
        "Chain 2 - Curated J Gene": "beta_j_gene",
        "Chain 2 - Protein Sequence": "beta_full_seq",
        "Chain 2 - CDR3 Curated": "beta_cdr3_curated",
        "Chain 2 - CDR3 Calculated": "beta_cdr3_calc",
        "Chain 2 - CDR1 Curated": "beta_cdr1_curated",
        "Chain 2 - CDR1 Calculated": "beta_cdr1_calc",
        "Chain 2 - CDR2 Curated": "beta_cdr2_curated",
        "Chain 2 - CDR2 Calculated": "beta_cdr2_calc",
    }

    # Keep only columns that exist in the data
    existing = {k: v for k, v in rename.items() if k in tcr_raw.columns}
    missing = [k for k in rename if k not in tcr_raw.columns]
    if missing:
        print(f"  ⚠️ Missing columns (skipped): {missing}")

    df = tcr_raw[list(existing.keys())].rename(columns=existing)

    # Convert receptor_id to numeric
    df["receptor_id"] = pd.to_numeric(df["receptor_id"], errors="coerce")

    # Extract structure_id from epitope_iri
    if "epitope_iri" in df.columns:
        df["structure_id"] = (
            df["epitope_iri"].str.extract(r"/epitope/(\d+)")[0].astype(float)
        )
        df = df.drop(columns=["epitope_iri"])

    # Consolidate CDR columns: prefer curated, fall back to calculated
    for chain, prefix in [("alpha", "alpha"), ("beta", "beta")]:
        for cdr in ["cdr1", "cdr2", "cdr3"]:
            curated = f"{prefix}_{cdr}_curated"
            calc = f"{prefix}_{cdr}_calc"
            final = f"{prefix}_{cdr}"
            if curated in df.columns and calc in df.columns:
                df[final] = df[curated].fillna(df[calc])
                df = df.drop(columns=[curated, calc])
            elif curated in df.columns:
                df = df.rename(columns={curated: final})
            elif calc in df.columns:
                df = df.rename(columns={calc: final})

    # Summary
    print(f"  Parsed TCR records: {len(df):,}")
    print(f"  Unique receptor IDs: {df['receptor_id'].nunique():,}")
    print(f"  Unique epitope sequences: {df['epitope_sequence'].nunique():,}")
    print(f"  With alpha CDR3: {df.get('alpha_cdr3', pd.Series()).notna().sum():,}")
    print(f"  With beta CDR3: {df.get('beta_cdr3', pd.Series()).notna().sum():,}")
    if "alpha_cdr3" in df.columns and "beta_cdr3" in df.columns:
        paired = (df["alpha_cdr3"].notna() & df["beta_cdr3"].notna()).sum()
        print(f"  With paired αβ CDR3: {paired:,}")

    # Save parsed TCR data
    out_path = os.path.join(IEDB_DIR, "all_receptor_data", "tcr_full_v3_parsed.csv.gz")
    df.to_csv(out_path, index=False, compression="gzip")
    print(f"\n  ✅ Parsed TCR data → {out_path}")

    return df


# ====================================================
# Step 2: Load epitope data
# ====================================================
def load_epitope_data(mhc_class_key):
    """Load epitope data, keep only relevant columns."""
    filepath = os.path.join(IEDB_DIR, mhc_class_key,
                            f"iedb_tcell_{mhc_class_key}_epitopes.csv")
    if not os.path.exists(filepath):
        print(f"  ⚠️ Not found: {filepath}")
        return pd.DataFrame()

    df = pd.read_csv(filepath, low_memory=False)
    print(f"  Loaded {len(df):,} epitope records")

    # Keep only what we need
    keep_cols = [
        "structure_id",
        "linear_sequence",
        "peptide_length",
        "parent_source_antigen_name",
        "parent_source_antigen_source_org_name",
        "source_organism_name",
        "mhc_allele_name",
        "mhc_class",
        "qualitative_measure",
        "assay_names",
        "pubmed_id",
        "receptor_ids",
        "tcell_id",
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    df = df[keep_cols]

    # Parse receptor_ids
    def parse_receptor_ids(val):
        if pd.isna(val):
            return []
        if isinstance(val, list):
            return val
        try:
            parsed = ast.literal_eval(str(val))
            return parsed if isinstance(parsed, list) else [parsed]
        except:
            return []

    df["receptor_ids_list"] = df["receptor_ids"].apply(parse_receptor_ids)
    df["has_receptor"] = df["receptor_ids_list"].apply(len) > 0

    return df


# ====================================================
# Step 3: Merge
# ====================================================
def merge_data(epitope_df, tcr_df, mhc_class_key):
    """Merge epitopes with TCR data via structure_id and receptor_id."""
    label = MHC_CLASSES[mhc_class_key]
    print(f"\n{'=' * 60}")
    print(f"Merging: {label}")
    print(f"{'=' * 60}")

    if epitope_df.empty or tcr_df.empty:
        print("  ⚠️ Empty input.")
        return pd.DataFrame(), epitope_df

    # --- Approach 1: Join on structure_id ---
    merged = pd.DataFrame()
    if "structure_id" in epitope_df.columns and "structure_id" in tcr_df.columns:
        merged = epitope_df.merge(
            tcr_df, on="structure_id", how="inner", suffixes=("", "_tcr")
        )
        print(f"  Structure_id join: {len(merged):,} records, "
              f"{merged['linear_sequence'].nunique():,} unique epitopes")

    # --- Approach 2: Join on receptor_id (for records not caught above) ---
    if epitope_df["has_receptor"].any():
        epi_with_receptor = epitope_df[epitope_df["has_receptor"]].copy()
        epi_exploded = epi_with_receptor.explode("receptor_ids_list")
        epi_exploded["receptor_id_match"] = pd.to_numeric(
            epi_exploded["receptor_ids_list"], errors="coerce"
        )

        merged2 = epi_exploded.merge(
            tcr_df, left_on="receptor_id_match", right_on="receptor_id",
            how="inner", suffixes=("", "_tcr")
        )
        print(f"  Receptor_id join: {len(merged2):,} records, "
              f"{merged2['linear_sequence'].nunique():,} unique epitopes")

        if merged.empty:
            merged = merged2
        elif not merged2.empty:
            # Add rows from approach 2 not already in approach 1
            merged = pd.concat([merged, merged2]).drop_duplicates(
                subset=["tcell_id", "receptor_id"], keep="first"
            )
            print(f"  Combined: {len(merged):,} records")

    # --- Clean up merged output ---
    if not merged.empty:
        # Drop helper columns
        drop_cols = ["receptor_ids", "receptor_ids_list", "has_receptor",
                     "receptor_id_match", "receptor_ids_list"]
        merged = merged.drop(columns=[c for c in drop_cols if c in merged.columns],
                             errors="ignore")

        # Remove duplicate _tcr columns where values are redundant
        tcr_suffix_cols = [c for c in merged.columns if c.endswith("_tcr")]
        merged = merged.drop(columns=tcr_suffix_cols, errors="ignore")

        # --- Flag DecoderTCR-ready pairs ---
        # DecoderTCR requires: paired alpha+beta CDR3, V and J genes 
        # for both chains (for Stitchr reconstruction), peptide, MHC allele
        has_alpha_cdr3 = merged.get("alpha_cdr3", pd.Series(dtype=str)).notna()
        has_beta_cdr3 = merged.get("beta_cdr3", pd.Series(dtype=str)).notna()
        has_alpha_v = merged.get("alpha_v_gene", pd.Series(dtype=str)).notna()
        has_alpha_j = merged.get("alpha_j_gene", pd.Series(dtype=str)).notna()
        has_beta_v = merged.get("beta_v_gene", pd.Series(dtype=str)).notna()
        has_beta_j = merged.get("beta_j_gene", pd.Series(dtype=str)).notna()
        has_mhc = merged.get("mhc_allele_name", pd.Series(dtype=str)).notna()
        has_peptide = merged.get("linear_sequence", pd.Series(dtype=str)).notna()

        # Minimum: paired CDR3 + peptide + MHC
        merged["decoderTCR_paired_cdr3"] = (
            has_alpha_cdr3 & has_beta_cdr3 & has_peptide & has_mhc
        )
        # Full: paired CDR3 + V/J genes for both chains (Stitchr-ready) + peptide + MHC
        merged["decoderTCR_full_ready"] = (
            has_alpha_cdr3 & has_beta_cdr3 &
            has_alpha_v & has_alpha_j &
            has_beta_v & has_beta_j &
            has_peptide & has_mhc
        )

        n_paired = merged["decoderTCR_paired_cdr3"].sum()
        n_full = merged["decoderTCR_full_ready"].sum()
        print(f"\n  DecoderTCR readiness:")
        print(f"    Paired αβ CDR3 + peptide + MHC: {n_paired:,}")
        print(f"    Full (+ V/J genes for Stitchr): {n_full:,}")

    # --- Epitopes without TCR ---
    if merged.empty:
        no_tcr = epitope_df.copy()
    else:
        matched_structures = set(merged["structure_id"].dropna().unique())
        no_tcr = epitope_df[~epitope_df["structure_id"].isin(matched_structures)]

    no_tcr = no_tcr.drop(columns=["receptor_ids", "receptor_ids_list", "has_receptor"],
                         errors="ignore")
    no_tcr = no_tcr.drop_duplicates(subset=["linear_sequence", "mhc_allele_name"])

    # --- Save ---
    out_dir = os.path.join(IEDB_DIR, mhc_class_key)

    if not merged.empty:
        merged_path = os.path.join(out_dir, f"iedb_epitopes_with_tcr_{mhc_class_key}.csv.gz")
        merged.to_csv(merged_path, index=False, compression="gzip")
        print(f"\n  ✅ Epitopes with TCR: {len(merged):,} → {merged_path}")

    no_tcr_path = os.path.join(out_dir, f"iedb_epitopes_no_tcr_{mhc_class_key}.csv.gz")
    no_tcr.to_csv(no_tcr_path, index=False, compression="gzip")
    print(f"  ✅ Epitopes without TCR: {len(no_tcr):,} → {no_tcr_path}")

    return merged, no_tcr

# ====================================================
# STEP 4: Write FASTA files for BLAST
# ====================================================
def write_epitope_fasta(results):
    """
    Write unique epitope sequences to FASTA files for BLAST against
    the human proteome. Produces:
      - Per MHC class FASTA files
      - Combined all-epitopes FASTA
    
    FASTA header format:
      >structure_id|sequence|organism|mhc_allele
    """
    print(f"\n{'=' * 60}")
    print("Writing FASTA files for BLAST")
    print(f"{'=' * 60}")

    fasta_dir = os.path.join(IEDB_DIR, "fasta")
    os.makedirs(fasta_dir, exist_ok=True)

    all_seqs = OrderedDict()  # sequence → header info, preserves order, deduplicates

    for mhc_key in ["mhc_i", "mhc_ii"]:
        label = MHC_CLASSES[mhc_key]
        class_seqs = OrderedDict()

        # Collect from both merged and no_tcr
        for source_name, source_df in [("merged", results[mhc_key]["merged"]),
                                        ("no_tcr", results[mhc_key]["no_tcr"])]:
            if source_df.empty:
                continue
            for _, row in source_df.iterrows():
                seq = row.get("linear_sequence", "")
                if pd.isna(seq) or not seq:
                    continue
                seq = str(seq).strip().upper()

                if seq not in class_seqs:
                    sid = row.get("structure_id", "NA")
                    org = row.get("parent_source_antigen_source_org_name", "NA")
                    mhc = row.get("mhc_allele_name", "NA")
                    # Clean special characters from header
                    org = str(org).replace(",", ";").replace(" ", "_") if not pd.isna(org) else "NA"
                    mhc = str(mhc).replace(",", ";").replace(" ", "_") if not pd.isna(mhc) else "NA"
                    sid = int(sid) if not pd.isna(sid) else "NA"
                    class_seqs[seq] = f"{sid}|{seq}|{org}|{mhc}"

        # Write per-class FASTA
        fasta_path = os.path.join(fasta_dir, f"iedb_epitopes_{mhc_key}.fasta")
        with open(fasta_path, "w") as f:
            for seq, header in class_seqs.items():
                f.write(f">{header}\n{seq}\n")

        print(f"  {label}: {len(class_seqs):,} unique sequences → {fasta_path}")
        all_seqs.update(class_seqs)

    # Write combined FASTA
    combined_path = os.path.join(fasta_dir, "iedb_epitopes_all.fasta")
    with open(combined_path, "w") as f:
        for seq, header in all_seqs.items():
            f.write(f">{header}\n{seq}\n")

    print(f"  Combined: {len(all_seqs):,} unique sequences → {combined_path}")

    return all_seqs


# ====================================================
# Summary
# ====================================================
def print_summary(merged_df, no_tcr_df, mhc_class_key):
    label = MHC_CLASSES[mhc_class_key]
    print(f"\n{'=' * 60}")
    print(f"Summary: {label}")
    print(f"{'=' * 60}")

    if not merged_df.empty:
        print(f"  Epitope-TCR records: {len(merged_df):,}")
        print(f"  Unique epitopes with TCR: {merged_df['linear_sequence'].nunique():,}")
        print(f"  Unique TCRs: {merged_df['receptor_id'].nunique():,}")

        for col, name in [("alpha_cdr3", "alpha CDR3"),
                          ("beta_cdr3", "beta CDR3")]:
            if col in merged_df.columns:
                print(f"  With {name}: {merged_df[col].notna().sum():,}")

        if "alpha_cdr3" in merged_df.columns and "beta_cdr3" in merged_df.columns:
            paired = (merged_df["alpha_cdr3"].notna() &
                      merged_df["beta_cdr3"].notna()).sum()
            print(f"  With paired αβ: {paired:,}")

        org_col = "parent_source_antigen_source_org_name"
        if org_col in merged_df.columns:
            print(f"\n  Top 10 organisms (with TCR):")
            for org, n in merged_df[org_col].value_counts().head(10).items():
                print(f"    {org}: {n}")

        if "mhc_allele_name" in merged_df.columns:
            print(f"\n  Top 10 MHC alleles (with TCR):")
            for allele, n in merged_df["mhc_allele_name"].value_counts().head(10).items():
                print(f"    {allele}: {n}")

        print(f"\n  Columns in output:")
        for col in merged_df.columns:
            print(f"    {col}")
    else:
        print("  No epitope-TCR pairs found.")

    print(f"\n  Epitopes without TCR (for VDJdb): "
          f"{no_tcr_df['linear_sequence'].nunique():,} unique sequences")


# ====================================================
# Main
# ====================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Merge IEDB Epitopes with TCR Bulk Export")
    print("=" * 60)

    # Load and parse TCR data once
    tcr_raw = load_tcr_export(TCR_FILE)
    tcr_df = parse_tcr_export(tcr_raw)

    # Process each MHC class
    results = {}
    for mhc_key in ["mhc_i", "mhc_ii"]:
        print(f"\n\n{'#' * 60}")
        print(f"# {MHC_CLASSES[mhc_key]}")
        print(f"{'#' * 60}")

        epitope_df = load_epitope_data(mhc_key)
        merged_df, no_tcr_df = merge_data(epitope_df, tcr_df, mhc_key)
        results[mhc_key] = {"merged": merged_df, "no_tcr": no_tcr_df}

    # Summaries
    print(f"\n\n{'#' * 60}")
    print("# FINAL SUMMARIES")
    print(f"{'#' * 60}")

    for mhc_key in ["mhc_i", "mhc_ii"]:
        r = results[mhc_key]
        print_summary(r["merged"], r["no_tcr"], mhc_key)

    # Combined
    print(f"\n{'=' * 60}")
    print("Combined")
    print(f"{'=' * 60}")
    total_with = sum(len(r["merged"]) for r in results.values()
                     if not r["merged"].empty)
    total_without = sum(r["no_tcr"]["linear_sequence"].nunique()
                        for r in results.values() if not r["no_tcr"].empty)
    print(f"  Epitope-TCR records: {total_with:,}")
    print(f"  Epitopes without TCR (for VDJdb): {total_without:,}")

    # Write FASTA files
    write_epitope_fasta(results)

    print(f"\n✅ Done.")


