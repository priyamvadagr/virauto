#!/usr/bin/env python3
"""
======================================================================
Script: parse_netmhcpan_results.py
Description:
    Unified parser for NetMHCpan (Class I) and NetMHCIIpan (Class II)
    xls output files for IEDB mimicry pairs.

    Handles column naming differences between Class I and II:
      Class I:  core, icore, BA_score, BA_Rank
      Class II: Core, Inverted, Score_BA, nM, Rank_BA

    ΔBA categories:
      Viral-dominant:    viral binds & ΔBA_score ≥ +threshold
      Human-dominant:    viral binds & ΔBA_score ≤ -threshold
      Equivalent:        viral binds & |ΔBA_score| < threshold
      Non-binder:        viral Rank > binder threshold

Input:
    - NetMHCpan/NetMHCIIpan xls files in results directory
    - Original pairs file for metadata
    - Pair ID mapping

Output:
    - Merged results with ΔBA scores
    - Categorized mimicry candidates

Usage:
    # Class I
    python parse_netmhcpan_results.py \
        --mhc-class I \
        --result-dir results/netmhcpan/iedb/mhc_i/ \
        --pairs-file data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla.csv.gz \
        --pair-map data/epitopes/iedb/netmhcpan/mhc_i/pair_id_mapping.csv.gz \
        --out-dir results/netmhcpan/iedb/mhc_i/parsed/ \
        --rank-threshold 2.0

    # Class II
    python parse_netmhcpan_results.py \
        --mhc-class II \
        --result-dir results/netmhcpan/iedb/mhc_ii/ \
        --pairs-file data/epitopes/iedb/blast/mhc_ii/iedb_mhc_ii_filtered_blast_hits.csv.gz \
        --pair-map data/epitopes/iedb/netmhcpan/mhc_ii/pair_id_mapping.csv.gz \
        --out-dir results/netmhcpan/iedb/mhc_ii/parsed/ \
        --rank-threshold 5.0
======================================================================
"""

import argparse
import pandas as pd
import os
import glob


# ====================================================
# Column name standardization
# ====================================================
# Map various column names to a standard set
COLUMN_ALIASES = {
    # BA score
    "BA_score": "ba_score",
    "Score_BA": "ba_score",
    # BA rank
    "BA_Rank": "ba_rank",
    "Rank_BA": "ba_rank",
    # EL score
    "Score": "el_score",
    # EL rank
    "Rank": "el_rank",
    # Binding core
    "core": "core",
    "Core": "core",
    # Inner core
    "icore": "icore",
    # Affinity
    "Aff(nM)": "affinity_nM",
    "nM": "affinity_nM",
    # Average
    "Ave": "ave",
    # Number of binders
    "NB": "nb",
    # Identity
    "Identity": "ID",
}


def standardize_columns(df):
    """Rename columns to standard names, keeping originals if no alias."""
    rename_map = {}
    for col in df.columns:
        if col in COLUMN_ALIASES:
            rename_map[col] = COLUMN_ALIASES[col]
    return df.rename(columns=rename_map)


# ====================================================
# Step 1: Parse a single xls file
# ====================================================
def parse_xls(filepath):
    """
    Parse one NetMHCpan or NetMHCIIpan xls output file.

    Format:
      Line 1: allele name (tab-prefixed)
      Line 2: column headers
      Lines 3+: data rows
    """
    rows = []
    header = None
    allele_name = None

    with open(filepath) as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            if not line:
                continue

            # Skip comment lines
            if line.startswith("#"):
                continue

            parts = line.split("\t")

            # Line 1: allele name
            if allele_name is None and header is None:
                # First non-empty, non-comment line — allele header
                allele_name = parts[-1].strip() if parts else None
                continue

            # Line 2: column headers
            if header is None:
                header = [h.strip() for h in parts]
                continue

            # Data rows
            if len(parts) >= len(header) - 2:
                rows.append(parts)

    if not rows or header is None:
        return pd.DataFrame(), allele_name

    df = pd.DataFrame(rows)
    n_cols = min(len(header), df.shape[1])
    df = df.iloc[:, :n_cols]
    df.columns = header[:n_cols]

    return df, allele_name


# ====================================================
# Step 2: Parse all result files
# ====================================================
def parse_all_results(result_dir):
    print("=" * 60)
    print("Step 1: Parsing results")
    print("=" * 60)

    xls_files = sorted(glob.glob(os.path.join(result_dir, "*.xls")))
    print(f"  Found {len(xls_files)} xls files")

    if not xls_files:
        print("  ⚠️ No xls files found!")
        return pd.DataFrame()

    all_results = []
    failed = []

    for i, xls_file in enumerate(xls_files):
        allele_file_name = os.path.basename(xls_file).replace(".xls", "")

        try:
            df, allele_from_header = parse_xls(xls_file)
            if df.empty:
                failed.append(xls_file)
                continue

            df["allele_file"] = allele_file_name
            df["mhc_allele"] = allele_from_header or allele_file_name
            all_results.append(df)
        except Exception as e:
            failed.append(xls_file)
            if i < 5:
                print(f"    ⚠️ Failed to parse {os.path.basename(xls_file)}: {e}")

        if (i + 1) % 50 == 0:
            print(f"    Parsed {i + 1}/{len(xls_files)} files...")

    if failed:
        print(f"  ⚠️ {len(failed)} files failed to parse")

    if not all_results:
        print("  ❌ No results parsed!")
        return pd.DataFrame()

    results = pd.concat(all_results, ignore_index=True)

    # Standardize column names
    results = standardize_columns(results)

    # Convert numeric columns
    for col in ["el_score", "el_rank", "ba_score", "ba_rank", "affinity_nM", "ave"]:
        if col in results.columns:
            results[col] = pd.to_numeric(results[col], errors="coerce")
    if "Pos" in results.columns:
        results["Pos"] = pd.to_numeric(results["Pos"], errors="coerce")
    if "nb" in results.columns:
        results["nb"] = pd.to_numeric(results["nb"], errors="coerce")

    print(f"  ✅ Parsed {len(results):,} total predictions from "
          f"{len(all_results)} files")
    print(f"  Standardized columns: {list(results.columns)}")

    return results


# ====================================================
# Step 3: Extract pair identity and keep best prediction
# ====================================================
def extract_pair_info(results):
    """
    Parse ID column (V_XXXX / H_XXXX) and keep best-binding
    prediction per (ID, allele).

    For Class I: best across length windows (-l 8,...,14)
    For Class II: best across binding core positions
    """
    print(f"\n{'=' * 60}")
    print("Step 2: Extracting pair identity")
    print("=" * 60)

    if "ID" not in results.columns:
        print("  ⚠️ No 'ID' column — available columns:")
        print(f"    {list(results.columns)}")
        return results

    # Parse peptide type
    results["peptide_type"] = results["ID"].apply(
        lambda x: "VIRAL" if str(x).startswith("V_") else
                  ("HUMAN" if str(x).startswith("H_") else "UNKNOWN")
    )

    # Extract pair_id
    results["pair_id"] = results["ID"].apply(
        lambda x: str(x).split("_")[1] if "_" in str(x) else None
    )

    # Keep best-binding prediction per (ID, allele_file)
    # Sort by BA rank (strongest binding affinity)
    results["peptide_len"] = results["Peptide"].str.len()

    n_before = len(results)
    results = results.sort_values("ba_rank").drop_duplicates(
        subset=["ID", "allele_file"], keep="first"
    )
    n_removed = n_before - len(results)
    if n_removed > 0:
        print(f"  Kept best-binding prediction per (ID, allele) "
              f"by BA rank (removed {n_removed:,} sub-optimal predictions)")

    # Summary
    type_counts = results["peptide_type"].value_counts()
    print(f"  Peptide types:")
    for ptype, count in type_counts.items():
        print(f"    {ptype}: {count:,}")
    print(f"  Unique pair_ids: {results['pair_id'].nunique():,}")
    print(f"  Unique alleles: {results['mhc_allele'].nunique()}")

    return results


# ====================================================
# Step 4: Match viral-human pairs and compute ΔBA
# ====================================================
def compute_delta_ba(results):
    print(f"\n{'=' * 60}")
    print("Step 3: Computing ΔBA for viral-human pairs")
    print("=" * 60)

    viral = results[results["peptide_type"] == "VIRAL"].copy()
    human = results[results["peptide_type"] == "HUMAN"].copy()

    print(f"  Viral predictions: {len(viral):,}")
    print(f"  Human predictions: {len(human):,}")

    if viral.empty or human.empty:
        print("  ⚠️ Missing viral or human predictions!")
        return pd.DataFrame()

    # Rename for merge
    viral_rename = {
        "Peptide": "viral_peptide",
        "ba_rank": "viral_rank_BA",
        "el_rank": "viral_rank_EL",
        "ba_score": "viral_score_BA",
        "el_score": "viral_score_EL",
        "affinity_nM": "viral_nM",
        "nb": "viral_NB",
        "core": "viral_core",
        "peptide_len": "viral_peptide_len",
    }
    viral_r = viral.rename(columns=viral_rename)

    human_rename = {
        "Peptide": "human_peptide",
        "ba_rank": "human_rank_BA",
        "el_rank": "human_rank_EL",
        "ba_score": "human_score_BA",
        "el_score": "human_score_EL",
        "affinity_nM": "human_nM",
        "nb": "human_NB",
        "core": "human_core",
        "peptide_len": "human_peptide_len",
    }
    human_r = human.rename(columns=human_rename)

    viral_keep = ["pair_id", "allele_file", "mhc_allele",
                  "viral_peptide", "viral_rank_BA", "viral_rank_EL",
                  "viral_score_BA", "viral_score_EL", "viral_nM",
                  "viral_NB", "viral_core", "viral_peptide_len"]
    viral_keep = [c for c in viral_keep if c in viral_r.columns]

    human_keep = ["pair_id", "allele_file",
                  "human_peptide", "human_rank_BA", "human_rank_EL",
                  "human_score_BA", "human_score_EL", "human_nM",
                  "human_NB", "human_core", "human_peptide_len"]
    human_keep = [c for c in human_keep if c in human_r.columns]

    merged = viral_r[viral_keep].merge(
        human_r[human_keep],
        on=["pair_id", "allele_file"],
        how="inner"
    )

    print(f"  Matched pairs: {len(merged):,}")

    if merged.empty:
        print("  ⚠️ No pairs matched!")
        return merged

    # Ensure numeric
    for col in ["viral_rank_BA", "viral_rank_EL", "viral_score_BA", "viral_score_EL",
                "human_rank_BA", "human_rank_EL", "human_score_BA", "human_score_EL",
                "viral_nM", "human_nM"]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    # ΔBA
    merged["delta_rank_BA"] = merged["viral_rank_BA"] - merged["human_rank_BA"]
    merged["delta_rank_EL"] = merged["viral_rank_EL"] - merged["human_rank_EL"]
    if "viral_score_BA" in merged.columns and "human_score_BA" in merged.columns:
        merged["delta_score_BA"] = merged["viral_score_BA"] - merged["human_score_BA"]

    return merged


# ====================================================
# Step 5: Categorize pairs
# ====================================================
def categorize_pairs(merged, rank_threshold, delta_threshold):
    print(f"\n{'=' * 60}")
    print("Step 4: Categorizing pairs")
    print(f"  Binder threshold: BA_Rank ≤ {rank_threshold}%")
    print(f"  ΔBA threshold: |ΔBA_score| < {delta_threshold}")
    print("=" * 60)

    def classify(row):
        viral_rank = row.get("viral_rank_BA")
        viral_score = row.get("viral_score_BA")
        human_score = row.get("human_score_BA")

        if pd.isna(viral_rank):
            return "Unknown"

        if viral_rank > rank_threshold:
            return "Non-binder"

        if pd.isna(viral_score) or pd.isna(human_score):
            return "Unknown"

        delta_ba = viral_score - human_score

        if delta_ba >= delta_threshold:
            return "Viral-dominant"
        elif delta_ba <= -delta_threshold:
            return "Human-dominant"
        else:
            return "Equivalent-binding"

    merged["category"] = merged.apply(classify, axis=1)
    merged["delta_BA_score"] = merged["viral_score_BA"] - merged["human_score_BA"]
    merged["human_is_binder"] = merged["human_rank_BA"] <= rank_threshold

    # Summary
    print(f"\n  Category distribution:")
    cats = merged["category"].value_counts()
    for cat, count in cats.items():
        pct = count / len(merged) * 100
        print(f"    {cat}: {count:,} ({pct:.1f}%)")

    print(f"\n  Human peptide binding (BA_Rank ≤ {rank_threshold}):")
    hb = merged["human_is_binder"].value_counts()
    for is_binder, count in hb.items():
        label = "Binder" if is_binder else "Non-binder"
        print(f"    {label}: {count:,}")

    strong = merged[
        (merged["category"] == "Equivalent-binding") &
        (merged["human_is_binder"])
    ]
    print(f"\n  Strong mimicry candidates (equivalent + both bind): {len(strong):,}")

    return merged


# ====================================================
# Step 6: Merge with original pair metadata
# ====================================================
def merge_metadata(merged, pairs_file, pair_map_file):
    print(f"\n{'=' * 60}")
    print("Step 5: Merging with original pair metadata")
    print("=" * 60)

    # Load pair ID map
    pair_map = pd.read_csv(pair_map_file, low_memory=False)
    print(f"  Pair ID map: {len(pair_map):,} entries")

    final = merged.merge(pair_map, on="pair_id", how="left")
    print(f"  After pair_id merge: {len(final):,}")

    # Load original pairs for additional metadata
    pairs = pd.read_csv(pairs_file, low_memory=False)
    print(f"  Original pairs: {len(pairs):,}")

    metadata_cols = [
        "structure_id", "hu_prot_id", "hu_prot_name",
        "source_organism", "mhc_allele", "n_mismatches",
        "pident", "coverage", "sstart", "send", "qlen"
    ]
    metadata_cols = [c for c in metadata_cols if c in pairs.columns]

    # Ensure types match for join
    for col in ["structure_id"]:
        if col in final.columns:
            final[col] = pd.to_numeric(final[col], errors="coerce")
        if col in pairs.columns:
            pairs[col] = pd.to_numeric(pairs[col], errors="coerce")

    join_cols = ["structure_id", "hu_prot_id"]
    join_cols = [c for c in join_cols if c in final.columns and c in pairs.columns]

    if join_cols:
        extra_cols = [c for c in metadata_cols if c not in join_cols and c not in final.columns]
        if extra_cols:
            final = final.merge(
                pairs[join_cols + extra_cols].drop_duplicates(join_cols),
                on=join_cols,
                how="left",
                suffixes=("", "_orig")
            )

    final = final.drop_duplicates(subset=["pair_id", "allele_file"])
    print(f"  Final merged records: {len(final):,}")

    return final


# ====================================================
# Step 7: Save results
# ====================================================
def save_results(final, out_dir, mhc_class, rank_threshold):
    print(f"\n{'=' * 60}")
    print("Step 6: Saving results")
    print("=" * 60)

    prefix = f"iedb_mhc{'i' if mhc_class == 'I' else 'ii'}"

    full_path = os.path.join(out_dir, f"{prefix}_netmhcpan_results.csv.gz")
    final.to_csv(full_path, index=False, compression="gzip")
    print(f"  ✅ Full results: {len(final):,} → {full_path}")

    binders = final[final["viral_rank_BA"] <= rank_threshold]
    binder_path = os.path.join(out_dir, f"{prefix}_binders.csv.gz")
    binders.to_csv(binder_path, index=False, compression="gzip")
    print(f"  ✅ Binders (BA_Rank ≤ {rank_threshold}): "
          f"{len(binders):,} → {binder_path}")

    mimicry = final[
        (final["viral_rank_BA"] <= rank_threshold) &
        (final["human_is_binder"])
    ]
    mimicry_path = os.path.join(out_dir, f"{prefix}_mimicry_candidates.csv.gz")
    mimicry.to_csv(mimicry_path, index=False, compression="gzip")
    print(f"  ✅ Mimicry candidates (both bind): {len(mimicry):,} → {mimicry_path}")

    strong = mimicry[mimicry["category"] == "Equivalent-binding"]
    strong_path = os.path.join(out_dir, f"{prefix}_strong_mimicry.csv.gz")
    strong.to_csv(strong_path, index=False, compression="gzip")
    print(f"  ✅ Strong mimicry (equivalent binding): {len(strong):,} → {strong_path}")

    return final


# ====================================================
# Summary
# ====================================================
def print_summary(final, mhc_class, rank_threshold):
    print(f"\n{'=' * 60}")
    print(f"Final Summary — MHC Class {mhc_class}")
    print("=" * 60)

    print(f"  Total pairs scored: {len(final):,}")
    print(f"  Unique viral epitopes: {final['viral_peptide'].nunique():,}")
    print(f"  Unique human mimics: {final['human_peptide'].nunique():,}")
    print(f"  Unique pair IDs: {final['pair_id'].nunique():,}")
    print(f"  Unique HLA alleles: {final['allele_file'].nunique()}")
    print(f"  Binder threshold: BA_Rank ≤ {rank_threshold}% + ΔBA_score for categorization")

    print(f"\n  Category breakdown:")
    for cat in ["Viral-dominant", "Human-dominant", "Equivalent-binding", "Non-binder", "Unknown"]:
        subset = final[final["category"] == cat]
        if not subset.empty:
            print(f"    {cat}: {len(subset):,} pairs, "
                  f"{subset['viral_peptide'].nunique():,} epitopes")

    if "source_organism" in final.columns:
        mimicry = final[final["human_is_binder"] & (final["viral_rank_BA"] <= rank_threshold)]
        if not mimicry.empty:
            print(f"\n  Top organisms in mimicry candidates:")
            orgs = mimicry["source_organism"].value_counts().head(10)
            for org, count in orgs.items():
                print(f"    {org}: {count}")


# ====================================================
# Main
# ====================================================
def main():
    parser = argparse.ArgumentParser(
        description="Parse NetMHCpan/NetMHCIIpan results for IEDB mimicry pairs"
    )
    parser.add_argument("--mhc-class", required=True, choices=["I", "II"],
                        help="MHC class (I or II)")
    parser.add_argument("--result-dir", required=True,
                        help="Directory containing .xls output files")
    parser.add_argument("--pairs-file", required=True,
                        help="Original filtered BLAST pairs CSV")
    parser.add_argument("--pair-map", required=True,
                        help="Pair ID mapping file (pair_id_mapping.csv.gz)")
    parser.add_argument("--out-dir", required=True,
                        help="Output directory for parsed results")
    parser.add_argument("--rank-threshold", type=float, default=None,
                        help="Binder rank threshold (default: 2.0 for Class I, 5.0 for Class II)")
    parser.add_argument("--delta-threshold", type=float, default=0.5,
                        help="ΔBA_score threshold for categorization (default: 0.5)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Set default rank threshold based on class
    # EL rank thresholds: Class I ≤ 0.5% strong / ≤ 2% weak
    #                     Class II ≤ 2% strong / ≤ 5% weak
    # Using weak binder cutoff as default to be inclusive
    rank_threshold = args.rank_threshold
    if rank_threshold is None:
        rank_threshold = 2.0 if args.mhc_class == "I" else 5.0

    print("=" * 60)
    print(f"Parse NetMHC{'II' if args.mhc_class == 'II' else ''}pan Results")
    print(f"MHC Class {args.mhc_class}")
    print(f"Binder threshold: BA_Rank ≤ {rank_threshold}%")
    print("=" * 60)

    # Step 1: Parse
    results = parse_all_results(args.result_dir)
    if results.empty:
        print("\n❌ No results to process.")
        exit(1)

    # Step 2: Extract pair identity
    results = extract_pair_info(results)

    # Step 3: Match pairs and compute ΔBA
    merged = compute_delta_ba(results)
    if merged.empty:
        print("\n❌ No pairs matched.")
        exit(1)

    # Step 4: Categorize
    merged = categorize_pairs(merged, rank_threshold, args.delta_threshold)

    # Step 5: Merge metadata
    final = merge_metadata(merged, args.pairs_file, args.pair_map)

    # Step 6: Save
    final = save_results(final, args.out_dir, args.mhc_class, rank_threshold)

    # Summary
    print_summary(final, args.mhc_class, rank_threshold)

    print(f"\n✅ Done.")


if __name__ == "__main__":
    main()