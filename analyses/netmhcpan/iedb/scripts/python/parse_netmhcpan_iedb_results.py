#!/usr/bin/env python3
"""
======================================================================
Script: parse_netmhcpan_iedb_results.py
Description:
    Parse NetMHCpan xls output files for IEDB mimicry pairs.
    Match viral-human peptide pairs, compute ΔBA, and categorize.

    NetMHCpan xls output columns (tab-separated, header rows start with #):
      Pos, MHC, Peptide, Core, Of, Gp, Gl, Ip, Il, Icore, Identity,
      Score_EL, %Rank_EL, Score_BA, %Rank_BA, Aff(nM), BindLevel

    FASTA headers encode pair identity:
      >VIRAL|structure_id|sequence
      >HUMAN|structure_id|hu_prot_id|sequence

    ΔBA categories:
      Viral-dominant:    viral binds (Rank ≤ 2) & ΔBA ≤ -0.5
      Human-dominant:    viral binds (Rank ≤ 2) & ΔBA ≥ +0.5
      Equivalent:        viral binds (Rank ≤ 2) & |ΔBA| < 0.5
      Non-binder:        viral Rank > 2

Input:
    - NetMHCpan xls files in results directory
    - Original pairs file for metadata
    - Allele manifest for mapping

Output:
    - Merged results with ΔBA scores
    - Categorized mimicry candidates

Dependencies:
    pandas, glob

Usage:
    python parse_netmhcpan_iedb_results.py
======================================================================
"""

import pandas as pd
import os
import glob
import re

# ====================================================
# Config
# ====================================================
BASE_DIR = "/ix/djishnu/Priyamvada/virauto"
RESULT_DIR = os.path.join(BASE_DIR, "results/netmhcpan/iedb/mhc_i")
PAIRS_FILE = os.path.join(BASE_DIR, "data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla.csv.gz")
MANIFEST_FILE = os.path.join(BASE_DIR, "data/epitopes/iedb/netmhcpan/allele_manifest.tsv")
PAIR_ID_MAP_FILE = os.path.join(BASE_DIR, "data/epitopes/iedb/netmhcpan/mhc_i/pair_id_mapping.csv.gz")

OUT_DIR = os.path.join(BASE_DIR, "results/netmhcpan/iedb/mhc_i/parsed")
os.makedirs(OUT_DIR, exist_ok=True)

# ΔBA thresholds
RANK_BINDER_THRESHOLD = 2.0    # %Rank_BA ≤ 2 = binder
DELTA_BA_THRESHOLD = 0.5       # |ΔBA| threshold for categorization


# ====================================================
# Step 1: Parse a single NetMHCpan xls file
# ====================================================
def parse_xls(filepath):
    """
    Parse one NetMHCpan xls output file.
    
    The xls format is tab-separated with:
      - Header lines starting with '#'
      - Column names in the first non-comment line
      - Data rows following
    
    Returns DataFrame with columns:
      Pos, MHC, Peptide, Core, Identity, Score_EL, %Rank_EL, 
      Score_BA, %Rank_BA, Aff(nM), BindLevel
    """
    rows = []
    header = None

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Skip comment lines but capture column header
            if line.startswith("#"):
                continue

            parts = line.split("\t")

            # First non-comment line with many fields is the header
            if header is None and len(parts) > 5:
                # Check if this looks like a header (contains known column names)
                if "Peptide" in parts or "MHC" in parts:
                    header = parts
                    continue
                # Sometimes the header doesn't have labels — 
                # use known NetMHCpan column order
                elif any(p.replace(".", "").replace("-", "").isdigit() for p in parts[5:8]):
                    # This is a data row, infer header from position
                    header = [
                        "Pos", "MHC", "Peptide", "Core", "Of", "Gp", "Gl",
                        "Ip", "Il", "Icore", "Identity", "Score_EL",
                        "%Rank_EL", "Score_BA", "%Rank_BA", "Aff(nM)", "BindLevel"
                    ]
                    # Don't skip this line — it's data
                    rows.append(parts)
                    continue

            if header is not None and len(parts) >= len(header) - 2:
                rows.append(parts)

    if not rows or header is None:
        return pd.DataFrame()

    # Build DataFrame, handling rows with fewer/more columns than header
    df = pd.DataFrame(rows)

    # Assign header names to available columns
    n_cols = min(len(header), df.shape[1])
    df = df.iloc[:, :n_cols]
    df.columns = header[:n_cols]

    # Convert numeric columns
    for col in ["%Rank_EL", "%Rank_BA", "Score_EL", "Score_BA", "Aff(nM)"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


# ====================================================
# Step 2: Parse all result files
# ====================================================
def parse_all_results(result_dir):
    print("=" * 60)
    print("Step 1: Parsing NetMHCpan results")
    print("=" * 60)

    xls_files = sorted(glob.glob(os.path.join(result_dir, "*.xls")))
    print(f"  Found {len(xls_files)} xls files")

    if not xls_files:
        print("  ⚠️ No xls files found!")
        return pd.DataFrame()

    all_results = []
    failed = []

    for i, xls_file in enumerate(xls_files):
        allele_name = os.path.basename(xls_file).replace(".xls", "")

        try:
            df = parse_xls(xls_file)
            if df.empty:
                failed.append(xls_file)
                continue

            df["allele_file"] = allele_name
            all_results.append(df)
        except Exception as e:
            failed.append(xls_file)
            if i < 5:  # Only print first few errors
                print(f"    ⚠️ Failed to parse {os.path.basename(xls_file)}: {e}")

        if (i + 1) % 50 == 0:
            print(f"    Parsed {i + 1}/{len(xls_files)} files...")

    if failed:
        print(f"  ⚠️ {len(failed)} files failed to parse")

    if not all_results:
        print("  ❌ No results parsed!")
        return pd.DataFrame()

    results = pd.concat(all_results, ignore_index=True)
    print(f"  ✅ Parsed {len(results):,} total predictions from "
          f"{len(all_results)} files")

    return results


# ====================================================
# Step 3: Extract pair identity from FASTA headers
# ====================================================
def extract_pair_info(results):
    """
    Parse the Identity column to identify VIRAL vs HUMAN peptides
    and extract structure_id for matching pairs.
    
    FASTA header formats:
      VIRAL|structure_id|sequence
      HUMAN|structure_id|hu_prot_id|sequence
    """
    print(f"\n{'=' * 60}")
    print("Step 2: Extracting pair identity")
    print("=" * 60)

    if "ID" not in results.columns:
        print("  ⚠️ No 'ID' column — checking column names:")
        print(f"    {list(results.columns)}")
        return results

    # Parse ID field
    # Format: V_XXXX or H_XXXX (4-char alphanumeric pair code)
    results["peptide_type"] = results["ID"].apply(
        lambda x: "VIRAL" if str(x).startswith("V_") else
                  ("HUMAN" if str(x).startswith("H_") else "UNKNOWN")
    )

    # Extract pair_id: the 4-char code after V_ or H_
    # V_A01B → A01B
    # H_A01B → A01B
    def extract_pair_id(id_str):
        parts = str(id_str).split("_")
        if len(parts) >= 2:
            return parts[1]
        return None

    results["pair_id"] = results["ID"].apply(extract_pair_id)

    # --- Keep only best-binding prediction per peptide ID ---
    # NetMHCpan with -l 8,9,10,11,12,13,14 scores every possible
    # sub-window and length for each input sequence. We want the 
    # strongest binding prediction (lowest BA_Rank) for each ID,
    # which represents the optimal binding register for that peptide.
    results["peptide_len"] = results["Peptide"].str.len()
    
    # Ensure BA_Rank is numeric for sorting
    results["BA_Rank"] = pd.to_numeric(results.get("BA_Rank", pd.Series()), errors="coerce")
    
    n_before = len(results)
    # For each ID + allele_file, keep the row with the lowest BA_Rank
    results = results.sort_values("BA_Rank").drop_duplicates(
        subset=["ID", "allele_file"], keep="first"
    )
    n_removed = n_before - len(results)
    if n_removed > 0:
        print(f"  Kept best-binding prediction per ID "
              f"(removed {n_removed:,} sub-optimal windows)")

    # Summary
    type_counts = results["peptide_type"].value_counts()
    print(f"  Peptide types:")
    for ptype, count in type_counts.items():
        print(f"    {ptype}: {count:,}")

    print(f"  Unique pair_ids: {results['pair_id'].nunique():,}")

    return results


# ====================================================
# Step 4: Match viral-human pairs and compute ΔBA
# ====================================================
def compute_delta_ba(results):
    """
    For each structure_id + allele combination, match the viral and 
    human predictions and compute ΔBA = Rank_BA(viral) - Rank_BA(human).
    
    ΔBA < 0 → viral binds stronger
    ΔBA > 0 → human binds stronger
    ΔBA ≈ 0 → equivalent binding
    """
    print(f"\n{'=' * 60}")
    print("Step 3: Computing ΔBA for viral-human pairs")
    print("=" * 60)

    # Split into viral and human
    viral = results[results["peptide_type"] == "VIRAL"].copy()
    human = results[results["peptide_type"] == "HUMAN"].copy()

    print(f"  Viral predictions: {len(viral):,}")
    print(f"  Human predictions: {len(human):,}")

    if viral.empty or human.empty:
        print("  ⚠️ Missing viral or human predictions!")
        return pd.DataFrame()

    # Rename columns for merge
    viral_cols = {
        "Peptide": "viral_peptide",
        "BA_Rank": "viral_rank_BA",
        "Rank": "viral_rank_EL",
        "BA_score": "viral_score_BA",
        "Score": "viral_score_EL",
        "NB": "viral_NB",
        "MHC": "mhc_allele",
        "peptide_len": "viral_peptide_len",
    }
    viral_renamed = viral.rename(columns=viral_cols)

    human_cols = {
        "Peptide": "human_peptide",
        "BA_Rank": "human_rank_BA",
        "Rank": "human_rank_EL",
        "BA_score": "human_score_BA",
        "Score": "human_score_EL",
        "NB": "human_NB",
        "peptide_len": "human_peptide_len",
    }
    human_renamed = human.rename(columns=human_cols)

    # Select columns for merge
    viral_keep = ["pair_id", "allele_file", "mhc_allele",
                  "viral_peptide", "viral_rank_BA", "viral_rank_EL",
                  "viral_score_BA", "viral_score_EL", "viral_NB",
                  "viral_peptide_len"]
    viral_keep = [c for c in viral_keep if c in viral_renamed.columns]

    human_keep = ["pair_id", "allele_file",
                  "human_peptide", "human_rank_BA", "human_rank_EL",
                  "human_score_BA", "human_score_EL", "human_NB",
                  "human_peptide_len"]
    human_keep = [c for c in human_keep if c in human_renamed.columns]

    # Merge on pair_id + allele_file
    # pair_id uniquely identifies a viral-human pair, so this is 
    # a clean one-to-one merge per allele
    merged = viral_renamed[viral_keep].merge(
        human_renamed[human_keep],
        on=["pair_id", "allele_file"],
        how="inner"
    )

    print(f"  Matched pairs: {len(merged):,}")

    if merged.empty:
        print("  ⚠️ No pairs matched!")
        return merged

    # Ensure numeric types after merge
    numeric_cols = [
        "viral_rank_BA", "viral_rank_EL", "viral_score_BA", "viral_score_EL",
        "human_rank_BA", "human_rank_EL", "human_score_BA", "human_score_EL",
    ]
    for col in numeric_cols:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    # Compute ΔBA (using BA_Rank)
    merged["delta_rank_BA"] = merged["viral_rank_BA"] - merged["human_rank_BA"]
    merged["delta_rank_EL"] = merged["viral_rank_EL"] - merged["human_rank_EL"]

    # Score difference (higher score = stronger binding)
    if "viral_score_BA" in merged.columns and "human_score_BA" in merged.columns:
        merged["delta_score_BA"] = merged["viral_score_BA"] - merged["human_score_BA"]

    return merged


# ====================================================
# Step 5: Categorize pairs
# ====================================================
def categorize_pairs(merged):
    """
    Categorize based on binding and ΔBA:
      Viral-dominant:  viral binds (Rank ≤ 2) & ΔBA ≤ -threshold
      Human-dominant:  viral binds (Rank ≤ 2) & ΔBA ≥ +threshold
      Equivalent:      viral binds (Rank ≤ 2) & |ΔBA| < threshold
      Non-binder:      viral Rank > 2
    """
    print(f"\n{'=' * 60}")
    print("Step 4: Categorizing pairs")
    print("=" * 60)

    def classify(row):
        viral_rank = row.get("viral_rank_BA")
        viral_score = row.get("viral_score_BA")
        human_score = row.get("human_score_BA")

        if pd.isna(viral_rank):
            return "Unknown"

        # Viral must bind (BA_Rank ≤ 2)
        if viral_rank > RANK_BINDER_THRESHOLD:
            return "Non-binder"

        if pd.isna(viral_score) or pd.isna(human_score):
            return "Unknown"

        # ΔBA = BA_score(viral) - BA_score(human)
        # Higher BA_score = stronger binding
        # ΔBA < 0 → human binds stronger
        # ΔBA > 0 → viral binds stronger
        delta_ba = viral_score - human_score

        if delta_ba >= DELTA_BA_THRESHOLD:
            return "Viral-dominant"
        elif delta_ba <= -DELTA_BA_THRESHOLD:
            return "Human-dominant"
        else:
            return "Equivalent-binding"

    merged["category"] = merged.apply(classify, axis=1)

    # Also compute ΔBA_score column for reference
    merged["delta_BA_score"] = merged["viral_score_BA"] - merged["human_score_BA"]

    # Also classify human binding independently
    merged["human_is_binder"] = merged["human_rank_BA"] <= RANK_BINDER_THRESHOLD

    # Summary
    print(f"  Category distribution:")
    cats = merged["category"].value_counts()
    for cat, count in cats.items():
        pct = count / len(merged) * 100
        print(f"    {cat}: {count:,} ({pct:.1f}%)")

    print(f"\n  Human peptide binding (BA_Rank ≤ {RANK_BINDER_THRESHOLD}):")
    hb = merged["human_is_binder"].value_counts()
    for is_binder, count in hb.items():
        label = "Binder" if is_binder else "Non-binder"
        print(f"    {label}: {count:,}")

    # Strongest mimicry candidates: equivalent binding where both bind
    strong = merged[
        (merged["category"] == "Equivalent-binding") &
        (merged["human_is_binder"])
    ]
    print(f"\n  Strong mimicry candidates (equivalent + both bind): {len(strong):,}")

    return merged


# ====================================================
# Step 6: Merge with original pair metadata
# ====================================================
def merge_metadata(merged, pairs_file):
    print(f"\n{'=' * 60}")
    print("Step 5: Merging with original pair metadata")
    print("=" * 60)

    # Load pair ID map (generated by prepare_netmhcpan_iedb.py)
    pair_map = pd.read_csv(PAIR_ID_MAP_FILE, low_memory=False)
    print(f"  Pair ID map: {len(pair_map):,} entries")

    # Merge pair_id → structure_id, hu_prot_id, sequences
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

    # Ensure types match
    for col in ["structure_id"]:
        if col in final.columns:
            final[col] = pd.to_numeric(final[col], errors="coerce")
        if col in pairs.columns:
            pairs[col] = pd.to_numeric(pairs[col], errors="coerce")

    # Merge on structure_id + hu_prot_id for precise matching
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

    # Deduplicate
    final = final.drop_duplicates(
        subset=["pair_id", "allele_file"]
    )

    print(f"  Final merged records: {len(final):,}")

    return final


# ====================================================
# Step 7: Save results
# ====================================================
def save_results(final):
    print(f"\n{'=' * 60}")
    print("Step 6: Saving results")
    print("=" * 60)

    # Full results
    full_path = os.path.join(OUT_DIR, "iedb_mhci_netmhcpan_results.csv.gz")
    final.to_csv(full_path, index=False, compression="gzip")
    print(f"  ✅ Full results: {len(final):,} → {full_path}")

    # Binders only (viral Rank ≤ 2)
    binders = final[final["viral_rank_BA"] <= RANK_BINDER_THRESHOLD]
    binder_path = os.path.join(OUT_DIR, "iedb_mhci_binders.csv.gz")
    binders.to_csv(binder_path, index=False, compression="gzip")
    print(f"  ✅ Binders (viral BA_Rank ≤ {RANK_BINDER_THRESHOLD}): "
          f"{len(binders):,} → {binder_path}")

    # Mimicry candidates (both bind)
    mimicry = final[
        (final["viral_rank_BA"] <= RANK_BINDER_THRESHOLD) &
        (final["human_is_binder"])
    ]
    mimicry_path = os.path.join(OUT_DIR, "iedb_mhci_mimicry_candidates.csv.gz")
    mimicry.to_csv(mimicry_path, index=False, compression="gzip")
    print(f"  ✅ Mimicry candidates (both bind): {len(mimicry):,} → {mimicry_path}")

    # Strong candidates (equivalent binding)
    strong = mimicry[mimicry["category"] == "Equivalent-binding"]
    strong_path = os.path.join(OUT_DIR, "iedb_mhci_strong_mimicry.csv.gz")
    strong.to_csv(strong_path, index=False, compression="gzip")
    print(f"  ✅ Strong mimicry (equivalent binding): {len(strong):,} → {strong_path}")

    return final


# ====================================================
# Summary
# ====================================================
def print_summary(final):
    print(f"\n{'=' * 60}")
    print("Final Summary")
    print("=" * 60)

    print(f"  Total pairs scored: {len(final):,}")
    print(f"  Unique viral epitopes: {final['viral_peptide'].nunique():,}")
    print(f"  Unique human mimics: {final['human_peptide'].nunique():,}")
    print(f"  Unique pair IDs: {final['pair_id'].nunique():,}")
    print(f"  Unique HLA alleles: {final['allele_file'].nunique()}")

    print(f"\n  Category breakdown:")
    for cat in ["Viral-dominant", "Human-dominant", "Equivalent-binding", "Non-binder", "Unknown"]:
        subset = final[final["category"] == cat]
        if not subset.empty:
            print(f"    {cat}: {len(subset):,} pairs, "
                  f"{subset['viral_peptide'].nunique():,} epitopes")

    # Top organisms among mimicry candidates
    if "source_organism" in final.columns:
        mimicry = final[final["human_is_binder"] & (final["viral_rank_BA"] <= RANK_BINDER_THRESHOLD)]
        if not mimicry.empty:
            print(f"\n  Top organisms in mimicry candidates:")
            orgs = mimicry["source_organism"].value_counts().head(10)
            for org, count in orgs.items():
                print(f"    {org}: {count}")

    print(f"\n  Next steps:")
    print(f"    1. Cross-reference mimicry candidates with TCR data")
    print(f"    2. Run DecoderTCR embedding analysis on candidates with paired TCRs")
    print(f"    3. Check human protein expression in autoimmune-relevant tissues")


# ====================================================
# Main
# ====================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Parse NetMHCpan Results — IEDB Mimicry Pairs")
    print("=" * 60)

    # Step 1: Parse all xls files
    results = parse_all_results(RESULT_DIR)

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
    merged = categorize_pairs(merged)

    # Step 5: Merge metadata
    final = merge_metadata(merged, PAIRS_FILE)

    # Step 6: Save
    final = save_results(final)

    # Summary
    print_summary(final)

    print(f"\n✅ Done.")