#!/usr/bin/env python3
"""
======================================================================
Script: filter_iedb_mhcii_blast_hits.py
Description:
    Filter BLAST results from IEDB MHC-II epitopes vs human proteome.
    Extract human mimic subsequences and prepare for NetMHCIIpan.

    Human mimic sequences are padded with flanking residues from the
    source protein to match the viral epitope length (minimum 15 aa),
    ensuring reliable NetMHCIIpan predictions.

Filtering steps:
    1. Epitope length: >= 15 aa (Class II standard)
    2. Coverage: alignment must cover >= 80% of epitope length
    3. Mismatches: 1 to qlen//3 mismatches
    4. Remove deprecated UniProt IDs
    5. Remove hits to immunoglobulin/TCR/MHC proteins
    6. Extract human mimic subsequences (padded to viral epitope length)
    7. Parse MHC allele from query header for downstream NetMHCIIpan

Input FASTA header format (from merge script):
    >structure_id|sequence|organism|mhc_allele

Dependencies:
    pandas, biopython

Usage:
    python filter_iedb_mhcii_blast_hits.py
======================================================================
"""

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import re

# ====================================================
# Config
# ====================================================
BLAST_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/mhc_ii/iedb_mhc_ii_vs_human_proteome.tsv"
HUMAN_FASTA = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_sprot.fasta"
UNIPROT_FILTER = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/proteins_to_remove_from_UniProtKB.txt"

OUT_DIR = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/mhc_ii"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_CSV = os.path.join(OUT_DIR, "iedb_mhc_ii_filtered_blast_hits.csv.gz")
OUT_FASTA = os.path.join(OUT_DIR, "iedb_mhc_ii_human_mimic_seqs.fasta")
OUT_NETMHCPAN = os.path.join(OUT_DIR, "iedb_mhc_ii_pairs_for_netmhcpan.csv.gz")

# Minimum epitope length for Class II
MIN_EPITOPE_LEN = 15

# Minimum human mimic length after padding (NetMHCIIpan default)
MIN_MIMIC_LEN = 15

# Mismatch range
MIN_MISMATCHES = 1  # exclude identical (0 mismatches)

# Minimum alignment coverage (alignment_length / query_length)
MIN_COVERAGE = 0.80

# Proteins to exclude (immunoglobulins, TCRs, MHC molecules)
EXCLUDE_PATTERNS = [
    r"HLA-",
    r"histocompatibility",
    r"immunoglobulin",
    r"T.cell.receptor",
    r"TCR\b",
    r"MHC\b",
]

# ====================================================
# Step 1: Load BLAST results
# ====================================================
def load_blast_results(filepath):
    print("=" * 60)
    print("Step 1: Loading BLAST results")
    print("=" * 60)

    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    df = pd.read_csv(filepath, sep="\t", names=cols)
    print(f"  Raw BLAST hits: {len(df):,}")
    return df


# ====================================================
# Step 2: Parse query header fields
# ====================================================
def parse_query_headers(df):
    """Parse FASTA header: structure_id|sequence|organism|mhc_allele"""
    print(f"\n{'=' * 60}")
    print("Step 2: Parsing query headers")
    print("=" * 60)

    parts = df["qseqid"].str.split("|", expand=True)

    df["structure_id"] = parts[0]
    df["viral_sequence"] = parts[1] if 1 in parts.columns else None
    df["source_organism"] = parts[2].str.replace("_", " ") if 2 in parts.columns else None
    df["mhc_allele"] = parts[3].str.replace("_", " ") if 3 in parts.columns else None

    df["hu_prot_id"] = df["sseqid"].apply(
        lambda x: x.split("|")[1] if "|" in str(x) else str(x)
    )
    df["hu_prot_name"] = df["sseqid"].apply(
        lambda x: x.split("|")[2] if str(x).count("|") >= 2 else ""
    )

    print(f"  Unique viral epitopes: {df['viral_sequence'].nunique():,}")
    print(f"  Unique MHC alleles: {df['mhc_allele'].nunique():,}")
    print(f"  Unique human proteins hit: {df['hu_prot_id'].nunique():,}")

    return df


# ====================================================
# Step 3: Apply filters
# ====================================================
def apply_filters(df):
    print(f"\n{'=' * 60}")
    print("Step 3: Applying filters")
    print("=" * 60)

    n_start = len(df)

    # --- 3a: Epitope length filter ---
    df = df[df["qlen"] >= MIN_EPITOPE_LEN]
    print(f"  After epitope length >= {MIN_EPITOPE_LEN}: {len(df):,} "
          f"(removed {n_start - len(df):,})")
    n_prev = len(df)

    # --- 3b: Coverage filter ---
    df = df.copy()
    df["coverage"] = df["length"] / df["qlen"]
    df = df[df["coverage"] >= MIN_COVERAGE]
    print(f"  After coverage >= {MIN_COVERAGE:.0%}: {len(df):,} "
          f"(removed {n_prev - len(df):,})")
    n_prev = len(df)

    # --- 3c: Mismatch filter ---
    df["n_identical"] = (df["pident"] * df["length"] / 100).round().astype(int)
    df["n_mismatches"] = df["length"] - df["n_identical"]
    df["max_mismatches"] = df["qlen"] // 3

    df = df[
        (df["n_mismatches"] >= MIN_MISMATCHES) &
        (df["n_mismatches"] <= df["max_mismatches"])
    ]
    print(f"  After {MIN_MISMATCHES} <= mismatches <= qlen//3: {len(df):,} "
          f"(removed {n_prev - len(df):,})")
    n_prev = len(df)

    # --- 3d: Remove deprecated UniProt IDs ---
    if os.path.exists(UNIPROT_FILTER):
        deprecated = pd.read_csv(UNIPROT_FILTER, sep="\t", names=["hu_prot_id"])
        df = df[~df["hu_prot_id"].isin(deprecated["hu_prot_id"])]
        print(f"  After removing deprecated UniProt: {len(df):,} "
              f"(removed {n_prev - len(df):,})")
        n_prev = len(df)
    else:
        print(f"  ⚠️ UniProt filter file not found, skipping")

    # --- 3e: Remove immune system proteins ---
    exclude_regex = "|".join(EXCLUDE_PATTERNS)
    immune_mask = df["hu_prot_name"].str.contains(
        exclude_regex, case=False, regex=True, na=False
    )
    n_immune = immune_mask.sum()
    df = df[~immune_mask]
    print(f"  After removing immune proteins (IG/TCR/MHC): {len(df):,} "
          f"(removed {n_immune:,})")

    # --- Summary ---
    print(f"\n  Mismatch distribution in filtered hits:")
    print(df["n_mismatches"].value_counts().sort_index().to_string())

    print(f"\n  Epitope length distribution:")
    print(df["qlen"].value_counts().sort_index().to_string())

    print(f"\n  Filtered results:")
    print(f"    Total hits: {len(df):,}")
    print(f"    Unique viral epitopes with mimics: {df['viral_sequence'].nunique():,}")
    print(f"    Unique human proteins: {df['hu_prot_id'].nunique():,}")
    print(f"    Unique viral-human pairs: "
          f"{df.drop_duplicates(['viral_sequence', 'hu_prot_id']).shape[0]:,}")

    return df


# ====================================================
# Step 4: Extract human mimic subsequences (with padding)
# ====================================================
def extract_human_sequences(df, human_fasta):
    print(f"\n{'=' * 60}")
    print("Step 4: Extracting human mimic subsequences")
    print(f"  Padding: human mimics padded to max(viral_len, {MIN_MIMIC_LEN})")
    print("=" * 60)

    # Load human proteome
    print("  Loading human proteome...")
    seq_dict = {}
    for rec in SeqIO.parse(human_fasta, "fasta"):
        acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
        seq_dict[acc] = rec
    print(f"  Loaded {len(seq_dict):,} human proteins")

    # Extract and pad subsequences
    records = []
    human_seqs = []
    padded_starts = []
    padded_ends = []
    not_found = 0
    n_padded = 0

    for _, row in df.iterrows():
        prot_id = row["hu_prot_id"]

        if prot_id not in seq_dict:
            not_found += 1
            human_seqs.append(None)
            padded_starts.append(None)
            padded_ends.append(None)
            continue

        full_seq = seq_dict[prot_id].seq
        prot_len = len(full_seq)
        start, end = int(row["sstart"]), int(row["send"])

        # Ensure start <= end
        if start > end:
            start, end = end, start

        aligned_len = end - start + 1
        viral_len = len(str(row["viral_sequence"]))
        target_len = max(viral_len, MIN_MIMIC_LEN)

        # Pad symmetrically if needed
        if aligned_len < target_len:
            pad_needed = target_len - aligned_len
            pad_left = pad_needed // 2
            pad_right = pad_needed - pad_left

            new_start = start - pad_left
            new_end = end + pad_right

            # Clamp to protein boundaries and rebalance
            if new_start < 1:
                new_end = min(prot_len, new_end + (1 - new_start))
                new_start = 1
            if new_end > prot_len:
                new_start = max(1, new_start - (new_end - prot_len))
                new_end = prot_len

            start, end = new_start, new_end
            n_padded += 1

        subseq = str(full_seq[start - 1:end])
        human_seqs.append(subseq)
        padded_starts.append(start)
        padded_ends.append(end)

        # FASTA record
        header = (
            f"{row['structure_id']}"
            f"|VIRAL_{row['viral_sequence']}"
            f"|HUMAN_{prot_id}"
            f"|{start}-{end}"
            f"|mm={row['n_mismatches']}"
            f"|{row['mhc_allele']}"
        )
        records.append(SeqRecord(Seq(subseq), id=header, description=""))

    df["human_mimic_sequence"] = human_seqs
    df["padded_sstart"] = padded_starts
    df["padded_send"] = padded_ends

    if not_found > 0:
        print(f"  ⚠️ {not_found} proteins not found in FASTA")

    # Remove rows where we couldn't extract the sequence
    df = df[df["human_mimic_sequence"].notna()]

    # Length diagnostics
    mimic_lens = df["human_mimic_sequence"].str.len()
    print(f"\n  Padding summary:")
    print(f"    Sequences padded: {n_padded:,} / {len(df):,}")
    print(f"    Human mimic length range: {mimic_lens.min()}-{mimic_lens.max()}")
    print(f"    Mimics < {MIN_MIMIC_LEN} aa after padding: "
          f"{(mimic_lens < MIN_MIMIC_LEN).sum()}")
    print(f"\n  Human mimic length distribution:")
    print(mimic_lens.value_counts().sort_index().to_string())

    # Write FASTA
    SeqIO.write(records, OUT_FASTA, "fasta")
    print(f"\n  ✅ Extracted {len(records):,} human mimic sequences → {OUT_FASTA}")

    return df


# ====================================================
# Step 5: Prepare NetMHCIIpan input
# ====================================================
def prepare_netmhcpan(df):
    """Create a table of viral-human pairs with MHC alleles for NetMHCIIpan."""
    print(f"\n{'=' * 60}")
    print("Step 5: Preparing NetMHCIIpan input")
    print("=" * 60)

    netmhcpan_cols = [
        "structure_id", "viral_sequence", "human_mimic_sequence",
        "hu_prot_id", "hu_prot_name", "source_organism", "mhc_allele",
        "n_mismatches", "n_identical", "pident", "coverage",
        "sstart", "send", "padded_sstart", "padded_send", "qlen"
    ]
    netmhcpan_df = df[[c for c in netmhcpan_cols if c in df.columns]].copy()

    # Classify MHC allele resolution
    # Matches any HLA allele with *XX:XX pattern
    four_digit = netmhcpan_df["mhc_allele"].str.contains(
        r"\*\d{2}:\d{2}", regex=True, na=False
    )
    netmhcpan_df["mhc_4digit"] = four_digit

    n_4digit = four_digit.sum()
    n_lowres = (~four_digit).sum()
    print(f"  Pairs with 4-digit HLA: {n_4digit:,}")
    print(f"  Pairs with low-res HLA: {n_lowres:,}")

    if n_lowres > 0:
        print(f"\n  Low-resolution HLA values (top 10):")
        lowres = netmhcpan_df[~four_digit]["mhc_allele"].value_counts().head(10)
        for allele, count in lowres.items():
            print(f"    {allele}: {count}")

    # Save full table
    netmhcpan_df.to_csv(OUT_NETMHCPAN, index=False, compression="gzip")
    print(f"\n  ✅ NetMHCIIpan input table: {len(netmhcpan_df):,} pairs → {OUT_NETMHCPAN}")

    # Also save 4-digit-only subset
    if n_4digit > 0:
        four_digit_path = os.path.join(OUT_DIR, "iedb_mhcii_pairs_4digit_hla.csv.gz")
        netmhcpan_df[four_digit].to_csv(four_digit_path, index=False, compression="gzip")
        print(f"  ✅ 4-digit HLA subset: {n_4digit:,} pairs → {four_digit_path}")

    return netmhcpan_df


# ====================================================
# Step 6: Summary
# ====================================================
def print_summary(df, netmhcpan_df):
    print(f"\n{'=' * 60}")
    print("Final Summary")
    print("=" * 60)

    print(f"  Total filtered viral-human pairs: {len(df):,}")
    print(f"  Unique viral epitopes: {df['viral_sequence'].nunique():,}")
    print(f"  Unique human mimic sequences: {df['human_mimic_sequence'].nunique():,}")
    print(f"  Unique human source proteins: {df['hu_prot_id'].nunique():,}")

    print(f"\n  By mismatch count:")
    for mm, group in df.groupby("n_mismatches"):
        print(f"    {mm} mismatches: {len(group):,} pairs, "
              f"{group['viral_sequence'].nunique():,} epitopes")

    print(f"\n  Top 10 source organisms:")
    orgs = df["source_organism"].value_counts().head(10)
    for org, count in orgs.items():
        print(f"    {org}: {count}")

    print(f"\n  Top 10 MHC alleles:")
    alleles = df["mhc_allele"].value_counts().head(10)
    for allele, count in alleles.items():
        print(f"    {allele}: {count}")

    if "mhc_4digit" in netmhcpan_df.columns:
        ready = netmhcpan_df["mhc_4digit"].sum()
        print(f"\n  Ready for NetMHCIIpan (4-digit HLA): {ready:,} pairs")

    print(f"\n  Output files:")
    print(f"    Filtered hits: {OUT_CSV}")
    print(f"    Human mimic FASTA: {OUT_FASTA}")
    print(f"    NetMHCIIpan input: {OUT_NETMHCPAN}")

    print(f"\n  Next steps:")
    print(f"    1. Run NetMHCIIpan on human mimics using paired HLA alleles")
    print(f"    2. Compute ΔBA (viral vs human binding affinity)")
    print(f"    3. Categorize: viral-dominant / human-dominant / equivalent")
    print(f"    4. Cross-reference with TCR data for DecoderTCR analysis")


# ====================================================
# Main
# ====================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Filter IEDB MHC-II BLAST hits")
    print("Extract human molecular mimicry candidates")
    print("  (with flanking-residue padding for NetMHCIIpan)")
    print("=" * 60)

    # Step 1: Load
    df = load_blast_results(BLAST_FILE)

    # Step 2: Parse headers
    df = parse_query_headers(df)

    # Step 3: Filter
    df = apply_filters(df)

    if df.empty:
        print("\n❌ No hits survived filtering.")
        exit(1)

    # Step 4: Extract human sequences (with padding)
    df = extract_human_sequences(df, HUMAN_FASTA)

    # Step 5: Save filtered hits
    print(f"\n  ✅ Filtered hits saved → {OUT_CSV}")
    df.to_csv(OUT_CSV, index=False, compression="gzip")

    # Step 6: Prepare NetMHCIIpan input
    netmhcpan_df = prepare_netmhcpan(df)

    # Step 7: Summary
    print_summary(df, netmhcpan_df)

    print(f"\n✅ Pipeline complete.")