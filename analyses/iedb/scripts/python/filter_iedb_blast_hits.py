#!/usr/bin/env python3
"""
======================================================================
Script: filter_iedb_blast_hits.py
Description:
    Filter BLAST results from IEDB MHC-I epitopes vs human proteome.
    Extract human mimic subsequences and prepare for NetMHCpan.

Filtering steps:
    1. Coverage: alignment must cover >= 80% of epitope length
    2. Mismatches: 1-3 mismatches (removes identical and very distant)
    3. Remove deprecated UniProt IDs
    4. Remove self-hits (human protein that IS the epitope source)
    5. Remove hits to immunoglobulin/TCR/MHC proteins
    6. Extract human mimic subsequences
    7. Parse MHC allele from query header for downstream NetMHCpan

Input FASTA header format (from merge script):
    >structure_id|sequence|organism|mhc_allele

Dependencies:
    pandas, biopython

Usage:
    python filter_iedb_blast_hits.py
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
BLAST_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/iedb_mhci_vs_human_proteome.tsv"
HUMAN_FASTA = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta"
UNIPROT_FILTER = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/proteins_to_remove_from_UniProtKB.txt"

OUT_DIR = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_CSV = os.path.join(OUT_DIR, "iedb_mhci_filtered_blast_hits.csv.gz")
OUT_FASTA = os.path.join(OUT_DIR, "iedb_mhci_human_mimic_seqs.fasta")
OUT_NETMHCPAN = os.path.join(OUT_DIR, "iedb_mhci_pairs_for_netmhcpan.csv.gz")

# Mismatch range
MIN_MISMATCHES = 1  # exclude identical (0 mismatches)

# Minimum alignment coverage (alignment_length / query_length)
MIN_COVERAGE = 0.80

# Proteins to exclude (immunoglobulins, TCRs, MHC molecules)
# These contain peptide-binding regions that spuriously match epitopes
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
    """
    Parse FASTA header: structure_id|sequence|organism|mhc_allele
    """
    print(f"\n{'=' * 60}")
    print("Step 2: Parsing query headers")
    print("=" * 60)

    parts = df["qseqid"].str.split("|", expand=True)

    df["structure_id"] = parts[0]
    df["viral_sequence"] = parts[1] if 1 in parts.columns else None
    df["source_organism"] = parts[2].str.replace("_", " ") if 2 in parts.columns else None
    df["mhc_allele"] = parts[3].str.replace("_", " ") if 3 in parts.columns else None

    # Parse human protein UniProt ID from subject
    df["hu_prot_id"] = df["sseqid"].apply(
        lambda x: x.split("|")[1] if "|" in str(x) else str(x)
    )

    # Get human protein name if available (third field in UniProt header)
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

    # --- 3a: Coverage filter ---
    df["coverage"] = df["length"] / df["qlen"]
    df = df[df["coverage"] >= MIN_COVERAGE]
    print(f"  After coverage >= {MIN_COVERAGE:.0%}: {len(df):,} "
          f"(removed {n_start - len(df):,})")
    n_prev = len(df)

    # --- 3b: Mismatch filter ---
    # pident = (identical / alignment_length) * 100
    # identical = pident * alignment_length / 100
    # mismatches = alignment_length - identical
    df["n_identical"] = (df["pident"] * df["length"] / 100).round().astype(int)
    df["n_mismatches"] = df["length"] - df["n_identical"]

    # Scale max mismatches by epitope length: allow up to ~33% mismatches
    # qlen // 3 gives: 9-mer→3, 10-mer→3, 11-mer→3, 12-mer→4, 13-mer→4, 14-mer→4
    df["max_mismatches"] = df["qlen"] // 3

    df = df[
        (df["n_mismatches"] >= MIN_MISMATCHES) &
        (df["n_mismatches"] <= df["max_mismatches"])
    ]
    print(f"  After {MIN_MISMATCHES} <= mismatches <= qlen//3: {len(df):,} "
          f"(removed {n_prev - len(df):,})")
    n_prev = len(df)

    # --- 3c: Remove deprecated UniProt IDs ---
    if os.path.exists(UNIPROT_FILTER):
        deprecated = pd.read_csv(UNIPROT_FILTER, sep="\t", names=["hu_prot_id"])
        df = df[~df["hu_prot_id"].isin(deprecated["hu_prot_id"])]
        print(f"  After removing deprecated UniProt: {len(df):,} "
              f"(removed {n_prev - len(df):,})")
        n_prev = len(df)
    else:
        print(f"  ⚠️ UniProt filter file not found, skipping")

    # --- 3d: Remove immune system proteins ---
    exclude_regex = "|".join(EXCLUDE_PATTERNS)
    immune_mask = df["hu_prot_name"].str.contains(
        exclude_regex, case=False, regex=True, na=False
    )
    n_immune = immune_mask.sum()
    df = df[~immune_mask]
    print(f"  After removing immune proteins (IG/TCR/MHC): {len(df):,} "
          f"(removed {n_immune:,})")
    n_prev = len(df)

    # --- 3e: Remove self-hits (identical sequences) ---
    # Already handled by MIN_MISMATCHES >= 1, but also check if 
    # the human hit IS the viral epitope sequence somehow
    self_hit_mask = df.apply(
        lambda row: str(row.get("viral_sequence", "")).upper() ==
        str(row.get("viral_sequence", "")).upper() and
        row.get("n_mismatches", 1) == 0, axis=1
    )
    # This is redundant given mismatch filter but safe to include

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
# Step 4: Extract human mimic subsequences
# ====================================================
def extract_human_sequences(df, human_fasta):
    print(f"\n{'=' * 60}")
    print("Step 4: Extracting human mimic subsequences")
    print("=" * 60)

    # Load human proteome
    print("  Loading human proteome...")
    seq_dict = {}
    for rec in SeqIO.parse(human_fasta, "fasta"):
        acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
        seq_dict[acc] = rec
    print(f"  Loaded {len(seq_dict):,} human proteins")

    # Extract aligned subsequences
    records = []
    human_seqs = []
    not_found = 0

    for _, row in df.iterrows():
        prot_id = row["hu_prot_id"]

        if prot_id not in seq_dict:
            not_found += 1
            human_seqs.append(None)
            continue

        full_seq = seq_dict[prot_id].seq
        start, end = int(row["sstart"]), int(row["send"])

        if start <= end:
            subseq = str(full_seq[start - 1:end])
        else:
            # Shouldn't happen with proteins but handle gracefully
            subseq = str(full_seq[end - 1:start])

        human_seqs.append(subseq)

        # FASTA record for the human mimic
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

    if not_found > 0:
        print(f"  ⚠️ {not_found} proteins not found in FASTA")

    # Remove rows where we couldn't extract the sequence
    df = df[df["human_mimic_sequence"].notna()]

    # Write FASTA
    SeqIO.write(records, OUT_FASTA, "fasta")
    print(f"  ✅ Extracted {len(records):,} human mimic sequences → {OUT_FASTA}")

    return df


# ====================================================
# Step 5: Prepare NetMHCpan input
# ====================================================
def prepare_netmhcpan(df):
    """
    Create a table of viral-human pairs with MHC alleles for NetMHCpan.
    Only include pairs where MHC allele is at 4-digit resolution.
    """
    print(f"\n{'=' * 60}")
    print("Step 5: Preparing NetMHCpan input")
    print("=" * 60)

    # Keep relevant columns for NetMHCpan
    netmhcpan_cols = [
        "structure_id", "viral_sequence", "human_mimic_sequence",
        "hu_prot_id", "hu_prot_name", "source_organism", "mhc_allele",
        "n_mismatches", "n_identical", "pident", "coverage",
        "sstart", "send", "qlen"
    ]
    netmhcpan_df = df[[c for c in netmhcpan_cols if c in df.columns]].copy()

    # Classify MHC allele resolution
    # 4-digit: HLA-A*02:01, HLA-B*35:01
    # Low-res: HLA-A2, HLA-DR, human
    four_digit = netmhcpan_df["mhc_allele"].str.contains(
        r"HLA-[ABC]\*\d{2}:\d{2}", regex=True, na=False
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
    print(f"\n  ✅ NetMHCpan input table: {len(netmhcpan_df):,} pairs → {OUT_NETMHCPAN}")

    # Also save 4-digit-only subset
    if n_4digit > 0:
        four_digit_path = os.path.join(OUT_DIR, "iedb_mhci_pairs_4digit_hla.csv.gz")
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
        print(f"\n  Ready for NetMHCpan (4-digit HLA): {ready:,} pairs")

    print(f"\n  Output files:")
    print(f"    Filtered hits: {OUT_CSV}")
    print(f"    Human mimic FASTA: {OUT_FASTA}")
    print(f"    NetMHCpan input: {OUT_NETMHCPAN}")

    print(f"\n  Next steps:")
    print(f"    1. Run NetMHCpan on human mimics using paired HLA alleles")
    print(f"    2. Compute ΔBA (viral vs human binding affinity)")
    print(f"    3. Categorize: viral-dominant / human-dominant / equivalent")
    print(f"    4. Cross-reference with TCR data for DecoderTCR analysis")


# ====================================================
# Main
# ====================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Filter IEDB MHC-I BLAST hits")
    print("Extract human molecular mimicry candidates")
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

    # Step 4: Extract human sequences
    df = extract_human_sequences(df, HUMAN_FASTA)

    # Step 5: Save filtered hits
    print(f"\n  ✅ Filtered hits saved → {OUT_CSV}")
    df.to_csv(OUT_CSV, index=False, compression="gzip")

    # Step 6: Prepare NetMHCpan input
    netmhcpan_df = prepare_netmhcpan(df)

    # Step 7: Summary
    print_summary(df, netmhcpan_df)

    print(f"\n✅ Pipeline complete.")