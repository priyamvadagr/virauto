#!/usr/bin/env python3
"""
======================================================================
Script: filter_iedb_blast_hits.py
Description:
    Filter BLAST results from IEDB MHC-I or MHC-II epitopes vs human proteome.
    Extract human mimic subsequences and prepare for NetMHCpan.

    For MHC-II: extends short human mimics to 15 aa by adding symmetric
    flanking residues from the human protein sequence. Drops pairs where
    the viral epitope is < 15 aa (NetMHCIIpan minimum).

Input FASTA header format:
    >structure_id|sequence|organism|mhc_allele

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
import argparse

# ====================================================
# Config
# ====================================================
UNIPROT_FILTER = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/proteins_to_remove_from_UniProtKB.txt"

MIN_MISMATCHES = 1
MIN_COVERAGE = 0.80
MIN_PEPTIDE_LEN_CLASS_II = 15  # NetMHCIIpan minimum

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

    df["coverage"] = df["length"] / df["qlen"]
    df = df[df["coverage"] >= MIN_COVERAGE].copy()
    print(f"  After coverage >= {MIN_COVERAGE:.0%}: {len(df):,} "
          f"(removed {n_start - len(df):,})")
    n_prev = len(df)

    df["n_identical"] = (df["pident"] * df["length"] / 100).round().astype(int)
    df["n_mismatches"] = df["length"] - df["n_identical"]
    df["max_mismatches"] = df["qlen"] // 3

    # Report exact matches before filtering them out
    n_exact = (df["n_mismatches"] == 0).sum()
    n_exact_epitopes = df[df["n_mismatches"] == 0]["viral_sequence"].nunique() if n_exact > 0 else 0
    print(f"  Exact matches (0 mismatches): {n_exact:,} hits "
          f"({n_exact_epitopes:,} unique epitopes) — will be excluded")

    df = df[
        (df["n_mismatches"] >= MIN_MISMATCHES) &
        (df["n_mismatches"] <= df["max_mismatches"])
    ].copy()
    print(f"  After {MIN_MISMATCHES} <= mismatches <= qlen//3: {len(df):,} "
          f"(removed {n_prev - len(df):,})")
    n_prev = len(df)

    if os.path.exists(UNIPROT_FILTER):
        deprecated = pd.read_csv(UNIPROT_FILTER, sep="\t", names=["hu_prot_id"])
        df = df[~df["hu_prot_id"].isin(deprecated["hu_prot_id"])].copy()
        print(f"  After removing deprecated UniProt: {len(df):,} "
              f"(removed {n_prev - len(df):,})")
        n_prev = len(df)

    exclude_regex = "|".join(EXCLUDE_PATTERNS)
    immune_mask = df["hu_prot_name"].str.contains(
        exclude_regex, case=False, regex=True, na=False
    )
    df = df[~immune_mask].copy()
    print(f"  After removing immune proteins: {len(df):,} "
          f"(removed {immune_mask.sum():,})")

    print(f"\n  Mismatch distribution:")
    print(df["n_mismatches"].value_counts().sort_index().to_string())
    print(f"\n  Epitope length distribution:")
    print(df["qlen"].value_counts().sort_index().to_string())

    return df


# ====================================================
# Step 4: Extract human mimic subsequences
#         with flanking extension for Class II
# ====================================================
def extract_human_sequences(df, human_fasta, mhc_class="I"):
    print(f"\n{'=' * 60}")
    print("Step 4: Extracting human mimic subsequences")
    if mhc_class == "II":
        print(f"  Class II mode: extending short mimics to {MIN_PEPTIDE_LEN_CLASS_II} aa")
    print("=" * 60)

    print("  Loading human proteome...")
    seq_dict = {}
    for rec in SeqIO.parse(human_fasta, "fasta"):
        acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
        seq_dict[acc] = str(rec.seq)
    print(f"  Loaded {len(seq_dict):,} human proteins")

    records = []
    human_seqs = []
    new_starts = []
    new_ends = []
    n_extended = 0
    n_too_short_viral = 0
    n_too_short_human = 0
    not_found = 0

    for _, row in df.iterrows():
        prot_id = row["hu_prot_id"]

        if prot_id not in seq_dict:
            not_found += 1
            human_seqs.append(None)
            new_starts.append(None)
            new_ends.append(None)
            continue

        full_seq = seq_dict[prot_id]
        prot_len = len(full_seq)
        start, end = int(row["sstart"]), int(row["send"])

        # Ensure start <= end
        if start > end:
            start, end = end, start

        subseq = full_seq[start - 1:end]
        cur_len = len(subseq)
        new_start = start
        new_end = end

        # --- Class II: extend short human mimics to 15 aa ---
        if mhc_class == "II":
            # Skip if viral epitope is < 15 aa
            viral_len = len(str(row.get("viral_sequence", "")))
            if viral_len < MIN_PEPTIDE_LEN_CLASS_II:
                n_too_short_viral += 1
                human_seqs.append(None)
                new_starts.append(None)
                new_ends.append(None)
                continue

            # Extend human mimic if < 15 aa
            if cur_len < MIN_PEPTIDE_LEN_CLASS_II:
                needed = MIN_PEPTIDE_LEN_CLASS_II - cur_len

                # Symmetric extension: add equally on both sides
                extend_left = needed // 2
                extend_right = needed - extend_left

                # Compute new boundaries (1-based)
                new_start = max(1, start - extend_left)
                new_end = min(prot_len, end + extend_right)

                # If one side hit the protein boundary, add more to the other
                actual_left = start - new_start
                actual_right = new_end - end
                total_added = actual_left + actual_right

                if total_added < needed:
                    remaining = needed - total_added
                    if new_start > 1:
                        new_start = max(1, new_start - remaining)
                    elif new_end < prot_len:
                        new_end = min(prot_len, new_end + remaining)

                subseq = full_seq[new_start - 1:new_end]
                n_extended += 1

                # If still < 15 after extension (very short protein), skip
                if len(subseq) < MIN_PEPTIDE_LEN_CLASS_II:
                    n_too_short_human += 1
                    human_seqs.append(None)
                    new_starts.append(None)
                    new_ends.append(None)
                    continue

        human_seqs.append(subseq)
        new_starts.append(new_start)
        new_ends.append(new_end)

        header = (
            f"{row['structure_id']}"
            f"|VIRAL_{row['viral_sequence']}"
            f"|HUMAN_{prot_id}"
            f"|{new_start}-{new_end}"
            f"|mm={row['n_mismatches']}"
            f"|{row['mhc_allele']}"
        )
        records.append(SeqRecord(Seq(subseq), id=header, description=""))

    df["human_mimic_sequence"] = human_seqs
    df["sstart"] = new_starts
    df["send"] = new_ends

    if not_found > 0:
        print(f"  ⚠️ {not_found} proteins not found in FASTA")

    if mhc_class == "II":
        print(f"  Extended to {MIN_PEPTIDE_LEN_CLASS_II} aa: {n_extended:,} human mimics")
        print(f"  Dropped (viral < {MIN_PEPTIDE_LEN_CLASS_II} aa): {n_too_short_viral:,}")
        print(f"  Dropped (human still < {MIN_PEPTIDE_LEN_CLASS_II} aa after extension): {n_too_short_human:,}")

    df = df[df["human_mimic_sequence"].notna()].copy()

    out_fasta = os.path.join(OUT_DIR, f"iedb_mhc_{'ii' if mhc_class == 'II' else 'i'}_human_mimic_seqs.fasta")
    SeqIO.write(records, out_fasta, "fasta")
    print(f"  ✅ Extracted {len(records):,} human mimic sequences → {out_fasta}")

    return df


# ====================================================
# Step 5: Prepare NetMHCpan input
# ====================================================
def prepare_netmhcpan(df, mhc_class="I"):
    print(f"\n{'=' * 60}")
    print("Step 5: Preparing NetMHCpan input")
    print("=" * 60)

    netmhcpan_cols = [
        "structure_id", "viral_sequence", "human_mimic_sequence",
        "hu_prot_id", "hu_prot_name", "source_organism", "mhc_allele",
        "n_mismatches", "n_identical", "pident", "coverage",
        "sstart", "send", "qlen"
    ]
    netmhcpan_df = df[[c for c in netmhcpan_cols if c in df.columns]].copy()

    # Classify MHC allele resolution
    if mhc_class == "I":
        four_digit = netmhcpan_df["mhc_allele"].str.contains(
            r"HLA-[ABC]\*\d{2}:\d{2}", regex=True, na=False
        )
    else:
        # Class II: multiple formats are valid
        four_digit = netmhcpan_df["mhc_allele"].str.contains(
            r"HLA-D[RQP]", regex=True, na=False
        )

    netmhcpan_df["mhc_4digit"] = four_digit

    n_4digit = four_digit.sum()
    n_lowres = (~four_digit).sum()
    print(f"  Pairs with resolved HLA: {n_4digit:,}")
    print(f"  Pairs with low-res/invalid HLA: {n_lowres:,}")

    if n_lowres > 0:
        print(f"\n  Low-resolution HLA values (top 10):")
        lowres = netmhcpan_df[~four_digit]["mhc_allele"].value_counts().head(10)
        for allele, count in lowres.items():
            print(f"    {allele}: {count}")

    # Save full table
    prefix = "mhc_ii" if mhc_class == "II" else "mhci"
    out_full = os.path.join(OUT_DIR, f"iedb_{prefix}_pairs_for_netmhcpan.csv.gz")
    netmhcpan_df.to_csv(out_full, index=False, compression="gzip")
    print(f"\n  ✅ Full table: {len(netmhcpan_df):,} pairs → {out_full}")

    if n_4digit > 0:
        out_4digit = os.path.join(OUT_DIR, f"iedb_{prefix}_pairs_4digit_hla.csv.gz")
        netmhcpan_df[four_digit].to_csv(out_4digit, index=False, compression="gzip")
        print(f"  ✅ Resolved HLA subset: {n_4digit:,} pairs → {out_4digit}")

    # Class II: report peptide length distribution after extension
    if mhc_class == "II":
        print(f"\n  Human mimic length distribution after extension:")
        hu_lens = netmhcpan_df["human_mimic_sequence"].str.len()
        print(hu_lens.describe().to_string())

    return netmhcpan_df


# ====================================================
# Summary
# ====================================================
def print_summary(df, netmhcpan_df, mhc_class):
    print(f"\n{'=' * 60}")
    print(f"Final Summary — MHC Class {'II' if mhc_class == 'II' else 'I'}")
    print("=" * 60)

    print(f"  Total filtered pairs: {len(df):,}")
    print(f"  Unique viral epitopes: {df['viral_sequence'].nunique():,}")
    print(f"  Unique human mimic sequences: {df['human_mimic_sequence'].nunique():,}")
    print(f"  Unique human proteins: {df['hu_prot_id'].nunique():,}")

    print(f"\n  By mismatch count:")
    for mm, group in df.groupby("n_mismatches"):
        print(f"    {mm} mismatches: {len(group):,} pairs")

    print(f"\n  Top 10 source organisms:")
    for org, count in df["source_organism"].value_counts().head(10).items():
        print(f"    {org}: {count}")

    print(f"\n  Top 10 MHC alleles:")
    for allele, count in df["mhc_allele"].value_counts().head(10).items():
        print(f"    {allele}: {count}")

    if "mhc_4digit" in netmhcpan_df.columns:
        ready = netmhcpan_df["mhc_4digit"].sum()
        print(f"\n  Ready for NetMHC{'II' if mhc_class == 'II' else ''}pan: {ready:,} pairs")


# ====================================================
# Main
# ====================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter IEDB BLAST hits for mimicry candidates"
    )
    parser.add_argument("--mhc-class", choices=["I", "II"], default="I",
                        help="MHC class (default: I)")
    parser.add_argument("--blast-file", required=True,
                        help="BLAST output TSV")
    parser.add_argument("--human-fasta", required=True,
                        help="Human proteome FASTA")
    parser.add_argument("--out-dir", required=True,
                        help="Output directory")
    args = parser.parse_args()

    OUT_DIR = args.out_dir
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 60)
    print(f"Filter IEDB MHC Class {'II' if args.mhc_class == 'II' else 'I'} BLAST hits")
    print("=" * 60)

    df = load_blast_results(args.blast_file)
    df = parse_query_headers(df)
    df = apply_filters(df)

    if df.empty:
        print("\n❌ No hits survived filtering.")
        exit(1)

    df = extract_human_sequences(df, args.human_fasta, mhc_class=args.mhc_class)

    out_csv = os.path.join(OUT_DIR, f"iedb_mhc_{'ii' if args.mhc_class == 'II' else 'i'}_filtered_blast_hits.csv.gz")
    df.to_csv(out_csv, index=False, compression="gzip")
    print(f"\n  ✅ Filtered hits saved → {out_csv}")

    netmhcpan_df = prepare_netmhcpan(df, mhc_class=args.mhc_class)
    print_summary(df, netmhcpan_df, args.mhc_class)

    print(f"\n✅ Pipeline complete.")