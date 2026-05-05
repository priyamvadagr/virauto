#!/usr/bin/env python3
"""
======================================================================
Script: fetch_mhc_i_sequences.py

Description:
    Build HLA_seq strings for DecoderTCR from locally downloaded
    IMGT/HLA protein FASTAs and B2M sequence.

    HLA_seq = mature_heavy_chain_aa + mature_b2m_aa

    Signal peptide stripped by finding conserved GSHSM motif.
    Resulting heavy chains are ~341 aa (1 residue shorter than
    DecoderTCR training data — deemed negligible for scoring).

Input files:
    /ix/djishnu/Priyamvada/virauto/data/refs/hla/A_prot.fasta
    /ix/djishnu/Priyamvada/virauto/data/refs/hla/B_prot.fasta
    /ix/djishnu/Priyamvada/virauto/data/refs/hla/C_prot.fasta
    /ix/djishnu/Priyamvada/virauto/data/refs/hla/b2m_human.fasta

Output:
    hla_sequences.tsv    columns: allele, source_allele, fallback,
                                  heavy_chain_len, b2m_len, hla_seq
    hla_sequences.fasta  for manual inspection
    hla_fallback_log.tsv alleles resolved by 2-digit fallback

Usage:
    python fetch_mhc_i_sequences.py
======================================================================
"""

import os
import re
import pandas as pd

# ====================================================================
# Config
# ====================================================================
MIMICRY_TCR_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/mimicry_candidates/mhc_i/tcr_mapping/mimicry_with_tcr.csv.gz"
HLA_DIR          = "/ix/djishnu/Priyamvada/virauto/data/refs/hla"
OUT_DIR          = '/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/mimicry_candidates/mhc_i/hla_sequences'

LOCUS_FASTAS = {
    "A": os.path.join(HLA_DIR, "A_prot.fasta"),
    "B": os.path.join(HLA_DIR, "B_prot.fasta"),
    "C": os.path.join(HLA_DIR, "C_prot.fasta"),
}
B2M_FASTA    = os.path.join(HLA_DIR, "b2m_human.fasta")
HLA_TSV      = os.path.join(OUT_DIR, "hla_sequences.tsv")
HLA_FASTA    = os.path.join(OUT_DIR, "hla_sequences.fasta")
FALLBACK_LOG = os.path.join(OUT_DIR, "hla_fallback_log.tsv")

B2M_SIGNAL_LEN = 20  # residues 1-20 of UniProt P61769 are signal peptide

# ====================================================================
# Step 1: Parse B2M
# ====================================================================
print("Parsing B2M (UniProt P61769)...")
b2m_full = ""
with open(B2M_FASTA) as fh:
    for line in fh:
        if not line.startswith(">"):
            b2m_full += line.strip()
b2m_mature = b2m_full[B2M_SIGNAL_LEN:]
print(f"  Full B2M   : {len(b2m_full)} aa")
print(f"  Mature B2M : {len(b2m_mature)} aa")

# ====================================================================
# Step 2: Parse IMGT/HLA protein FASTAs
# ====================================================================
def parse_imgt_fasta(fasta_path):
    """
    Parse an IMGT/HLA protein FASTA.
    Header: >HLA:HLA00001 A*01:01:01:01 365 bp
    Returns dict: 4-digit allele (e.g. 'A*01:01') → longest sequence.
    """
    allele_seqs    = {}
    current_allele = None
    current_seq    = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_allele and current_seq:
                    allele_seqs.setdefault(current_allele, []).append(
                        "".join(current_seq)
                    )
                parts = line[1:].split()
                if len(parts) >= 2:
                    fields         = parts[1].split(":")
                    current_allele = ":".join(fields[:2])
                else:
                    current_allele = None
                current_seq = []
            else:
                current_seq.append(line.strip())

    if current_allele and current_seq:
        allele_seqs.setdefault(current_allele, []).append("".join(current_seq))

    return {a: max(seqs, key=len) for a, seqs in allele_seqs.items()}


def strip_signal_peptide(seq):
    """
    Strip HLA heavy chain signal peptide by finding the conserved
    mature chain start motif (GSHSM or variants).
    Returns mature sequence only — no terminal residue appended.
    """
    for motif in ["GSHSM", "GSHSL", "GSHDM"]:
        pos = seq.find(motif)
        if pos > 0:
            return seq[pos:]
    return seq  # no signal peptide detected


print("\nParsing IMGT/HLA protein FASTAs...")
locus_seqs = {}
for locus, path in LOCUS_FASTAS.items():
    locus_seqs[locus] = parse_imgt_fasta(path)
    print(f"  Locus {locus}: {len(locus_seqs[locus]):,} alleles")

# ====================================================================
# Step 3: Get alleles needed from mimicry data
# ====================================================================
print("\nLoading mimicry TCR data...")
tcr = pd.read_csv(MIMICRY_TCR_FILE, low_memory=False)
both = tcr["alpha_cdr3"].notna() & tcr["beta_cdr3"].notna()
alleles_needed = sorted(tcr.loc[both, "mhc_allele"].dropna().unique())
print(f"  Alleles needed: {len(alleles_needed)}")
for a in alleles_needed:
    print(f"    {a}")

# ====================================================================
# Step 4: Look up sequences for each allele
# ====================================================================
def lookup_allele(allele, locus_seqs):
    """
    Parse allele string (e.g. 'HLA-A*02:01'), look up in locus_seqs.
    Returns (mature_sequence, source_allele, is_fallback).
    Falls back to *:01 of the same 2-digit group if exact match missing.
    """
    m = re.match(r"HLA-([ABC])\*(\d+):(\d+)", allele)
    if not m:
        return None, None, False

    locus  = m.group(1)
    field1 = m.group(2)
    field2 = m.group(3)
    imgt   = f"{locus}*{field1}:{field2}"

    db = locus_seqs.get(locus, {})

    # Exact match
    if imgt in db:
        return strip_signal_peptide(db[imgt]), imgt, False

    # 2-digit fallback
    group      = f"{locus}*{field1}"
    candidates = {k: v for k, v in db.items() if k.startswith(group + ":")}
    if candidates:
        preferred = f"{group}:01"
        source    = preferred if preferred in candidates else sorted(candidates)[0]
        return strip_signal_peptide(candidates[source]), source, True

    return None, None, False


print("\nLooking up HLA sequences...")
results   = []
fallbacks = []

for allele in alleles_needed:
    heavy_seq, source, is_fallback = lookup_allele(allele, locus_seqs)

    if heavy_seq is None:
        print(f"  ❌ {allele} — not found")
        continue

    hla_seq = heavy_seq + b2m_mature

    tag = "⚠️ " if is_fallback else "✅"
    fb  = f"→ fallback {source}" if is_fallback else f"source={source}"
    print(f"  {tag} {allele:22s}  {fb:20s}  "
          f"heavy={len(heavy_seq)} aa  total={len(hla_seq)} aa")

    if is_fallback:
        fallbacks.append({"requested_allele": allele, "fallback_allele": source})

    results.append({
        "allele"         : allele,
        "source_allele"  : source,
        "fallback"       : is_fallback,
        "heavy_chain_len": len(heavy_seq),
        "b2m_len"        : len(b2m_mature),
        "total_len"      : len(hla_seq),
        "hla_seq"        : hla_seq,
    })

# ====================================================================
# Save
# ====================================================================
results_df = pd.DataFrame(results)
results_df.to_csv(HLA_TSV, sep="\t", index=False)
print(f"\n  HLA TSV   → {HLA_TSV}  ({len(results_df)} alleles)")

with open(HLA_FASTA, "w") as fh:
    for _, row in results_df.iterrows():
        fh.write(f">{row['allele']} source={row['source_allele']} "
                 f"heavy={row['heavy_chain_len']} total={row['total_len']}\n")
        seq = row["hla_seq"]
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + "\n")
print(f"  HLA FASTA → {HLA_FASTA}")

if fallbacks:
    pd.DataFrame(fallbacks).to_csv(FALLBACK_LOG, sep="\t", index=False)
    print(f"  Fallback log → {FALLBACK_LOG}  ({len(fallbacks)} fallbacks)")

print(f"\n{'=' * 55}")
print(f"  Summary")
print(f"{'=' * 55}")
print(f"  Alleles requested : {len(alleles_needed)}")
print(f"  Alleles resolved  : {len(results_df)}")
print(f"  Exact matches     : {(~results_df['fallback']).sum()}")
print(f"  Fallbacks used    : {results_df['fallback'].sum()}")
print(f"\n  Heavy chain lengths (expect ~341 aa):")
print(results_df[["allele", "source_allele", "heavy_chain_len", "total_len"]]
      .to_string(index=False))
print(f"\n✅ Done.")