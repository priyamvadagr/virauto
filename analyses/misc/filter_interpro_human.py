#!/usr/bin/env python3
"""
filter_interpro_human.py

Filter protein2ipr.dat.gz to only human UniProt ACs.
Produces a much smaller file for downstream analyses.

Input:
    - protein2ipr.dat.gz (full InterPro, all organisms)
    - uniprot_human_all.fasta (human proteome)

Output:
    - protein2ipr_human.tsv.gz (human-only, with inferred source_db)

Usage:
    python filter_interpro_human.py
"""

import gzip
import os

# ====================================================
# Config
# ====================================================
PROTEOME_FASTA = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta"
INTERPRO_FILE = "/ix/djishnu/Priyamvada/virauto/data/refs/interpro/protein2ipr.dat.gz"
OUTPUT_FILE = "/ix/djishnu/Priyamvada/virauto/data/refs/interpro/protein2ipr_human.tsv.gz"

# ====================================================
# Step 1: Parse human UniProt ACs from FASTA
# ====================================================
print("Parsing human UniProt ACs from FASTA...")
human_acs = set()
with open(PROTEOME_FASTA) as fh:
    for line in fh:
        if line.startswith(">"):
            parts = line[1:].split("|")
            if len(parts) >= 3:
                human_acs.add(parts[1].strip())

print(f"  {len(human_acs):,} human UniProt ACs")

# ====================================================
# Step 2: Infer source DB from source ID prefix
# ====================================================
def infer_source_db(source_id):
    if source_id.startswith("PF"):
        return "Pfam"
    elif source_id.startswith("SM"):
        return "SMART"
    elif source_id.startswith("PS"):
        return "PROSITE"
    elif source_id.startswith("SSF"):
        return "SUPERFAMILY"
    elif source_id.startswith("G3DSA"):
        return "Gene3D"
    elif source_id.startswith("PTHR"):
        return "PANTHER"
    elif source_id.startswith("TIGR"):
        return "TIGRFAMs"
    elif source_id.startswith("cd"):
        return "CDD"
    elif source_id.startswith("MF"):
        return "HAMAP"
    elif source_id.startswith("PIRSF"):
        return "PIRSF"
    else:
        return source_id.split(":")[0] if ":" in source_id else source_id[:4]

# ====================================================
# Step 3: Stream InterPro, filter to human, write output
# ====================================================
print(f"\nFiltering {INTERPRO_FILE} to human proteins...")

n_total = 0
n_kept = 0
n_skipped_format = 0

with gzip.open(INTERPRO_FILE, "rt") as fin, \
     gzip.open(OUTPUT_FILE, "wt") as fout:

    # Write header
    fout.write("uniprot_ac\tipr_id\tipr_name\tsource_id\tsource_db\tstart\tend\n")

    for line in fin:
        n_total += 1

        if n_total % 50_000_000 == 0:
            print(f"  Processed {n_total:,} lines, kept {n_kept:,} ...")

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            n_skipped_format += 1
            continue

        ac = parts[0]
        if ac not in human_acs:
            continue

        ipr_id = parts[1]
        ipr_name = parts[2]
        source_id = parts[3]
        source_db = infer_source_db(source_id)

        try:
            start = int(parts[4])
            end = int(parts[5])
        except ValueError:
            n_skipped_format += 1
            continue

        fout.write(f"{ac}\t{ipr_id}\t{ipr_name}\t{source_id}\t{source_db}\t{start}\t{end}\n")
        n_kept += 1

print(f"\n  Total lines processed: {n_total:,}")
print(f"  Lines kept (human): {n_kept:,}")
print(f"  Lines skipped (format): {n_skipped_format:,}")

# ====================================================
# Step 4: Quick summary of output
# ====================================================
print(f"\nReading back output for summary...")
import pandas as pd
df = pd.read_csv(OUTPUT_FILE, sep="\t")
print(f"  Output rows: {len(df):,}")
print(f"  Unique proteins: {df['uniprot_ac'].nunique():,}")
print(f"  Unique domains (ipr_name): {df['ipr_name'].nunique():,}")
print(f"  Source DBs: {df['source_db'].value_counts().to_dict()}")
print(f"\n  Output: {OUTPUT_FILE}")
print(f"  Size: {os.path.getsize(OUTPUT_FILE) / 1e6:.1f} MB")

print(f"\n✅ Done.")