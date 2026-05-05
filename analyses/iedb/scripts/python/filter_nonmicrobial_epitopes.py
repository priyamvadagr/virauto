#!/usr/bin/env python3
"""
filter_nonmicrobial_epitopes.py

Filter IEDB epitope FASTA files to remove non-microbial organisms
(allergens, animal self-antigens, insect venoms, plants).
Keeps only bacteria, viruses, fungi, and parasites.

Applied at the earliest pipeline stage (pre-BLAST).

Input:
    IEDB epitope FASTA (Class I or Class II)

Output:
    Filtered FASTA (microbial only)

Usage:
    python filter_nonmicrobial_epitopes.py \
        --input iedb_epitopes_mhc_ii.fasta \
        --output iedb_epitopes_mhc_ii_microbial.fasta
"""

import argparse
from collections import Counter

# ====================================================
# Non-microbial organism patterns
# ====================================================
# These are matched case-insensitively against the organism name
# from the FASTA header (field 3, underscores replaced with spaces)

EXCLUDE_PATTERNS = [
    # --- Mammals (self-antigens, model organisms) ---
    "mus musculus",
    "homo sapiens",
    "bos taurus",
    "bovine",
    "canis lupus",
    "canis familiaris",
    "felis catus",
    "equus caballus",
    "rattus norvegicus",
    "ovis aries",
    "sus scrofa",
    "macaca",
    "pan troglodytes",
    "gallus gallus",
    # --- Plant allergens ---
    "phleum pratense",       # timothy grass
    "lolium perenne",        # ryegrass
    "cynodon dactylon",      # bermuda grass
    "poa pratensis",         # bluegrass
    "secale cereale",        # rye
    "triticum aestivum",     # wheat
    "triticum turgidum",     # durum wheat
    "arachis hypogaea",      # peanut
    "betula pendula",        # birch
    "betula verrucosa",      # birch
    "alnus glutinosa",       # alder
    "corylus avellana",      # hazelnut
    "cryptomeria japonica",  # cedar
    "chamaecyparis obtusa",  # cypress
    "cupressus",             # cypress
    "ambrosia artemisiifolia",  # ragweed
    "artemisia vulgaris",    # mugwort
    "plantago lanceolata",   # plantain
    "daucus carota",         # carrot
    "malus domestica",       # apple
    "prunus persica",        # peach
    "prunus avium",          # cherry
    "hevea brasiliensis",    # rubber tree
    "glycine max",           # soybean
    "brassica",              # mustard/cabbage family
    "olea europaea",         # olive
    "parietaria judaica",    # pellitory
    # --- Arthropod allergens (mites, cockroach) ---
    "dermatophagoides",      # dust mites
    "euroglyphus",           # dust mites
    "blomia",                # storage mites
    "blattella germanica",   # cockroach
    "periplaneta americana", # cockroach
    # --- Insect venoms ---
    "vespula",               # wasp
    "apis mellifera",        # honey bee
    "polistes",              # paper wasp
    "solenopsis",            # fire ant
    # --- Other non-microbial ---
    "anisakis",              # parasitic nematode in fish (sometimes kept)
    "caenorhabditis",        # model nematode
    "drosophila",            # fruit fly
]


def is_excluded(organism_name):
    """Check if organism matches any exclusion pattern."""
    name_lower = organism_name.lower()
    for pattern in EXCLUDE_PATTERNS:
        if pattern in name_lower:
            return True
    return False


def main():
    parser = argparse.ArgumentParser(
        description="Filter IEDB FASTA to microbial organisms only"
    )
    parser.add_argument("--input", required=True,
                        help="Input IEDB epitope FASTA")
    parser.add_argument("--output", required=True,
                        help="Output filtered FASTA")
    parser.add_argument("--report", default=None,
                        help="Optional: save exclusion report TSV")
    args = parser.parse_args()

    # --- Parse and filter ---
    print(f"Filtering {args.input} ...")

    kept = 0
    excluded = 0
    excluded_orgs = Counter()
    kept_orgs = Counter()

    with open(args.input) as fin, open(args.output, "w") as fout:
        write_current = False
        current_org = None

        for line in fin:
            if line.startswith(">"):
                parts = line[1:].strip().split("|")
                if len(parts) >= 3:
                    org = parts[2].replace("_", " ")
                else:
                    org = "unknown"

                if is_excluded(org):
                    write_current = False
                    excluded += 1
                    excluded_orgs[org] += 1
                else:
                    write_current = True
                    kept += 1
                    kept_orgs[org] += 1
                    fout.write(line)
            else:
                if write_current:
                    fout.write(line)

    total = kept + excluded
    print(f"\n  Total epitopes: {total:,}")
    print(f"  Kept (microbial): {kept:,} ({kept/total*100:.1f}%)")
    print(f"  Excluded (non-microbial): {excluded:,} ({excluded/total*100:.1f}%)")
    print(f"\n  Kept organisms: {len(kept_orgs)}")
    print(f"  Excluded organisms: {len(excluded_orgs)}")

    # Show excluded organisms
    print(f"\n  Excluded organisms by count:")
    for org, n in excluded_orgs.most_common():
        print(f"    {n:>5}  {org}")

    # Show top kept organisms
    print(f"\n  Top 15 kept organisms:")
    for org, n in kept_orgs.most_common(15):
        print(f"    {n:>5}  {org}")

    # Save report
    if args.report:
        import pandas as pd
        rows = []
        for org, n in kept_orgs.items():
            rows.append({"organism": org, "count": n, "status": "kept"})
        for org, n in excluded_orgs.items():
            rows.append({"organism": org, "count": n, "status": "excluded"})
        pd.DataFrame(rows).sort_values(["status", "count"], ascending=[True, False]).to_csv(
            args.report, sep="\t", index=False
        )
        print(f"\n  Report: {args.report}")

    print(f"\n  Output: {args.output}")
    print(f"\n✅ Done.")


if __name__ == "__main__":
    main()