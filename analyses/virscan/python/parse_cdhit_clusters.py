#!/usr/bin/env python3
"""
===============================================================================
Script: parse_cdhit_clusters.py
Description:
    Parses CD-HIT cluster output (.clstr) and creates:
      1. Mapping file: representative → equivalent alleles (with readable names)
      2. Representative allele list with readable names

Workflow:
    1. Load HLA allele list (internal ID → human-readable name).
    2. Parse CD-HIT .clstr file to extract cluster representatives and members.
    3. Write:
        - Mapping file (TSV): Representative_ID, Representative_Name,
                              Equivalent_IDs, Equivalent_Names
        - Allele list file (TSV): Representative_ID, Representative_Name
===============================================================================
"""

# ====================================================
# Imports
# ====================================================
import sys
import argparse
from pathlib import Path

# ====================================================
# Helper functions
# ====================================================

def load_allele_names(allele_list_file):
    """
    Load mapping of internal allele IDs to human-readable names.
    If no readable name is provided, use the ID itself.
    """
    name_map = {}
    with open(allele_list_file, "r") as f:
        for line in f:
            parts = line.strip().split(None, 1)
            if len(parts) == 2:
                name_map[parts[0]] = parts[1]
            elif len(parts) == 1:
                name_map[parts[0]] = parts[0]
    return name_map


def parse_cluster_file(cluster_file):
    """
    Parse CD-HIT .clstr file and extract representative sequences and members.
    Returns a dictionary: {representative: [members]}.
    """
    clusters = {}
    current_cluster = None
    rep_allele = None
    members = []

    with open(cluster_file, "r") as f:
        for line in f:
            line = line.strip()

            # Start of a new cluster
            if line.startswith(">Cluster"):
                if current_cluster is not None and rep_allele:
                    clusters[rep_allele] = members.copy()
                current_cluster = line.split()[1]
                members = []
                rep_allele = None
                continue

            # Member line (sequence entry)
            if line and line[0].isdigit():
                if ">" in line:
                    allele_part = line.split(">")[1]
                    allele = allele_part.split("...")[0]
                    members.append(allele)
                    if line.rstrip().endswith("*"):
                        rep_allele = allele

    # Save last cluster
    if current_cluster is not None and rep_allele:
        clusters[rep_allele] = members

    return clusters


def write_mapping_file(clusters, name_map, output_file):
    """
    Write mapping file linking representatives to all equivalent alleles.
    """
    with open(output_file, "w") as f:
        f.write("Representative_ID\tRepresentative_Name\tEquivalent_IDs\tEquivalent_Names\n")
        for rep_id in sorted(clusters.keys()):
            rep_name = name_map.get(rep_id, "NA")
            eq_ids = ",".join(clusters[rep_id])
            eq_names = ",".join([name_map.get(i, "NA") for i in clusters[rep_id]])
            f.write(f"{rep_id}\t{rep_name}\t{eq_ids}\t{eq_names}\n")


def write_alleles_file(clusters, name_map, output_file):
    """
    Write list of representative alleles and their human-readable names.
    """
    with open(output_file, "w") as f:
        for rep_id in sorted(clusters.keys()):
            rep_name = name_map.get(rep_id, "NA")
            f.write(f"{rep_id}\t{rep_name}\n")

# ====================================================
# Main
# ====================================================

def main():
    parser = argparse.ArgumentParser(
        description="Parse CD-HIT cluster output and generate mapping files."
    )
    parser.add_argument("--cluster-file", required=True, help="Path to CD-HIT .clstr file")
    parser.add_argument("--allele-list", required=True, help="Path to HLA allele list file (ID + readable name)")
    parser.add_argument("--out-mapping", required=True, help="Output TSV for representative→member mapping")
    parser.add_argument("--out-alleles", required=True, help="Output TSV for representative allele list")

    args = parser.parse_args()

    # -------------------------------
    # Validate input files
    # -------------------------------
    if not Path(args.cluster_file).exists():
        print(f"[ERROR] Cluster file not found: {args.cluster_file}", file=sys.stderr)
        sys.exit(1)
    if not Path(args.allele_list).exists():
        print(f"[ERROR] Allele list file not found: {args.allele_list}", file=sys.stderr)
        sys.exit(1)

    # -------------------------------
    # Load and parse data
    # -------------------------------
    print("[INFO] Loading allele name mappings...")
    name_map = load_allele_names(args.allele_list)
    print(f"[INFO] Loaded {len(name_map):,} allele names")

    print("[INFO] Parsing CD-HIT cluster file...")
    clusters = parse_cluster_file(args.cluster_file)
    print(f"[INFO] Found {len(clusters):,} representative alleles")

    # -------------------------------
    # Write outputs
    # -------------------------------
    print("[INFO] Writing mapping file...")
    write_mapping_file(clusters, name_map, args.out_mapping)
    print(f"[INFO] Mapping file written: {args.out_mapping}")

    print("[INFO] Writing representative allele list...")
    write_alleles_file(clusters, name_map, args.out_alleles)
    print(f"[INFO] Allele list written: {args.out_alleles}")

    # -------------------------------
    # Summary
    # -------------------------------
    print("------------------------------------------------------------")
    print(f"[DONE] Retained {len(clusters):,} representative pseudo-sequences")
    print("------------------------------------------------------------")


# ====================================================
# Entry point
# ====================================================
if __name__ == "__main__":
    main()
