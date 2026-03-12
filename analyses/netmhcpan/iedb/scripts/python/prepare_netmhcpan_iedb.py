#!/usr/bin/env python3
"""
======================================================================
Script: prepare_netmhcpan_iedb.py
Description:
    Prepare IEDB mimicry pairs for NetMHCpan by grouping pairs by HLA 
    allele, writing per-allele peptide FASTA files, and generating a 
    SLURM array job submission script.

    Strategy:
      - Group pairs by HLA allele (each allele = one NetMHCpan run)
      - Write FASTA with both viral and human peptides per allele
      - Generate SLURM array job script

    This runs both the viral and human peptide through NetMHCpan so 
    you can compute ΔBA directly from the output.

Input:
    iedb_mhci_pairs_4digit_hla.csv.gz from filter_iedb_blast_hits.py

Output:
    - Per-allele FASTA files in chunks/
    - SLURM submission script
    - Allele manifest (maps array task ID → allele + FASTA)
    - Pair ID mapping with BLAST coordinates and UniProt ACs

Dependencies:
    pandas

Usage:
    python prepare_netmhcpan_iedb.py
======================================================================
"""

import pandas as pd
import os
import re
import hashlib
import string

# ====================================================
# Config
# ====================================================
INPUT_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla.csv.gz"

BASE_DIR = "/ix/djishnu/Priyamvada/virauto"
FASTA_DIR = os.path.join(BASE_DIR, "data/epitopes/iedb/netmhcpan/mhc_i/fasta_by_allele")
SCRIPT_DIR = os.path.join(BASE_DIR, "analyses/netmhcpan/iedb/scripts")
RESULT_DIR = os.path.join(BASE_DIR, "results/netmhcpan/iedb/mhc_i")
LOG_DIR = os.path.join(BASE_DIR, "logs")
PAIR_ID_MAP_FILE = os.path.join(BASE_DIR, "data/epitopes/iedb/netmhcpan/mhc_i/pair_id_mapping.csv.gz")
MANIFEST_FILE = os.path.join(BASE_DIR, "data/epitopes/iedb/netmhcpan/mhc_i/allele_manifest.tsv")
SUBMIT_SCRIPT = os.path.join(SCRIPT_DIR, "slurm/submit_netmhcpan_iedb.sh")

for d in [FASTA_DIR, SCRIPT_DIR, RESULT_DIR, LOG_DIR]:
    os.makedirs(d, exist_ok=True)


# ====================================================
# Step 1: Load and prepare data
# ====================================================
def load_data(filepath):
    print("=" * 60)
    print("Step 1: Loading pairs")
    print("=" * 60)

    df = pd.read_csv(filepath, low_memory=False)
    print(f"  Total pairs: {len(df):,}")

    # Validate required columns
    required = [
        "viral_sequence", "human_mimic_sequence", "mhc_allele",
        "structure_id", "sstart", "send", "hu_prot_id"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    print("  Human UniProt AC column: 'hu_prot_id'")

    # Drop rows with missing sequences / allele / human protein ID
    df = df.dropna(subset=["viral_sequence", "human_mimic_sequence", "mhc_allele", "hu_prot_id"])
    print(f"  After dropping NAs: {len(df):,}")

    # Keep allele string for downstream NetMHCpan formatting
    df["netmhcpan_allele"] = df["mhc_allele"].astype(str)

    print(f"  Unique HLA alleles: {df['netmhcpan_allele'].nunique()}")
    print(f"  Unique viral epitopes: {df['viral_sequence'].nunique():,}")
    print(f"  Unique human mimics: {df['human_mimic_sequence'].nunique():,}")
    print(f"  Unique human proteins: {df['hu_prot_id'].nunique():,}")

    # Ensure sstart/send are integers
    df["sstart"] = df["sstart"].astype(int)
    df["send"] = df["send"].astype(int)
    print(f"  BLAST coordinates (sstart/send): present for {df['sstart'].notna().sum():,} rows")

    return df


# ====================================================
# Pair ID generator
# ====================================================
def generate_pair_ids(df):
    """
    Generate unique 4-character alphanumeric codes for each
    viral-human pair. A pair is defined by:
      structure_id + hu_prot_id + viral_sequence + human_mimic_sequence

    Returns df with 'pair_id' column and a mapping table.
    """
    print(f"\n{'=' * 60}")
    print("Generating unique pair IDs")
    print(f"{'=' * 60}")

    df["pair_key"] = (
        df["structure_id"].astype(str) + "_" +
        df["hu_prot_id"].astype(str) + "_" +
        df["viral_sequence"].astype(str) + "_" +
        df["human_mimic_sequence"].astype(str)
    )

    unique_keys = df["pair_key"].unique()
    print(f"  Unique viral-human pairs: {len(unique_keys):,}")

    chars = string.ascii_uppercase + string.digits

    def key_to_code(key, index):
        code = ""
        n = index
        for _ in range(4):
            code = chars[n % 36] + code
            n //= 36
        return code

    key_to_id = {}
    for i, key in enumerate(sorted(unique_keys)):
        code = key_to_code(key, i)
        key_to_id[key] = code

    df["pair_id"] = df["pair_key"].map(key_to_id)

    map_cols = [
        "pair_id", "pair_key", "structure_id", "hu_prot_id",
        "viral_sequence", "human_mimic_sequence",
        "sstart", "send",
    ]

    pair_map = (
        df[map_cols]
        .drop_duplicates(subset=["pair_id", "sstart", "send"])
        .sort_values(["pair_id", "sstart"])
    )

    pair_map.to_csv(PAIR_ID_MAP_FILE, index=False, compression="gzip")
    print(f"  Generated {df['pair_id'].nunique():,} unique pair IDs")
    print(f"  Mapping rows (including multi-hit positions): {len(pair_map):,}")
    print(f"  Pair ID map: {PAIR_ID_MAP_FILE}")
    print(f"  Columns: {list(pair_map.columns)}")

    hits_per_pair = pair_map.groupby("pair_id").size()
    multi_hit = (hits_per_pair > 1).sum()
    if multi_hit > 0:
        print(f"  Pairs with multiple BLAST hit positions: {multi_hit:,}")
        print(f"    (max {hits_per_pair.max()} positions for one pair)")

    df = df.drop(columns=["pair_key"])

    return df, pair_map


# ====================================================
# Step 2: Group by allele, write FASTA files
# ====================================================
def write_fasta_by_allele(df):
    print(f"\n{'=' * 60}")
    print("Step 2: Writing per-allele FASTA files")
    print(f"{'=' * 60}")

    manifest = []
    allele_groups = df.groupby("netmhcpan_allele")

    for allele, group in allele_groups:
        # Safe filename
        safe_allele = allele.replace("*", "_").replace(":", "_").replace("-", "_")
        fasta_path = os.path.join(FASTA_DIR, f"{safe_allele}.fasta")

        # Collect unique peptides (both viral and human)
        # Headers kept short to avoid NetMHCpan truncation
        # Format: V_structureid or H_structureid_protid
        # The peptide sequence itself is NOT in the header — 
        # NetMHCpan reports it in the Peptide column
        seen = set()
        records = []

        for _, row in group.iterrows():
            sid = row["structure_id"]
            viral_seq = str(row["viral_sequence"]).strip().upper()
            human_seq = str(row["human_mimic_sequence"]).strip().upper()
            pair_id = row["pair_id"]

            # Skip sequences shorter than 8 aa (NetMHCpan minimum)
            if len(viral_seq) < 8 or len(human_seq) < 8:
                continue

            # Skip sequences with non-standard amino acids
            valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
            if not set(viral_seq).issubset(valid_aa) or not set(human_seq).issubset(valid_aa):
                continue

            # Viral peptide: V_XXXX (7 chars total, well within NetMHCpan limit)
            viral_id = f"V_{pair_id}"
            viral_key = f"{viral_id}_{viral_seq}"
            if viral_key not in seen:
                records.append((viral_id, viral_seq))
                seen.add(viral_key)

            # Human mimic: H_XXXX (7 chars total)
            human_id = f"H_{pair_id}"
            human_key = f"{human_id}_{human_seq}"
            if human_key not in seen:
                records.append((human_id, human_seq))
                seen.add(human_key)

        # Write FASTA
        with open(fasta_path, "w") as f:
            for header, seq in records:
                f.write(f">{header}\n{seq}\n")

        manifest.append({
            "task_id": len(manifest) + 1,
            "allele": allele,
            "safe_allele": safe_allele,
            "fasta_path": fasta_path,
            "n_pairs": len(group),
            "n_peptides": len(records),
        })

    # Save manifest
    manifest_df = pd.DataFrame(manifest)
    manifest_df.to_csv(MANIFEST_FILE, sep="\t", index=False)

    print(f"  Created {len(manifest)} FASTA files")
    print(f"  Total peptides across all files: {manifest_df['n_peptides'].sum():,}")
    print(f"  Manifest: {MANIFEST_FILE}")

    # Size distribution
    print(f"\n  Peptides per allele:")
    print(f"    Min: {manifest_df['n_peptides'].min()}")
    print(f"    Median: {manifest_df['n_peptides'].median():.0f}")
    print(f"    Max: {manifest_df['n_peptides'].max()}")
    print(f"    Mean: {manifest_df['n_peptides'].mean():.0f}")

    return manifest_df


# ====================================================
# Step 3: Generate SLURM submission script
# ====================================================
def write_submit_script(manifest_df):
    print(f"\n{'=' * 60}")
    print("Step 3: Generating SLURM submission script")
    print(f"{'=' * 60}")

    n_tasks = len(manifest_df)
    max_array = 100  # cluster limit

    # Split into batches of max_array tasks
    n_batches = (n_tasks + max_array - 1) // max_array
    
    batch_scripts = []
    for batch_idx in range(n_batches):
        start_task = batch_idx * max_array + 1
        end_task = min((batch_idx + 1) * max_array, n_tasks)
        batch_size = end_task - start_task + 1

        if n_batches > 1:
            batch_script = SUBMIT_SCRIPT.replace(".sh", f"_batch{batch_idx + 1}.sh")
        else:
            batch_script = SUBMIT_SCRIPT

        script_content = f"""#!/bin/bash
#SBATCH -J netmhcpan_iedb_b{batch_idx + 1}
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-02:00:00
#SBATCH --mem=8G
#SBATCH --array=1-{batch_size}
#SBATCH --output={LOG_DIR}/netmhcpan_b{batch_idx + 1}_%A_%a.out
#SBATCH --error={LOG_DIR}/netmhcpan_b{batch_idx + 1}_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=prg65@pitt.edu

############################################################
# NetMHCpan for IEDB mimicry pairs
# Batch {batch_idx + 1}/{n_batches} (tasks {start_task}-{end_task} of {n_tasks})
# Each array task = one HLA allele with all its peptides
# Manifest: {MANIFEST_FILE}
############################################################

MANIFEST="{MANIFEST_FILE}"
RESULT_DIR="{RESULT_DIR}"
mkdir -p "$RESULT_DIR"

# Map array task ID (1-{batch_size}) to manifest task ID ({start_task}-{end_task})
MANIFEST_TASK_ID=$(( SLURM_ARRAY_TASK_ID + {start_task - 1} ))

# Read task info from manifest
TASK_LINE=$(awk -F'\\t' -v id="$MANIFEST_TASK_ID" 'NR>1 && $1==id' "$MANIFEST")

if [ -z "$TASK_LINE" ]; then
    echo "[ERROR] No manifest entry for task $MANIFEST_TASK_ID"
    exit 1
fi

ALLELE=$(echo "$TASK_LINE" | cut -f2)
SAFE_ALLELE=$(echo "$TASK_LINE" | cut -f3)
FASTA_FILE=$(echo "$TASK_LINE" | cut -f4)
N_PAIRS=$(echo "$TASK_LINE" | cut -f5)
N_PEPTIDES=$(echo "$TASK_LINE" | cut -f6)

echo "=========================================="
echo "  NetMHCpan — IEDB Mimicry Pairs"
echo "=========================================="
echo "  Batch: {batch_idx + 1}/{n_batches}"
echo "  Array task: $SLURM_ARRAY_TASK_ID → Manifest task: $MANIFEST_TASK_ID"
echo "  Allele: $ALLELE"
echo "  FASTA: $FASTA_FILE"
echo "  Pairs: $N_PAIRS"
echo "  Peptides: $N_PEPTIDES"
echo "  Start: $(date)"
echo "=========================================="

# Validate input
if [ ! -f "$FASTA_FILE" ]; then
    echo "[ERROR] FASTA not found: $FASTA_FILE"
    exit 1
fi

# Setup temp directory
TMPDIR="${{SLURM_SCRATCH:-/tmp}}/netmhcpan_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"
mkdir -p "$TMPDIR"
export TMPDIR

# Output files
OUT_XLS="${{RESULT_DIR}}/${{SAFE_ALLELE}}.xls"
OUT_TXT="${{RESULT_DIR}}/${{SAFE_ALLELE}}.txt"

# Skip if already complete
if [ -f "$OUT_XLS" ] && [ -s "$OUT_XLS" ]; then
    echo "[SKIP] Output already exists: $OUT_XLS"
    rm -rf "$TMPDIR"
    exit 0
fi

# Convert allele format for netMHCpan
# NetMHCpan accepts HLA-A02:01 format (no asterisk)
NETMHCPAN_ALLELE=$(echo "$ALLELE" | sed 's/HLA-\\([ABC]\\)\\*/HLA-\\1/g')

echo "[RUN] netMHCpan -a $NETMHCPAN_ALLELE -f $FASTA_FILE -BA -xls -xlsfile $OUT_XLS"

# Run NetMHCpan
# -l 8,9,10,11,12,13,14 scores all MHC-I peptide lengths
# Without this, NetMHCpan defaults to 9-mers and slides a window
# across longer peptides, producing spurious sub-peptide predictions
netMHCpan \\
    -a "$NETMHCPAN_ALLELE" \\
    -f "$FASTA_FILE" \\
    -l 8,9,10,11,12,13,14 \\
    -BA \\
    -xls -xlsfile "$OUT_XLS" \\
    > "$OUT_TXT" 2>&1

EXIT_CODE=$?

# Cleanup
rm -rf "$TMPDIR"

if [ $EXIT_CODE -eq 0 ] && [ -f "$OUT_XLS" ] && [ -s "$OUT_XLS" ]; then
    echo ""
    echo "[SUCCESS] Output: $OUT_XLS"
    echo "[INFO] Lines in XLS: $(wc -l < "$OUT_XLS")"
else
    echo ""
    echo "[FAILED] Exit code: $EXIT_CODE"
    echo "[FAILED] Check: $OUT_TXT"
    # Show last 20 lines of output for debugging
    tail -20 "$OUT_TXT" 2>/dev/null
    exit 1
fi

echo ""
echo "=========================================="
echo "  End: $(date)"
echo "=========================================="
"""

        with open(batch_script, "w") as f:
            f.write(script_content)
        os.chmod(batch_script, 0o755)
        batch_scripts.append(batch_script)

        print(f"  Batch {batch_idx + 1}: tasks {start_task}-{end_task} "
              f"({batch_size} tasks) → {os.path.basename(batch_script)}")

    # Write a master submit script that chains all batches
    if n_batches > 1:
        master_script = SUBMIT_SCRIPT.replace(".sh", "_all.sh")
        with open(master_script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"# Submit all {n_batches} batches for NetMHCpan IEDB\n")
            f.write(f"# Total tasks: {n_tasks} across {n_batches} batches of max {max_array}\n\n")
            f.write("PREV_JOB_ID=\"\"\n\n")
            for i, script in enumerate(batch_scripts):
                f.write(f"echo \"Submitting batch {i + 1}/{n_batches}...\"\n")
                f.write("if [ -n \"$PREV_JOB_ID\" ]; then\n")
                f.write(f"    JOB_OUTPUT=$(sbatch --dependency=afterany:$PREV_JOB_ID {script})\n")
                f.write("else\n")
                f.write(f"    JOB_OUTPUT=$(sbatch {script})\n")
                f.write("fi\n")
                f.write("PREV_JOB_ID=$(echo $JOB_OUTPUT | grep -oP '\\d+')\n")
                f.write("echo \"  Job ID: $PREV_JOB_ID\"\n")
                f.write("sleep 1\n\n")
            f.write(f"echo \"\"\n")
            f.write(f"echo \"All {n_batches} batches submitted.\"\n")
            f.write(f"echo \"Monitor: squeue -u $USER\"\n")
        os.chmod(master_script, 0o755)
        print(f"\n  Master submit script: {master_script}")
        print(f"  Run: bash {master_script}")
    else:
        print(f"\n  Single batch — submit directly:")
        print(f"  sbatch {batch_scripts[0]}")

    print(f"\n  Total array tasks: {n_tasks}")
    print(f"  Batches: {n_batches} (max {max_array} per batch)")
    print(f"  Log directory: {LOG_DIR}")
    print(f"  Results directory: {RESULT_DIR}")

    return SUBMIT_SCRIPT


# ====================================================
# Summary
# ====================================================
def print_summary(df, manifest_df):
    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")

    print(f"  Input pairs: {len(df):,}")
    print(f"  Unique HLA alleles: {df['netmhcpan_allele'].nunique()}")
    print(f"  SLURM array tasks: {len(manifest_df)}")
    print(f"  Total peptides to score: {manifest_df['n_peptides'].sum():,}")

    print(f"\n  Top 10 alleles by pair count:")
    top = manifest_df.nlargest(10, "n_pairs")
    for _, row in top.iterrows():
        print(f"    {row['allele']}: {row['n_pairs']} pairs, {row['n_peptides']} peptides")

    print(f"\n  Estimated runtime:")
    print(f"    ~{len(manifest_df)} array tasks")
    print(f"    ~1-5 min per task (depends on peptide count)")
    print(f"    With full parallelization: ~5-10 min total")

    print(f"\n  To submit:")
    n_batches = (len(manifest_df) + 99) // 100
    if n_batches > 1:
        master = SUBMIT_SCRIPT.replace(".sh", "_all.sh")
        print(f"    bash {master}")
    else:
        print(f"    sbatch {SUBMIT_SCRIPT}")

    print(f"\n  After completion, parse results with:")
    print(f"    python parse_netmhcpan_iedb_results.py")


# ====================================================
# Main
# ====================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Prepare NetMHCpan for IEDB mimicry pairs")
    print("=" * 60)

    df = load_data(INPUT_FILE)
    df, pair_map = generate_pair_ids(df)
    manifest_df = write_fasta_by_allele(df)
    submit_script = write_submit_script(manifest_df)
    print_summary(df, manifest_df)

    print(f"\n✅ Done.")