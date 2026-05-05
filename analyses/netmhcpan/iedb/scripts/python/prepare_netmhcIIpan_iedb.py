#!/usr/bin/env python3
"""
======================================================================
Script: prepare_netmhciipan_iedb.py
Description:
    Prepare IEDB MHC Class II mimicry pairs for NetMHCIIpan.
    
    Handles Class II allele complexity:
      - DR alleles: single beta chain (DRA is monomorphic)
      - DQ alleles: alpha/beta pairs or beta-only with default alpha
      - DP alleles: alpha/beta pairs only (weak LD, no default alpha)
      - Serotype → representative 4-digit allele mapping
      - Drops: 'human', 'mouse', H2-*, locus-only, gene-only

    NetMHCIIpan allele format:
      DR:  DRB1_0101  (or DRB3_0202, DRB4_0101, DRB5_0101)
      DQ:  HLA-DQA10501-DQB10201
      DP:  HLA-DPA10103-DPB10401

Input:
    Filtered BLAST hits CSV from filter_iedb_blast_hits.py

Output:
    - Per-allele FASTA files
    - SLURM submission script
    - Allele manifest
    - Pair ID mapping

Usage:
    python prepare_netmhciipan_iedb.py
======================================================================
"""

import pandas as pd
import os
import re
import string

# ====================================================
# Config
# ====================================================
INPUT_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/mhc_ii/iedb_mhc_ii_filtered_blast_hits.csv.gz"

BASE_DIR = "/ix/djishnu/Priyamvada/virauto"
FASTA_DIR = os.path.join(BASE_DIR, "data/epitopes/iedb/netmhcpan/mhc_ii/fasta_by_allele")
SCRIPT_DIR = os.path.join(BASE_DIR, "analyses/netmhcpan/iedb/scripts")
RESULT_DIR = os.path.join(BASE_DIR, "results/netmhcpan/iedb/mhc_ii")
LOG_DIR = os.path.join(BASE_DIR, "analyses/netmhcpan/iedb/logs")
PAIR_ID_MAP_FILE = os.path.join(BASE_DIR, "data/epitopes/iedb/netmhcpan/mhc_ii/pair_id_mapping.csv.gz")
MANIFEST_FILE = os.path.join(BASE_DIR, "data/epitopes/iedb/netmhcpan/mhc_ii/allele_manifest.tsv")
SUBMIT_SCRIPT = os.path.join(SCRIPT_DIR, "slurm/submit_netmhciipan_iedb.sh")
ALLELE_REPORT = os.path.join(BASE_DIR, "data/epitopes/iedb/netmhcpan/mhc_ii/allele_conversion_report.tsv")

for d in [FASTA_DIR, SCRIPT_DIR, RESULT_DIR, LOG_DIR,
          os.path.dirname(PAIR_ID_MAP_FILE), os.path.dirname(SUBMIT_SCRIPT)]:
    os.makedirs(d, exist_ok=True)


# ====================================================
# Serotype → representative 4-digit allele mappings
# ====================================================
SEROTYPE_TO_ALLELE = {
    # DR serotypes
    "HLA-DR1":  "HLA-DRB1*01:01",
    "HLA-DR2":  "HLA-DRB1*15:01",
    "HLA-DR3":  "HLA-DRB1*03:01",
    "HLA-DR4":  "HLA-DRB1*04:01",
    "HLA-DR5":  "HLA-DRB1*11:01",
    "HLA-DR6":  "HLA-DRB1*13:01",
    "HLA-DR7":  "HLA-DRB1*07:01",
    "HLA-DR8":  "HLA-DRB1*08:01",
    "HLA-DR9":  "HLA-DRB1*09:01",
    "HLA-DR11": "HLA-DRB1*11:01",
    "HLA-DR12": "HLA-DRB1*12:01",
    "HLA-DR13": "HLA-DRB1*13:01",
    "HLA-DR14": "HLA-DRB1*14:01",
    "HLA-DR15": "HLA-DRB1*15:01",
    "HLA-DR16": "HLA-DRB1*16:01",
    "HLA-DR17": "HLA-DRB1*03:01",   # DR17 = DRB1*03:01 split
    "HLA-DR51": "HLA-DRB5*01:01",
    "HLA-DR52": "HLA-DRB3*01:01",
    "HLA-DR53": "HLA-DRB4*01:01",
    # DQ serotypes
    "HLA-DQ1":  "HLA-DQA1*01:01/DQB1*05:01",
    "HLA-DQ2":  "HLA-DQA1*05:01/DQB1*02:01",
    "HLA-DQ3":  "HLA-DQA1*03:01/DQB1*03:02",
    "HLA-DQ5":  "HLA-DQA1*01:01/DQB1*05:01",
    "HLA-DQ6":  "HLA-DQA1*01:02/DQB1*06:02",
    "HLA-DQ7":  "HLA-DQA1*03:01/DQB1*03:01",
    "HLA-DQ8":  "HLA-DQA1*03:01/DQB1*03:02",
    # DP serotypes
    "HLA-DPw2": "HLA-DPA1*01:03/DPB1*02:01",
    "HLA-DPw4": "HLA-DPA1*01:03/DPB1*04:01",
}

# DQB1 → most common DQA1 pairing (strong LD)
DQB1_TO_DQA1 = {
    "DQB1*02:01": "DQA1*05:01",
    "DQB1*02:02": "DQA1*02:01",
    "DQB1*03:01": "DQA1*03:01",
    "DQB1*03:02": "DQA1*03:01",
    "DQB1*03:03": "DQA1*03:02",
    "DQB1*04:02": "DQA1*04:01",
    "DQB1*05:01": "DQA1*01:01",
    "DQB1*05:02": "DQA1*01:01",
    "DQB1*05:03": "DQA1*01:01",
    "DQB1*05:04": "DQA1*01:01",
    "DQB1*06:01": "DQA1*01:03",
    "DQB1*06:02": "DQA1*01:02",
    "DQB1*06:03": "DQA1*01:03",
    "DQB1*06:04": "DQA1*01:02",
    "DQB1*06:09": "DQA1*01:02",
}


# ====================================================
# Allele classification and conversion
# ====================================================

def classify_and_convert_allele(allele):
    """
    Classify a Class II allele string and convert to NetMHCIIpan format.
    
    Returns:
        (netmhciipan_allele, resolution, notes)
        netmhciipan_allele = None if unusable
    """
    allele = allele.strip()

    # --- Drop non-HLA ---
    if allele.lower() in ("human", "mouse"):
        return None, "non_hla", f"host organism: {allele}"
    if allele.startswith("H2-"):
        return None, "non_hla", f"mouse MHC: {allele}"

    # --- Serotype mapping ---
    if allele in SEROTYPE_TO_ALLELE:
        mapped = SEROTYPE_TO_ALLELE[allele]
        result, _, _ = classify_and_convert_allele(mapped)
        return result, "serotype_mapped", f"{allele} → {mapped}"

    # --- Locus-only (HLA-DR, HLA-DQ, HLA-DP) ---
    if re.match(r"^HLA-D[RQP]$", allele):
        return None, "locus_only", f"no allele specified: {allele}"

    # --- Gene-only (HLA-DRB1, HLA-DQB1, etc.) ---
    if re.match(r"^HLA-D[RQP][AB]\d$", allele):
        return None, "gene_only", f"no allele specified: {allele}"
    if re.match(r"^HLA-D[RQP][AB]\d\*?$", allele):
        return None, "gene_only", f"no allele specified: {allele}"

    # --- DRA/DRB pair (e.g., HLA-DRA*01:01/DRB1*01:01) ---
    m = re.match(r"^HLA-DRA\*\d+:\d+/DRB(\d)\*(\d+:\d+)$", allele)
    if m:
        drb_gene = m.group(1)
        drb_allele = m.group(2).replace(":", "")
        return f"DRB{drb_gene}_{drb_allele}", "4digit", allele

    # --- Single DR beta chain (e.g., HLA-DRB1*01:01) ---
    m = re.match(r"^HLA-DRB(\d)\*(\d+:\d+)$", allele)
    if m:
        drb_gene = m.group(1)
        drb_allele = m.group(2).replace(":", "")
        return f"DRB{drb_gene}_{drb_allele}", "4digit", allele

    # --- DQ alpha/beta pair (e.g., HLA-DQA1*05:01/DQB1*02:01) ---
    m = re.match(r"^HLA-DQA1\*(\d+:\d+)/DQB1\*(\d+:\d+)$", allele)
    if m:
        dqa = m.group(1).replace(":", "")
        dqb = m.group(2).replace(":", "")
        return f"HLA-DQA1{dqa}-DQB1{dqb}", "4digit", allele

    # --- DQ beta-only (e.g., HLA-DQB1*03:02) — pair with default alpha ---
    m = re.match(r"^HLA-DQB1\*(\d+:\d+)$", allele)
    if m:
        dqb_key = f"DQB1*{m.group(1)}"
        dqa = DQB1_TO_DQA1.get(dqb_key)
        if dqa:
            dqa_digits = dqa.replace("DQA1*", "").replace(":", "")
            dqb_digits = m.group(1).replace(":", "")
            return f"HLA-DQA1{dqa_digits}-DQB1{dqb_digits}", "beta_only_mapped", f"{allele} + {dqa}"
        else:
            return None, "beta_only_unmapped", f"no default alpha for {allele}"

    # --- DP alpha/beta pair (e.g., HLA-DPA1*01:03/DPB1*04:01) ---
    m = re.match(r"^HLA-DPA1\*(\d+:\d+)/DPB1\*(\d+:\d+)$", allele)
    if m:
        dpa = m.group(1).replace(":", "")
        dpb = m.group(2).replace(":", "")
        return f"HLA-DPA1{dpa}-DPB1{dpb}", "4digit", allele

    # --- DP beta-only — skip (weak LD) ---
    m = re.match(r"^HLA-DPB1\*(\d+:\d+)$", allele)
    if m:
        return None, "dp_beta_only_skipped", f"weak LD, no default alpha: {allele}"

    # --- DP alpha-only — skip ---
    m = re.match(r"^HLA-DPA1\*(\d+:\d+)$", allele)
    if m:
        return None, "dp_alpha_only_skipped", f"no beta chain: {allele}"

    # --- DQ alpha-only — skip ---
    m = re.match(r"^HLA-DQA1\*(\d+:\d+)$", allele)
    if m:
        return None, "dq_alpha_only_skipped", f"no beta chain: {allele}"

    # --- Catch-all ---
    return None, "unrecognized", f"cannot parse: {allele}"


# ====================================================
# Step 1: Load and classify alleles
# ====================================================
def load_data(filepath):
    print("=" * 60)
    print("Step 1: Loading pairs and classifying alleles")
    print("=" * 60)

    df = pd.read_csv(filepath, low_memory=False)
    print(f"  Total pairs: {len(df):,}")

    required = ["viral_sequence", "human_mimic_sequence", "mhc_allele", "structure_id"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df = df.dropna(subset=["viral_sequence", "human_mimic_sequence", "mhc_allele"])
    print(f"  After dropping NAs: {len(df):,}")

    # Classify each allele
    conversion_results = df["mhc_allele"].apply(classify_and_convert_allele)
    df["netmhciipan_allele"] = [r[0] for r in conversion_results]
    df["allele_resolution"] = [r[1] for r in conversion_results]
    df["allele_notes"] = [r[2] for r in conversion_results]

    # Report
    print(f"\n  Allele classification:")
    for res, count in df["allele_resolution"].value_counts().items():
        print(f"    {res}: {count:,}")

    # Filter to usable alleles
    usable = df[df["netmhciipan_allele"].notna()].copy()
    dropped = df[df["netmhciipan_allele"].isna()]

    print(f"\n  Usable pairs: {usable.shape[0]:,}")
    print(f"  Dropped pairs: {dropped.shape[0]:,}")
    print(f"  Unique NetMHCIIpan alleles: {usable['netmhciipan_allele'].nunique()}")

    print(f"\n  Top 10 NetMHCIIpan alleles:")
    for allele, count in usable["netmhciipan_allele"].value_counts().head(10).items():
        print(f"    {allele}: {count:,}")

    return usable, df


# ====================================================
# Pair ID generator (same as Class I)
# ====================================================
def generate_pair_ids(df):
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

    # Build mapping table
    map_cols = [
        "pair_id", "pair_key", "structure_id", "hu_prot_id",
        "viral_sequence", "human_mimic_sequence",
        "sstart", "send",
    ]
    map_cols = [c for c in map_cols if c in df.columns]

    pair_map = (
        df[map_cols]
        .drop_duplicates(subset=["pair_id", "sstart", "send"])
        .sort_values(["pair_id", "sstart"])
    )

    pair_map.to_csv(PAIR_ID_MAP_FILE, index=False, compression="gzip")
    print(f"  Generated {df['pair_id'].nunique():,} unique pair IDs")
    print(f"  Mapping rows: {len(pair_map):,}")
    print(f"  Pair ID map: {PAIR_ID_MAP_FILE}")

    df = df.drop(columns=["pair_key"])

    return df, pair_map


# ====================================================
# Step 2: Write per-allele FASTA files
# ====================================================
def write_fasta_by_allele(df):
    print(f"\n{'=' * 60}")
    print("Step 2: Writing per-allele FASTA files")
    print(f"{'=' * 60}")

    manifest = []
    allele_groups = df.groupby("netmhciipan_allele")

    for allele, group in allele_groups:
        safe_allele = (
            allele.replace("*", "_").replace(":", "_")
            .replace("-", "_").replace("/", "_")
        )
        fasta_path = os.path.join(FASTA_DIR, f"{safe_allele}.fasta")

        seen = set()
        records = []

        for _, row in group.iterrows():
            viral_seq = str(row["viral_sequence"]).strip().upper()
            human_seq = str(row["human_mimic_sequence"]).strip().upper()
            pair_id = row["pair_id"]

            # NetMHCIIpan minimum length is 9
            if len(viral_seq) < 9 or len(human_seq) < 9:
                continue

            valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
            if not set(viral_seq).issubset(valid_aa) or not set(human_seq).issubset(valid_aa):
                continue

            viral_id = f"V_{pair_id}"
            viral_key = f"{viral_id}_{viral_seq}"
            if viral_key not in seen:
                records.append((viral_id, viral_seq))
                seen.add(viral_key)

            human_id = f"H_{pair_id}"
            human_key = f"{human_id}_{human_seq}"
            if human_key not in seen:
                records.append((human_id, human_seq))
                seen.add(human_key)

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

    manifest_df = pd.DataFrame(manifest)
    manifest_df.to_csv(MANIFEST_FILE, sep="\t", index=False)

    print(f"  Created {len(manifest)} FASTA files")
    print(f"  Total peptides: {manifest_df['n_peptides'].sum():,}")
    print(f"  Manifest: {MANIFEST_FILE}")

    print(f"\n  Peptides per allele:")
    print(f"    Min: {manifest_df['n_peptides'].min()}")
    print(f"    Median: {manifest_df['n_peptides'].median():.0f}")
    print(f"    Max: {manifest_df['n_peptides'].max()}")

    return manifest_df


# ====================================================
# Step 3: Generate SLURM submission script
# ====================================================
def write_submit_script(manifest_df):
    print(f"\n{'=' * 60}")
    print("Step 3: Generating SLURM submission script")
    print(f"{'=' * 60}")

    n_tasks = len(manifest_df)
    max_array = 100

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
#SBATCH -J netmhciipan_iedb_b{batch_idx + 1}
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-04:00:00
#SBATCH --mem=8G
#SBATCH --array=1-{batch_size}
#SBATCH --output={LOG_DIR}/netmhciipan_b{batch_idx + 1}_%A_%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=prg65@pitt.edu

############################################################
# NetMHCIIpan for IEDB MHC Class II mimicry pairs
# Batch {batch_idx + 1}/{n_batches} (tasks {start_task}-{end_task} of {n_tasks})
############################################################

set -u

MANIFEST="{MANIFEST_FILE}"
RESULT_DIR="{RESULT_DIR}"
LOG_DIR="{LOG_DIR}"
mkdir -p "$RESULT_DIR" "$LOG_DIR"

ERR_FILE="${{LOG_DIR}}/netmhciipan_b{batch_idx + 1}_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}.err"

MANIFEST_TASK_ID=$(( SLURM_ARRAY_TASK_ID + {start_task - 1} ))

TASK_LINE=$(awk -F'\\t' -v id="$MANIFEST_TASK_ID" 'NR>1 && $1==id' "$MANIFEST")

if [ -z "$TASK_LINE" ]; then
    echo "[ERROR] No manifest entry for task $MANIFEST_TASK_ID"
    {{
        echo "[ERROR] No manifest entry for task $MANIFEST_TASK_ID"
    }} > "$ERR_FILE"
    exit 1
fi

ALLELE=$(echo "$TASK_LINE" | cut -f2)
SAFE_ALLELE=$(echo "$TASK_LINE" | cut -f3)
FASTA_FILE=$(echo "$TASK_LINE" | cut -f4)
N_PAIRS=$(echo "$TASK_LINE" | cut -f5)
N_PEPTIDES=$(echo "$TASK_LINE" | cut -f6)

echo "=========================================="
echo "  NetMHCIIpan — IEDB MHC-II Mimicry Pairs"
echo "=========================================="
echo "  Batch: {batch_idx + 1}/{n_batches}"
echo "  Array task: $SLURM_ARRAY_TASK_ID → Manifest task: $MANIFEST_TASK_ID"
echo "  Allele: $ALLELE"
echo "  FASTA: $FASTA_FILE"
echo "  Pairs: $N_PAIRS"
echo "  Peptides: $N_PEPTIDES"
echo "  Start: $(date)"
echo "=========================================="

if [ ! -f "$FASTA_FILE" ]; then
    echo "[ERROR] FASTA not found: $FASTA_FILE"
    {{
        echo "[ERROR] FASTA not found: $FASTA_FILE"
        echo "[ERROR] Allele: $ALLELE"
        echo "[ERROR] Manifest task: $MANIFEST_TASK_ID"
    }} > "$ERR_FILE"
    exit 1
fi

TMPDIR="${{SLURM_SCRATCH:-/tmp}}/netmhciipan_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"
mkdir -p "$TMPDIR"
export TMPDIR

OUT_XLS="${{RESULT_DIR}}/${{SAFE_ALLELE}}.xls"
OUT_TXT="${{RESULT_DIR}}/${{SAFE_ALLELE}}.txt"

if [ -f "$OUT_XLS" ] && [ -s "$OUT_XLS" ]; then
    echo "[SKIP] Output already exists: $OUT_XLS"
    rm -rf "$TMPDIR"
    exit 0
fi

echo "[RUN] netMHCIIpan -a $ALLELE -f $FASTA_FILE -BA -xls -xlsfile $OUT_XLS"

# NetMHCIIpan for Class II:
#   - No -l flag (Class II has open binding groove, scores full peptide)
#   - -BA for binding affinity predictions
#   - Allele format: DRB1_0101, HLA-DQA10501-DQB10201, HLA-DPA10103-DPB10401
netMHCIIpan \\
    -a "$ALLELE" \\
    -f "$FASTA_FILE" \\
    -BA \\
    -xls -xlsfile "$OUT_XLS" \\
    > "$OUT_TXT" 2>&1

EXIT_CODE=$?

rm -rf "$TMPDIR"

if [ $EXIT_CODE -eq 0 ] && [ -f "$OUT_XLS" ] && [ -s "$OUT_XLS" ]; then
    echo ""
    echo "[SUCCESS] Output: $OUT_XLS"
    echo "[INFO] Lines in XLS: $(wc -l < "$OUT_XLS")"
else
    echo ""
    echo "[FAILED] Exit code: $EXIT_CODE"
    echo "[FAILED] Check: $OUT_TXT"
    tail -20 "$OUT_TXT" 2>/dev/null

    {{
        echo "[FAILED] Exit code: $EXIT_CODE"
        echo "[FAILED] Allele: $ALLELE"
        echo "[FAILED] FASTA: $FASTA_FILE"
        echo "[FAILED] Output text log: $OUT_TXT"
        echo ""
        echo "---- Last 50 lines of $OUT_TXT ----"
        tail -50 "$OUT_TXT" 2>/dev/null
    }} > "$ERR_FILE"

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

    # Master submit script
    if n_batches > 1:
        master_script = SUBMIT_SCRIPT.replace(".sh", "_all.sh")
        with open(master_script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"# Submit all {n_batches} batches for NetMHCIIpan IEDB\n\n")
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
            f.write(f"echo \"All {n_batches} batches submitted.\"\n")
        os.chmod(master_script, 0o755)
        print(f"\n  Master script: {master_script}")
    else:
        print(f"\n  Single batch — submit directly:")
        print(f"  sbatch {batch_scripts[0]}")

    return SUBMIT_SCRIPT

# ====================================================
# Main
# ====================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Prepare NetMHCIIpan for IEDB MHC-II mimicry pairs")
    print("=" * 60)

    usable, full_df = load_data(INPUT_FILE)

    # Save allele conversion report
    allele_report = (
        full_df[["mhc_allele", "netmhciipan_allele", "allele_resolution", "allele_notes"]]
        .drop_duplicates("mhc_allele")
        .sort_values("allele_resolution")
    )
    allele_report.to_csv(ALLELE_REPORT, sep="\t", index=False)
    print(f"\n  Allele conversion report: {ALLELE_REPORT}")

    usable, pair_map = generate_pair_ids(usable)
    manifest_df = write_fasta_by_allele(usable)
    submit_script = write_submit_script(manifest_df)

    # Summary
    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")
    print(f"  Input pairs: {len(full_df):,}")
    print(f"  Usable pairs (4-digit resolved): {len(usable):,}")
    print(f"  Dropped pairs: {len(full_df) - len(usable):,}")
    print(f"  Unique NetMHCIIpan alleles: {usable['netmhciipan_allele'].nunique()}")
    print(f"  SLURM array tasks: {len(manifest_df)}")
    print(f"  Total peptides to score: {manifest_df['n_peptides'].sum():,}")

    print(f"\n  Allele type breakdown:")
    for locus in ["DRB", "DQA", "DPA"]:
        n = usable["netmhciipan_allele"].str.startswith(locus if locus != "DPA" else "HLA-DPA").sum()
        if locus == "DRB":
            n = usable["netmhciipan_allele"].str.startswith("DRB").sum()
        elif locus == "DQA":
            n = usable["netmhciipan_allele"].str.startswith("HLA-DQA").sum()
        elif locus == "DPA":
            n = usable["netmhciipan_allele"].str.startswith("HLA-DPA").sum()
        print(f"    {locus}: {n:,} pairs")

    print(f"\n✅ Done.")