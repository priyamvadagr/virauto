"""
================================================================================
Molecular Mimicry Risk Classification
================================================================================
Goal:
    Aggregate NetMHCpan prediction files (viral–human peptide pairs) into long-format data,
    compute Δ binding affinity (ΔBA_score = viral - human),
    and write to Parquet partitioned by Allele.
Binding Behavior and Interpretation:
    - Equivalent_binding:
        Both peptides bind well and with similar affinity
        (|ΔBA_score| ≤ 0.5).
    - Viral_dominant:
        Viral peptide binds significantly stronger (ΔBA_score ≤ -0.5).
    - Human_dominant:
        Human peptide binds significantly stronger (ΔBA_score ≥ +0.5).
    - Non_binder:
        Viral peptide fails to meet binding threshold (Rank > 2).
        → Immunologically silent; not expected to contribute to mimicry.
Quantitative thresholds:
    - Binding threshold: 
        BA_Rank ≤ 2  (NetMHCpan standard for binders)
    - Affinity similarity: 
        |ΔBA_score| ≤ 0.5  (log-scale window)
CLI Arguments:
    --class   : HLA class (e.g., Class_I or Class_II)
    --type  : MHC type (e.g., A, B, C)
    --kmer   : Peptide length (e.g., 8, 9, 10)

Example Usage:
    python merge_netmhcpan_predictions.py --type A --class I --kmer 9
======================================================================
"""
#!/usr/bin/env python3
"""
Aggregate NetMHCpan prediction files (viral–human peptide pairs) into long-format data,
compute Δ binding affinity (ΔBA_score = viral - human),
and write to Parquet partitioned by Allele.
"""

# ===================================================
# Imports
# ===================================================
import pandas as pd
import os, re, glob, argparse
from tqdm import tqdm
import pyarrow as pa
import pyarrow.parquet as pq

# ===================================================
# Arguments
# ===================================================

parser = argparse.ArgumentParser(
    description="Aggregate NetMHCpan viral–human peptide binding predictions into long-format Parquet files partitioned by allele."
)

parser.add_argument(
    "--class",
    dest="class_",
    type=str,
    required=True,
    default="Class_I",
    help="MHC class letter (e.g., Class_I or Class_II). Default: Class_I"
)

parser.add_argument(
    "--hla_type",
    dest="hla_type_",
    type=str,
    required=False,
    help="MHC class letter (e.g., A, B, or C)"
)

parser.add_argument(
    "--kmer",
    dest="k_mer",
    type=int,
    required=True,
    help="Peptide length (number of residues), e.g., 8, 9, or 10."
)

parser.add_argument(
    "--netmhc_dir",
    type=str,
    default="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/",
    help="Path to the directory containing chunked NetMHCpan prediction files."
)

parser.add_argument(
    "--out_dir",
    type=str,
    default="/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/",
    help="Output directory for aggregated results."
)

parser.add_argument(
    "--batch_size",
    type=int,
    default=15,
    help="Number of files to process before writing each Parquet batch."
)

parser.add_argument(
    "--force_reprocess",
    action="store_true",
    help="Reprocess all files even if they appear in the processed manifest."
)

parser.add_argument(
    "--workers",
    type=int,
    default=4,
    help="Number of parallel workers to use."
)

args = parser.parse_args()
hla_class = args.class_
hla_type = f'HLA-{args.hla_type_}'
k_mer = f'{args.k_mer}_mers'
BATCH_SIZE = args.batch_size
MAX_WORKERS = args.workers
base_results = args.out_dir
netmhc_dir = args.netmhc_dir



# ===================================================
# Setup
# ===================================================
data_dir = os.path.join(netmhc_dir, f"{k_mer}/{hla_class}/{hla_type}_processed")
print(f"Input data directory: {data_dir}")
out_dir = os.path.join(base_results, f"{k_mer}/{hla_class}/{hla_type}")
os.makedirs(os.path.join(out_dir, "data"), exist_ok=True)
manifest_file = os.path.join(out_dir, "processed_manifest.txt")

processed = set()
if os.path.exists(manifest_file) and not args.force_reprocess:
    with open(manifest_file) as f:
        processed = set(line.strip() for line in f)

pattern = os.path.join(data_dir, f"{hla_type}_*_predictions_part*.txt.gz")
files = sorted(glob.glob(pattern))
files_to_process = [f for f in files if f not in processed]
print(f"Found {len(files_to_process)} new files to process")

# ===================================================
# Parameters
# ===================================================
BA_rank_thresh = 2.0
BA_score_delta = 0.5

# ===================================================
# Helper
# ===================================================
def process_file(fpath):
    df = pd.read_csv(fpath, sep="\t", compression="gzip")

    # Δ binding affinity
    df["BA_score_diff"] = df["BA_score_viral"] - df["BA_score_human"]

    # Classification
    vir_bind = (df["BA_Rank_viral"] <= BA_rank_thresh) #only viral peptide needs to bind
    df["Binding_Class"] = "Non_binder"
    df.loc[vir_bind & (df["BA_score_diff"].abs() <= BA_score_delta), "Binding_Class"] = "Equivalent_binding"
    df.loc[vir_bind & (df["BA_score_diff"] <= -BA_score_delta), "Binding_Class"] = "Viral_dominant"
    df.loc[vir_bind & (df["BA_score_diff"] >= BA_score_delta), "Binding_Class"] = "Human_dominant"

    # Select relevant columns
    df_long = df[[
        "pair_id", "uniprot_viral", "uniprot_human", "Peptide_viral", "Peptide_human",
        "source_file_human_human", "Allele",
        "BA_score_viral", "BA_score_human", "BA_score_diff",
        "BA_Rank_viral", "BA_Rank_human", "Binding_Class"
    ]].copy()

    df_long.rename(columns={
        "pair_id_uniprot_viral": "pair_id",
        "Peptide_viral": "peptide_viral",
        "Peptide_human": "peptide_human",
        "source_file_human_human": "prot2_origin",
        "Allele": "allele"
    }, inplace=True)

    # Add metadata
    df_long["allele_safe"] = df_long["allele"].str.replace('*', '_').replace(':', '_')  
    df_long["chunk_source"] = os.path.basename(fpath)
    df_long["k_mer"] = args.k_mer
    df_long["mhc_type"] = df_long["hla_class"] = df_long["allele"].apply(
    lambda x: (m.group(1)
                   if (m := re.search(r"HLA[-_: ]?([A-C])", str(x), flags=re.IGNORECASE))
                   else "Unknown")
    )
    return df_long

# ===================================================
# Main loop — partitioned write
# ===================================================
batch = []
for fpath in tqdm(files_to_process, desc="Processing chunks"):
    df_chunk = process_file(fpath)
    batch.append(df_chunk)
    processed.add(fpath)

    if len(batch) >= args.batch_size:
        df_batch = pd.concat(batch, ignore_index=True)
        table = pa.Table.from_pandas(df_batch)
        pq.write_to_dataset(
            table,
            root_path=os.path.join(out_dir, "data"),
            partition_cols=["allele_safe"],          # <<<<< partitioning by allele
            compression="snappy"
        )
        with open(manifest_file, "a") as mf:
            for f in processed:
                mf.write(f + "\n")
        batch = []

# Write remaining
if batch:
    df_batch = pd.concat(batch, ignore_index=True)
    df_batch.drop_duplicates(subset=["pair_id", "allele_safe"], inplace=True)
    table = pa.Table.from_pandas(df_batch)
    pq.write_to_dataset(
        table,
        root_path=os.path.join(out_dir, "data"),
        partition_cols=["allele_safe"],
        compression="snappy"
    )
    with open(manifest_file, "a") as mf:
        for f in processed:
            mf.write(f + "\n")

print(f"✅ Completed. Partitioned Parquet dataset stored at: {os.path.join(out_dir, 'data')}")
