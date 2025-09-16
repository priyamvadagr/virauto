#!/usr/bin/env python
import os, sys
import hail as hl
from hail.linalg import BlockMatrix
import numpy as np
import pandas as pd

# -------------------------------
# Parse command line arguments
# -------------------------------
# Usage: panukbb_ld_extract.py <variant_list.txt> <out_prefix> [POP]
#   <variant_list.txt> : SNP list (CHR:POS:REF:ALT, GRCh37), one per line
#   <out_prefix>       : prefix for output files
#   [POP]              : population (EUR=default, or AFR/AMR/CSA/EAS/MID)
if len(sys.argv) < 3:
    print("Usage: panukbb_ld_extract.py <variant_list.txt> <out_prefix> [POP]")
    print("POP: EUR (default), AFR, AMR, CSA, EAS, MID")
    sys.exit(1)

in_variants = sys.argv[1]
out_prefix  =  sys.argv[2]
pop         = sys.argv[3] if len(sys.argv) > 3 else "EUR"

# -------------------------------
# Initialize Hail with Spark
# (Assumes your Spark has S3A wired; if not, add spark.jars.packages as discussed)
# -------------------------------
spark_conf = {
    "spark.hadoop.fs.s3a.impl": "org.apache.hadoop.fs.s3a.S3AFileSystem",
    "spark.hadoop.fs.s3a.aws.credentials.provider": "org.apache.hadoop.fs.s3a.AnonymousAWSCredentialsProvider",
    "spark.hadoop.fs.s3a.path.style.access": "true",

    # More resilience to S3 blips
    "spark.hadoop.fs.s3a.connection.maximum": "200",
    "spark.hadoop.fs.s3a.attempts.maximum": "20",
    "spark.hadoop.fs.s3a.connection.establish.timeout": "60000",
    "spark.hadoop.fs.s3a.connection.timeout": "60000",
    "spark.hadoop.fs.s3a.socket.timeout": "120000",
    "spark.network.timeout": "600s",
    "spark.executor.heartbeatInterval": "60s",

    # Give the JVM room (use a compute node)
    "spark.driver.memory": "12g",

    # 👉 Pick the matching pair for *your* Hadoop:
    # Example for Hadoop 3.3.1:
    "spark.jars.packages": "org.apache.hadoop:hadoop-aws:3.3.4,com.amazonaws:aws-java-sdk-bundle:1.12.367",
}
hl.init(app_name="virauto_ld", master="local[*]", spark_conf=spark_conf)
# -------------------------------
# Paths to Pan-UKBB LD resources
# -------------------------------

bm_path  = f"s3a://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.bm"

# If you have a local copy of the variant index, point to it:
idx_path = f"file:///ix/djishnu/Priyamvada/virauto/data/refs/panukbb/UKBB.{pop}.ldadj.variant.ht"

# If you want to keep a fallback to S3, you could do:
# local_idx = f"/ix/djishnu/Priyamvada/virauto/data/refs/panukbb/UKBB.{pop}.ldadj.variant.ht"
# if os.path.exists(local_idx):
#     idx_path = f"file://{local_idx}"
# else:
#     idx_path = f"s3a://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.variant.ht"


print(f"[info] Reading BlockMatrix: {bm_path}")
bm = BlockMatrix.read(bm_path)
print("BM shape:", bm.shape)

print(f"[info] Reading variant index: {idx_path}")
ht_idx = hl.read_table(idx_path)  # contains 'idx' mapping variants -> BM indices

# -------------------------------
# Load input variants into Hail
# -------------------------------
qt = hl.import_table(in_variants, no_header=True).rename({'f0': 'variant'})
qt = qt.annotate(v = hl.parse_variant(qt.variant, reference_genome='GRCh37'))
qt = qt.key_by(locus=qt.v.locus, alleles=qt.v.alleles).drop('v')

# -------------------------------
# Join input variants to index
# -------------------------------
qt_join = qt.join(ht_idx, how='left')
n_total = qt_join.count()
n_missing = qt_join.filter(hl.is_missing(qt_join.idx)).count()
print(f"[info] Matched indices for {n_total - n_missing}/{n_total} variants in {pop}.")

# Collect matched variants + indices
rows = qt_join.select('variant','idx').collect()
kept_map = {r['variant']: r['idx'] for r in rows if r['idx'] is not None}
if not kept_map:
    raise RuntimeError("No input variants matched the Pan-UKBB index.")

# Preserve original file order for labels/rows/cols
with open(in_variants) as f:
    input_ids = [ln.strip() for ln in f if ln.strip()]
ids  = [v for v in input_ids if v in kept_map]
idxs = [int(kept_map[v]) for v in ids]

# -------------------------------
# Subset BlockMatrix to variants
# (Hail API: filter_rows/cols take one argument)
# -------------------------------
# make sure it's a plain Python list of ints
idxs = [int(x) for x in idxs]
sub = bm.filter_rows(idxs).filter_cols(idxs)

R = sub.to_numpy()  # upper triangle filled

# -------------------------------
# Post-process LD matrix
# -------------------------------
R = np.triu(R) + np.triu(R, 1).T
np.fill_diagonal(R, 1.0)
R2 = R ** 2

# -------------------------------
# Write outputs
# -------------------------------
os.makedirs(os.path.dirname(f"{out_prefix}.ld_r.tsv.gz") or ".", exist_ok=True)
pd.DataFrame(R,  index=ids, columns=ids).to_csv(f"{out_prefix}.ld_r.tsv.gz",  sep="\t", compression="gzip")
pd.DataFrame(R2, index=ids, columns=ids).to_csv(f"{out_prefix}.ld_r2.tsv.gz", sep="\t", compression="gzip")
print(f"[done] Wrote {out_prefix}.ld_r.tsv.gz and {out_prefix}.ld_r2.tsv.gz")

