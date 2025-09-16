#!/usr/bin/env python
import os, sys
import hail as hl
from pyspark.sql import SparkSession

# --- pin Java & Python from your conda env (safe no-ops if already set) ---
os.environ.setdefault("PYSPARK_PYTHON", sys.executable)
os.environ.pop("_JAVA_OPTIONS", None)
os.environ.pop("JAVA_TOOL_OPTIONS", None)
if "CONDA_PREFIX" in os.environ:
    os.environ["JAVA_HOME"] = os.environ["CONDA_PREFIX"]
    os.environ["PATH"] = os.pathsep.join([os.path.join(os.environ["JAVA_HOME"], "bin"),
                                          os.environ["PATH"]])
    ...

# --- detect Hadoop version with a tiny Spark session ---
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
bm_path = "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm"

# Step 1: just list the directory to confirm connectivity
try:
    listing = hl.hadoop_ls(bm_path)[:5]
    print("List OK (first few entries):", listing)
except Exception as e:
    print("[FAIL] Could not list BM path; S3A/jars/connectivity issue:")
    raise

# Step 2: try reading metadata (shape)
from hail.linalg import BlockMatrix
try:
    bm = BlockMatrix.read(bm_path)
    print("BM read OK; fetching shape…")
    print("BM shape:", bm.shape)
    print("[SUCCESS] S3A + BM metadata path is working.")
except Exception as e:
    print("[FAIL] BM read/shape failed; JVM likely exited due to S3/jar/memory issue.")
    raise
