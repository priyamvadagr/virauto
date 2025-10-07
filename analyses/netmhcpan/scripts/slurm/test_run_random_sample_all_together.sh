#!/bin/bash
#SBATCH -J netmhcpan_test_all
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -t 04:00:00
#SBATCH --mem=8G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_runtime/netmhcpan_test_all_%j.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_runtime/netmhcpan_test_all_%j.err

# --- Pick one peptide file and one allele list ---
ALLELE_FILE=/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_random_sampling/HLA_A_random200.txt
PEP_FILE=/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks/matched_pairs_chunk_001.fasta

OUTDIR=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_runtime
mkdir -p "$OUTDIR"

# --- Optional temporary directory handling ---
if [ -d "$SLURM_TMPDIR" ]; then
    TMPDIR="$SLURM_TMPDIR"
else
    TMPDIR="/ix/djishnu/Priyamvada/virauto/tmp/netmhcpan_test_${SLURM_JOB_ID}"
    mkdir -p "$TMPDIR"
fi
export TMPDIR

echo "[INFO] Running test job on $(hostname) at $(date)"
echo "[INFO] Allele file: $ALLELE_FILE"
echo "[INFO] Peptide file: $PEP_FILE"
echo "[INFO] TMPDIR: $TMPDIR"

#############################################################
#   LOG JOB START
#############################################################
echo "[INFO] Starting NetMHCpan test - ALL ALLELES IN SINGLE JOB"
echo "[INFO] Peptide file: $PEP_FILE"
echo "[INFO] Allele list: $ALLELE_FILE"
echo "[INFO] Output directory: $OUTDIR"
echo "[INFO] Host: $(hostname)"
echo "[INFO] Start time: $(date)"

#############################################################
#  PREPARE ALLELE LIST
#############################################################
echo "[INFO] Reading allele file and preparing allele list..."

# Read the first column (short names) from allele file and create comma-separated list
ALLELE_LIST=$(awk '{print $1}' "$ALLELE_FILE" | grep -v '^$' | tr '\n' ',' | sed 's/,$//')

# Count alleles
NUM_ALLELES=$(echo "$ALLELE_LIST" | tr ',' '\n' | wc -l)
echo "[INFO] Number of alleles to process: $NUM_ALLELES"
echo "[INFO] Allele list (first 5): $(echo "$ALLELE_LIST" | cut -d',' -f1-5)..."

#############################################################
#  RUN NETMHCPAN FOR ALL ALLELES IN ONE JOB
#############################################################
echo "[INFO] Running netMHCpan with all $NUM_ALLELES alleles..."

ALLELE_BASENAME=$(basename "$ALLELE_FILE" .txt)
PEPTIDE_BASENAME=$(basename "$PEP_FILE" .fasta)

OUT_TXT="${OUTDIR}/${ALLELE_BASENAME}_${PEPTIDE_BASENAME}_ALL.txt"
OUT_XLS="${OUTDIR}/${ALLELE_BASENAME}_${PEPTIDE_BASENAME}_ALL.xls"

start_time=$(date +%s)

/usr/bin/time -v netMHCpan \
    -a "$ALLELE_LIST" \
    -f "$PEP_FILE" \
    -BA \
    -xls \
    -xlsfile "$OUT_XLS" \
    > "$OUT_TXT" 2>&1

exit_code=$?
end_time=$(date +%s)
runtime=$((end_time - start_time))

#############################################################
#  SUMMARY
#############################################################
echo ""
echo "=========================================="
echo "           JOB SUMMARY"
echo "=========================================="
echo "[INFO] Exit code: $exit_code"
echo "[INFO] Runtime: $runtime seconds ($(echo "scale=2; $runtime/60" | bc) minutes)"
echo "[INFO] Alleles processed: $NUM_ALLELES"
echo "[INFO] Output text file: $OUT_TXT"
echo "[INFO] Output Excel file: $OUT_XLS"

if [ $exit_code -eq 0 ]; then
    echo "[SUCCESS] netMHCpan completed successfully!"
    
    # Check output file size
    if [ -f "$OUT_TXT" ]; then
        filesize=$(du -h "$OUT_TXT" | cut -f1)
        lines=$(wc -l < "$OUT_TXT")
        echo "[INFO] Output file size: $filesize ($lines lines)"
    fi
else
    echo "[ERROR] netMHCpan failed with exit code $exit_code"
    echo "[ERROR] Check the output file for details: $OUT_TXT"
fi

echo "[INFO] End time: $(date)"
echo "=========================================="