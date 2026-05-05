#!/bin/bash
#SBATCH -J netmhciipan_iedb_b1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-04:00:00
#SBATCH --mem=8G
#SBATCH --array=1-72
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/iedb/logs/netmhciipan/netmhciipan_b1_%A_%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=prg65@pitt.edu

############################################################
# NetMHCIIpan for IEDB MHC Class II mimicry pairs
# Batch 1/1 (tasks 1-72 of 72)
############################################################

set -u

MANIFEST="/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/netmhcpan/mhc_ii/allele_manifest.tsv"
RESULT_DIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_ii"
LOG_DIR="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/iedb/logs/netmhciipan"
mkdir -p "$RESULT_DIR" "$LOG_DIR"

ERR_FILE="${LOG_DIR}/netmhciipan_b1_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

MANIFEST_TASK_ID=$(( SLURM_ARRAY_TASK_ID + 0 ))

TASK_LINE=$(awk -F'\t' -v id="$MANIFEST_TASK_ID" 'NR>1 && $1==id' "$MANIFEST")

if [ -z "$TASK_LINE" ]; then
    echo "[ERROR] No manifest entry for task $MANIFEST_TASK_ID"
    {
        echo "[ERROR] No manifest entry for task $MANIFEST_TASK_ID"
    } > "$ERR_FILE"
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
echo "  Batch: 1/1"
echo "  Array task: $SLURM_ARRAY_TASK_ID → Manifest task: $MANIFEST_TASK_ID"
echo "  Allele: $ALLELE"
echo "  FASTA: $FASTA_FILE"
echo "  Pairs: $N_PAIRS"
echo "  Peptides: $N_PEPTIDES"
echo "  Start: $(date)"
echo "=========================================="

if [ ! -f "$FASTA_FILE" ]; then
    echo "[ERROR] FASTA not found: $FASTA_FILE"
    {
        echo "[ERROR] FASTA not found: $FASTA_FILE"
        echo "[ERROR] Allele: $ALLELE"
        echo "[ERROR] Manifest task: $MANIFEST_TASK_ID"
    } > "$ERR_FILE"
    exit 1
fi

TMPDIR="${SLURM_SCRATCH:-/tmp}/netmhciipan_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
export TMPDIR

OUT_XLS="${RESULT_DIR}/${SAFE_ALLELE}.xls"
OUT_TXT="${RESULT_DIR}/${SAFE_ALLELE}.txt"

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
netMHCIIpan \
    -a "$ALLELE" \
    -f "$FASTA_FILE" \
    -BA \
    -xls -xlsfile "$OUT_XLS" \
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

    {
        echo "[FAILED] Exit code: $EXIT_CODE"
        echo "[FAILED] Allele: $ALLELE"
        echo "[FAILED] FASTA: $FASTA_FILE"
        echo "[FAILED] Output text log: $OUT_TXT"
        echo ""
        echo "---- Last 50 lines of $OUT_TXT ----"
        tail -50 "$OUT_TXT" 2>/dev/null
    } > "$ERR_FILE"

    exit 1
fi

echo ""
echo "=========================================="
echo "  End: $(date)"
echo "=========================================="
