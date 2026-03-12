#!/bin/bash
#SBATCH -J netmhcpan_iedb_b1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-02:00:00
#SBATCH --mem=8G
#SBATCH --array=1-85
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/logs/netmhcpan_b1_%A_%a.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/logs/netmhcpan_b1_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=prg65@pitt.edu

############################################################
# NetMHCpan for IEDB mimicry pairs
# Batch 1/1 (tasks 1-85 of 85)
# Each array task = one HLA allele with all its peptides
# Manifest: /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/netmhcpan/mhc_i/allele_manifest.tsv
############################################################

MANIFEST="/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/netmhcpan/mhc_i/allele_manifest.tsv"
RESULT_DIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_i"
mkdir -p "$RESULT_DIR"

# Map array task ID (1-85) to manifest task ID (1-85)
MANIFEST_TASK_ID=$(( SLURM_ARRAY_TASK_ID + 0 ))

# Read task info from manifest
TASK_LINE=$(awk -F'\t' -v id="$MANIFEST_TASK_ID" 'NR>1 && $1==id' "$MANIFEST")

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
echo "  Batch: 1/1"
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
TMPDIR="${SLURM_SCRATCH:-/tmp}/netmhcpan_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
export TMPDIR

# Output files
OUT_XLS="${RESULT_DIR}/${SAFE_ALLELE}.xls"
OUT_TXT="${RESULT_DIR}/${SAFE_ALLELE}.txt"

# Skip if already complete
if [ -f "$OUT_XLS" ] && [ -s "$OUT_XLS" ]; then
    echo "[SKIP] Output already exists: $OUT_XLS"
    rm -rf "$TMPDIR"
    exit 0
fi

# Convert allele format for netMHCpan
# NetMHCpan accepts HLA-A02:01 format (no asterisk)
NETMHCPAN_ALLELE=$(echo "$ALLELE" | sed 's/HLA-\([ABC]\)\*/HLA-\1/g')

echo "[RUN] netMHCpan -a $NETMHCPAN_ALLELE -f $FASTA_FILE -BA -xls -xlsfile $OUT_XLS"

# Run NetMHCpan
# -l 8,9,10,11,12,13,14 scores all MHC-I peptide lengths
# Without this, NetMHCpan defaults to 9-mers and slides a window
# across longer peptides, producing spurious sub-peptide predictions
netMHCpan \
    -a "$NETMHCPAN_ALLELE" \
    -f "$FASTA_FILE" \
    -l 8,9,10,11,12,13,14 \
    -BA \
    -xls -xlsfile "$OUT_XLS" \
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
