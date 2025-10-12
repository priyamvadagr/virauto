#!/bin/bash

############################################################
# TEST VERSION - NetMHCpan Submission
# Only processes 2 peptides × 2 allele chunks = 4 jobs
# Perfect for testing before full run
############################################################

set -uo pipefail

# --- Configuration ---
PEPTIDE_DIR="/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers"
ALLELE_CHUNK_DIR="/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_NR_alleles/Type1_chunks"
OUTDIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/Type1_NR/test_run"
LOG_DIR="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/Type1_NR_chunks/test"

# TEST MODE: Only process this many
TEST_NUM_PEPTIDES=2
TEST_NUM_ALLELE_CHUNKS=2

mkdir -p "$OUTDIR"
mkdir -p "$LOG_DIR"

echo "=========================================="
echo "  NetMHCpan TEST Submission"
echo "=========================================="
echo ""
echo "This will submit a small test job to verify everything works."
echo ""

echo "=========================================="
echo "  NetMHCpan TEST Submission"
echo "=========================================="
echo ""
echo "This will submit a small test job to verify everything works."
echo ""

############################################################
# VALIDATE INPUTS
############################################################
echo "[1/4] Validating inputs..."

# Check allele chunks
if [ ! -d "$ALLELE_CHUNK_DIR" ]; then
    echo "ERROR: Allele chunk directory not found: $ALLELE_CHUNK_DIR"
    echo ""
    echo "Please run split_alleles.sh first!"
    exit 1
fi

mapfile -t ALLELE_CHUNKS < <(find "$ALLELE_CHUNK_DIR" -name "allele_chunk_*.txt" -type f | sort | head -n $TEST_NUM_ALLELE_CHUNKS)
if [ ${#ALLELE_CHUNKS[@]} -eq 0 ]; then
    echo "ERROR: No allele chunks found"
    exit 1
fi

# Check peptide files
if [ ! -d "$PEPTIDE_DIR" ]; then
    echo "ERROR: Peptide directory not found: $PEPTIDE_DIR"
    exit 1
fi

mapfile -t PEPTIDE_FILES < <(find "$PEPTIDE_DIR" -name "*.fasta" -type f | sort | head -n $TEST_NUM_PEPTIDES)
if [ ${#PEPTIDE_FILES[@]} -eq 0 ]; then
    echo "ERROR: No FASTA files found"
    exit 1
fi

echo "  ✓ Peptide files: ${#PEPTIDE_FILES[@]} (TEST MODE)"
echo "  ✓ Allele chunks: ${#ALLELE_CHUNKS[@]} (TEST MODE)"
echo ""

############################################################
# SHOW TEST PLAN
############################################################
echo "[2/4] Test plan..."
echo ""
echo "Will process:"

for pep_file in "${PEPTIDE_FILES[@]}"; do
    pep_name=$(basename "$pep_file" .fasta)
    echo "  Peptide: $pep_name"
    
    for allele_chunk in "${ALLELE_CHUNKS[@]}"; do
        chunk_name=$(basename "$allele_chunk" .txt)
        allele_count=$(wc -l < "$allele_chunk")
        echo "    × $chunk_name ($allele_count alleles)"
    done
done

TOTAL_JOBS=$((${#PEPTIDE_FILES[@]} * ${#ALLELE_CHUNKS[@]}))
echo ""
echo "Total test jobs: $TOTAL_JOBS"
echo "Expected outputs: ~$((TOTAL_JOBS * 30)) .xls files (if 30 alleles per chunk)"
echo "Estimated time: 5-15 minutes"
echo ""

############################################################
# CREATE TEST BATCH
############################################################
echo "[3/4] Creating test batch..."

TEST_BATCH_DIR="$LOG_DIR/tmp/netmhcpan_test_$$"
mkdir -p "$TEST_BATCH_DIR"

TEST_BATCH_FILE="$TEST_BATCH_DIR/test_batch.txt"

# Create combinations
for pep_file in "${PEPTIDE_FILES[@]}"; do
    for allele_chunk in "${ALLELE_CHUNKS[@]}"; do
        echo "${pep_file}|${allele_chunk}" >> "$TEST_BATCH_FILE"
    done
done

echo "  ✓ Created test batch with $TOTAL_JOBS combinations"
echo "  Location: $TEST_BATCH_FILE"
echo ""

############################################################
# CONFIRM SUBMISSION
############################################################
echo "[4/4] Ready to submit test job"
echo ""

read -p "Submit test job? (yes/no): " confirm
if [[ ! $confirm =~ ^[Yy][Ee][Ss]$ ]]; then
    echo "Aborted"
    rm -rf "$TEST_BATCH_DIR"
    exit 0
fi

############################################################
# CREATE WRAPPER SCRIPT
############################################################
echo ""
echo "Creating test wrapper..."

WRAPPER_SCRIPT="$TEST_BATCH_DIR/test_wrapper.sh"

cat > "$WRAPPER_SCRIPT" << 'EOF'
#!/bin/bash
#SBATCH -J netmhcpan_test
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-01:00:00
#SBATCH --mem=8G
#SBATCH --constraint=amd
#SBATCH --error=TEST_LOG_DIR/test_%A_%a.err
#SBATCH --output=TEST_LOG_DIR/test_%A_%a.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=smp

# Get the combination for this array task
COMBO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "TEST_BATCH_FILE")

# Split into peptide and allele files
PEPTIDE_FILE=$(echo "$COMBO" | cut -d'|' -f1)
ALLELE_FILE=$(echo "$COMBO" | cut -d'|' -f2)

echo "=========================================="
echo "TEST JOB - Task ${SLURM_ARRAY_TASK_ID}"
echo "=========================================="
echo "Peptide: $(basename $PEPTIDE_FILE)"
echo "Allele chunk: $(basename $ALLELE_FILE)"
echo "Start: $(date)"
echo "=========================================="

# Setup temp directory
if [ -d "$SLURM_SCRATCH" ] && [ -w "$SLURM_SCRATCH" ]; then
    TMPDIR="$SLURM_SCRATCH/test_netmhcpan_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
else
    TMPDIR="/tmp/test_netmhcpan_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
fi
mkdir -p "$TMPDIR"
export TMPDIR

# Output directory
OUTDIR="TEST_OUTDIR"
mkdir -p "$OUTDIR"

# Process each allele in this chunk
PEPTIDE_NAME=$(basename "$PEPTIDE_FILE" .fasta)
success=0
failed=0

while read -r short long; do
    [[ -z "$short" || -z "$long" ]] && continue
    
    SAFE_NAME=$(echo "$long" | tr '*' '_' | tr ':' '_')
    OUT_TXT="${OUTDIR}/${SAFE_NAME}_${PEPTIDE_NAME}.txt"
    OUT_XLS="${OUTDIR}/${SAFE_NAME}_${PEPTIDE_NAME}.xls"
    
    # Skip if exists
    if [ -f "$OUT_XLS" ] && [ -s "$OUT_XLS" ]; then
        ((success++))
        continue
    fi
    
    echo "Processing: $long"
    netMHCpan -a "$short" -f "$PEPTIDE_FILE" -BA -xls -xlsfile "$OUT_XLS" > "$OUT_TXT" 2>&1
    
    if [ $? -eq 0 ]; then
        ((success++))
    else
        ((failed++))
    fi
done < "$ALLELE_FILE"

# Cleanup
rm -rf "$TMPDIR"

echo ""
echo "=========================================="
echo "TEST SUMMARY"
echo "=========================================="
echo "Successful: $success"
echo "Failed: $failed"
echo "End: $(date)"
echo "=========================================="

if [ $failed -gt 0 ]; then
    exit 1
fi
EOF

# Replace placeholders
sed -i "s|TEST_BATCH_FILE|$TEST_BATCH_FILE|g" "$WRAPPER_SCRIPT"
sed -i "s|TEST_LOG_DIR|$LOG_DIR|g" "$WRAPPER_SCRIPT"
sed -i "s|TEST_OUTDIR|$OUTDIR|g" "$WRAPPER_SCRIPT"

chmod +x "$WRAPPER_SCRIPT"

############################################################
# SUBMIT TEST JOB
############################################################
echo "Submitting test job..."

job_output=$(sbatch --array=1-${TOTAL_JOBS} "$WRAPPER_SCRIPT")
job_id=$(echo "$job_output" | grep -oP '\d+')

if [ -z "$job_id" ]; then
    echo "✗ Failed to submit"
    echo "Output: $job_output"
    exit 1
fi

echo ""
echo "=========================================="
echo "  Test Job Submitted!"
echo "=========================================="
echo ""
echo "Job ID: $job_id"
echo "Array tasks: 1-${TOTAL_JOBS}"
echo ""
echo "Monitor:"
echo "  squeue -j $job_id"
echo "  watch -n 5 'squeue -j $job_id'"
echo ""
echo "Check logs:"
echo "  ls $LOG_DIR/"
echo "  tail -f $LOG_DIR/test_${job_id}_1.out"
echo ""
echo "Check outputs:"
echo "  ls $OUTDIR/"
echo "  ls $OUTDIR/*.xls | wc -l"
echo ""
echo "Cancel if needed:"
echo "  scancel $job_id"
echo ""
echo "Expected completion: 5-15 minutes"
echo "=========================================="
echo ""
echo "After test completes successfully, run the full version:"
echo "  ./auto_submit_netmhcpan.sh"
echo ""

# Save info
echo "Job ID: $job_id" > "$TEST_BATCH_DIR/job_info.txt"
echo "Test batch dir: $TEST_BATCH_DIR" >> "$TEST_BATCH_DIR/job_info.txt"