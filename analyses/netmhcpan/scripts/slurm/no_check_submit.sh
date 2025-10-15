#!/bin/bash

############################################################
# NO-CHECK SUBMIT - For early stage runs (< 20% complete)
# Submits all combinations, relies on NetMHCpan to skip existing
############################################################

set -uo pipefail

# --- Configuration ---
PEPTIDE_DIR="/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/split_chunks/chunk_3"
ALLELE_CHUNK_DIR="/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_NR_alleles/Type1_chunks"
OUTDIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/Type1_NR/All_types_chunks/chunk_3"
mkdir -p "$OUTDIR"
LOG_DIR="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/Type1_NR_chunks/"
JOB_OUT_DIR="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/Type1_NR_chunks/job_outputs/chunk_3"
mkdir -p "$JOB_OUT_DIR"

MAX_JOBS_PER_BATCH=100
MAX_BATCHES_TO_SUBMIT=10  # Higher for early stage

# Temporary directory for submission script 
TMP_DIR_SUB="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/Type1_NR_chunks/tmp"

# Temporary directory for job combinations
COMBO_DIR="$TMP_DIR_SUB/netmhcpan_combos_$$"
mkdir -p "$COMBO_DIR"
mkdir -p "$OUTDIR"
mkdir -p "$LOG_DIR"

echo "=========================================="
echo "  NetMHCpan No-Check Submit"
echo "=========================================="
echo ""
echo "This mode skips all checking and submits everything."
echo "Good for early-stage runs (< 20% complete)."
echo ""

start_time=$(date +%s)

############################################################
# STEP 1: Validate and count
############################################################
echo "[1/3] Scanning input files..."

mapfile -t PEPTIDE_FILES < <(find "$PEPTIDE_DIR" -name "*.fasta" -type f | sort)
mapfile -t ALLELE_CHUNKS < <(find "$ALLELE_CHUNK_DIR" -name "allele_chunk_*.txt" -type f | sort)

if [ ${#PEPTIDE_FILES[@]} -eq 0 ] || [ ${#ALLELE_CHUNKS[@]} -eq 0 ]; then
    echo "ERROR: No input files found"
    exit 1
fi

echo "  ✓ Peptide files: ${#PEPTIDE_FILES[@]}"
echo "  ✓ Allele chunks: ${#ALLELE_CHUNKS[@]}"
echo ""

TOTAL_COMBINATIONS=$((${#PEPTIDE_FILES[@]} * ${#ALLELE_CHUNKS[@]}))

# Quick completion estimate (optional, takes 2 seconds)
if [ -d "$OUTDIR" ]; then
    existing_count=$(find "$OUTDIR" -name "*.xls" -type f 2>/dev/null | wc -l)
    
    # Rough estimate of expected total
    first_allele_file="${ALLELE_CHUNKS[0]}"
    alleles_per_chunk=$(wc -l < "$first_allele_file")
    expected_total=$((${#PEPTIDE_FILES[@]} * ${#ALLELE_CHUNKS[@]} * alleles_per_chunk))
    
    completion=$(awk "BEGIN {printf \"%.1f\", ($existing_count/$expected_total)*100}")
    echo "Quick estimate: ~$completion% complete ($existing_count/$expected_total files)"
    echo ""
    
    if (( $(echo "$completion > 20" | bc -l) )); then
        echo "⚠  WARNING: You're >20% complete."
        echo "   Consider using optimized_check_files.sh instead to avoid duplicate work."
        echo ""
        read -p "Continue anyway? (yes/no): " proceed
        [[ ! $proceed =~ ^[Yy][Ee][Ss]$ ]] && exit 0
        echo ""
    fi
fi

############################################################
# STEP 2: Create all combinations (no checking!)
############################################################
echo "[2/3] Creating all combinations (no checking)..."

ALL_COMBOS=()
for pep_file in "${PEPTIDE_FILES[@]}"; do
    for allele_file in "${ALLELE_CHUNKS[@]}"; do
        ALL_COMBOS+=("${pep_file}|${allele_file}")
    done
done

echo "  ✓ Total combinations: ${#ALL_COMBOS[@]}"
echo ""

############################################################
# STEP 3: Create batches
############################################################
echo "[3/3] Creating batch files..."

BATCH_FILES=()
batch_num=1

for ((i=0; i<${#ALL_COMBOS[@]}; i+=MAX_JOBS_PER_BATCH)); do
    batch_file="$COMBO_DIR/batch_${batch_num}.txt"
    
    end=$((i + MAX_JOBS_PER_BATCH))
    [ $end -gt ${#ALL_COMBOS[@]} ] && end=${#ALL_COMBOS[@]}
    
    for ((j=i; j<end; j++)); do
        echo "${ALL_COMBOS[$j]}" >> "$batch_file"
    done
    
    BATCH_FILES+=("$batch_file")
    num_in_batch=$((end - i))
    echo "  Batch $batch_num: $num_in_batch combinations"
    ((batch_num++))
    
    [ $batch_num -gt $MAX_BATCHES_TO_SUBMIT ] && break
done

TOTAL_BATCHES=${#BATCH_FILES[@]}

setup_time=$(($(date +%s) - start_time))

############################################################
# CONFIRM SUBMISSION
############################################################
echo ""
echo "=========================================="
echo "  Ready to Submit"
echo "=========================================="
echo ""
echo "Summary:"
echo "  Total combinations: ${#ALL_COMBOS[@]}"
echo "  Combinations per batch: $MAX_JOBS_PER_BATCH"
echo "  Batches to submit: $TOTAL_BATCHES"
if [ ${#ALL_COMBOS[@]} -gt $((TOTAL_BATCHES * MAX_JOBS_PER_BATCH)) ]; then
    remaining=$((${#ALL_COMBOS[@]} - TOTAL_BATCHES * MAX_JOBS_PER_BATCH))
    echo "  Remaining (run again): $remaining combinations"
fi
echo ""
echo "Performance:"
echo "  Setup time: ${setup_time}s (no checking!)"
echo "  Estimated job time: ~3-20 min per task"
echo "  Total runtime: ~1-2 days"
echo ""

read -p "Submit $TOTAL_BATCHES batches? (yes/no): " confirm
[[ ! $confirm =~ ^[Yy][Ee][Ss]$ ]] && exit 0

cd "$JOB_OUT_DIR" 
############################################################
# SUBMIT JOBS
############################################################

echo pwd
echo ""
echo "Submitting jobs..."
echo ""

PREV_JOB_ID=""
SUBMITTED_JOBS=()

for ((i=0; i<${#BATCH_FILES[@]}; i++)); do
    batch_num=$((i + 1))
    batch_file="${BATCH_FILES[$i]}"
    num_combos=$(wc -l < "$batch_file")
    
    wrapper_script="$COMBO_DIR/wrapper_batch_${batch_num}.sh"
    
    cat > "$wrapper_script" << 'EOFWRAPPER'
#!/bin/bash
#SBATCH -J netmhcpan_batch
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-02:00:00
#SBATCH --mem=8G
#SBATCH --constraint=amd
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=smp

# Get combination
COMBO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${COMBO_FILE}")
PEPTIDE_FILE=$(echo "$COMBO" | cut -d'|' -f1)
ALLELE_FILE=$(echo "$COMBO" | cut -d'|' -f2)

[ ! -f "$PEPTIDE_FILE" ] || [ ! -f "$ALLELE_FILE" ] && exit 1

OUTDIR="PLACEHOLDER_OUTDIR"
mkdir -p "$OUTDIR"

TMPDIR="${SLURM_SCRATCH:-/tmp}/netmhcpan_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
export TMPDIR

PEPTIDE_NAME=$(basename "$PEPTIDE_FILE" .fasta)
ALLELE_NAME=$(basename "$ALLELE_FILE" .txt)

echo "[INFO] =========================================="
echo "[INFO] Task ${SLURM_ARRAY_TASK_ID} - Batch PLACEHOLDER_BATCH/PLACEHOLDER_TOTAL"
echo "[INFO] Peptide: $PEPTIDE_NAME"
echo "[INFO] Allele chunk: $ALLELE_NAME"
echo "[INFO] Start: $(date)"
echo "[INFO] =========================================="

success=0
failed=0
skipped=0

while read -r short long; do
    [[ -z "$short" || -z "$long" ]] && continue
    
    ((allele_count++))
    
    SAFE_NAME=$(echo "$long" | tr '*' '_' | tr ':' '_')
    OUT_XLS="${OUTDIR}/${SAFE_NAME}_${PEPTIDE_NAME}.xls"
    
    # Skip if output exists
    if [ -f "$OUT_XLS" ] && [ -s "$OUT_XLS" ]; then
        ((skip_count++))
        ((success_count++))
        continue
    fi
    
    echo "[RUN] ($allele_count) $long → $(basename "$OUT_XLS")"
    
    # Run NetMHCpan and suppress stdout/stderr
    netMHCpan -a "$short" -f "$PEPTIDE_FILE" -BA -xls -xlsfile "$OUT_XLS" > /dev/null 2>&1
    
    if [ $? -eq 0 ]; then
        ((success_count++))
        echo "      ✓ Success"
    else
        ((fail_count++))
        echo "      ✗ Failed"
    fi
done < "$ALLELE_FILE"

rm -rf "$TMPDIR"

echo ""
echo "[INFO] =========================================="
echo "[INFO] Summary"
echo "[INFO] Success: $success | Skipped: $skipped | Failed: $failed"
echo "[INFO] End: $(date)"
echo "[INFO] =========================================="

[ $failed -gt 0 ] && exit 1
exit 0
EOFWRAPPER

    sed -i "s|PLACEHOLDER_OUTDIR|$OUTDIR|g" "$wrapper_script"
    sed -i "s|PLACEHOLDER_BATCH|$batch_num|g" "$wrapper_script"
    sed -i "s|PLACEHOLDER_TOTAL|$TOTAL_BATCHES|g" "$wrapper_script"
    chmod +x "$wrapper_script"
    
    SBATCH_CMD="sbatch --array=1-${num_combos}"
    SBATCH_CMD="$SBATCH_CMD --export=ALL,COMBO_FILE=${batch_file}"
    SBATCH_CMD="$SBATCH_CMD --output=$JOB_OUT_DIR/batch_${batch_num}_%A.out"
    
    [ -n "$PREV_JOB_ID" ] && SBATCH_CMD="$SBATCH_CMD --dependency=afterok:${PREV_JOB_ID}"
    
    SBATCH_CMD="$SBATCH_CMD $wrapper_script"
    
    job_output=$(eval $SBATCH_CMD)
    job_id=$(echo "$job_output" | grep -oP '\d+')
    
    if [ -z "$job_id" ]; then
        echo "✗ Failed to submit batch $batch_num"
        exit 1
    fi
    
    echo "[$batch_num/$TOTAL_BATCHES] Job $job_id ($num_combos combos)"
    SUBMITTED_JOBS+=("$job_id")
    PREV_JOB_ID="$job_id"
    sleep 0.5
done

end_time=$(date +%s)
total_time=$((end_time - start_time))

echo ""
echo "=========================================="
echo "  Submission Complete!"
echo "=========================================="
echo ""
echo "Submitted $TOTAL_BATCHES batches in ${total_time}s"
echo "Job IDs: ${SUBMITTED_JOBS[@]}"
echo ""
echo "Monitor: squeue -u $USER"
echo "Cancel all: scancel ${SUBMITTED_JOBS[@]}"
echo ""
echo "After these complete, run this script again to submit more."
echo "=========================================="