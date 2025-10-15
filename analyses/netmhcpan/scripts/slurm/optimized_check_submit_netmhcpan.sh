#!/bin/bash

############################################################
# OPTIMIZED NETMHCPAN BATCH SUBMISSION
# Processes peptide × allele chunk combinations
# Much faster: processes 30 alleles per task instead of 3,000
# Optimized checking strategies for existing outputs
############################################################

# Exit on error but show helpful messages
set -uo pipefail

# --- Configuration ---
PEPTIDE_DIR="/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/split_chunks/chunk_2"
ALLELE_CHUNK_DIR="/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_NR_alleles/Type1_chunks"
OUTDIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/Type1_NR/All_types_chunks/chunk_2"
LOG_DIR="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/Type1_NR_chunks/"
JOB_OUT_DIR="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/Type1_NR_chunks/job_outputs"


# Job batching parameters
MAX_JOBS_PER_BATCH=100  # CRC cluster limit
MAX_BATCHES_TO_SUBMIT=13  # Safety limit for initial run

# Temporary directory for submission script 
TMP_DIR_SUB="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/Type1_NR_chunks/tmp"

# Temporary directory for job combinations
COMBO_DIR="$TMP_DIR_SUB/netmhcpan_combos_$$"
mkdir -p "$COMBO_DIR"
mkdir -p "$OUTDIR"
mkdir -p "$LOG_DIR"
mkdir -p "$JOB_OUT_DIR"

############################################################
# STEP 1: VALIDATE INPUTS
############################################################
echo "=========================================="
echo "  NetMHCpan Optimized Submission"
echo "=========================================="
echo ""
echo "[1/5] Validating inputs..."

# Check if allele chunks exist
if [ ! -d "$ALLELE_CHUNK_DIR" ]; then
    echo ""
    echo "ERROR: Allele chunk directory not found: $ALLELE_CHUNK_DIR"
    echo ""
    echo "Please run split_alleles.sh first to create allele chunks:"
    echo "  ./split_alleles.sh"
    echo ""
    exit 1
fi

# Check peptide directory exists
if [ ! -d "$PEPTIDE_DIR" ]; then
    echo ""
    echo "ERROR: Peptide directory not found: $PEPTIDE_DIR"
    echo ""
    echo "Please update PEPTIDE_DIR path in this script."
    echo ""
    exit 1
fi

# Count allele chunks
mapfile -t ALLELE_CHUNKS < <(find "$ALLELE_CHUNK_DIR" -name "allele_chunk_*.txt" -type f | sort)
if [ ${#ALLELE_CHUNKS[@]} -eq 0 ]; then
    echo "ERROR: No allele chunk files found in $ALLELE_CHUNK_DIR"
    exit 1
fi

# Count peptide files
mapfile -t PEPTIDE_FILES < <(find "$PEPTIDE_DIR" -name "*.fasta" -type f | sort)
if [ ${#PEPTIDE_FILES[@]} -eq 0 ]; then
    echo "ERROR: No FASTA files found in $PEPTIDE_DIR"
    exit 1
fi

echo "  ✓ Peptide files: ${#PEPTIDE_FILES[@]}"
echo "  ✓ Allele chunks: ${#ALLELE_CHUNKS[@]}"

# Calculate total combinations
TOTAL_COMBINATIONS=$((${#PEPTIDE_FILES[@]} * ${#ALLELE_CHUNKS[@]}))
echo "  ✓ Total combinations: $TOTAL_COMBINATIONS"
echo ""

############################################################
# STEP 2: ULTRA-OPTIMIZED CHECK FOR EXISTING OUTPUTS
# Multiple strategies to speed up when many files exist
############################################################
############################################################
# STEP 2: ULTRA-OPTIMIZED CHECK FOR EXISTING OUTPUTS
############################################################
echo "[2/5] Checking for existing outputs (optimized)..."
check_start=$(date +%s)

TMP_DIR="$TMP_DIR_SUB/netmhcpan_check_$$"
mkdir -p "$TMP_DIR"

# Build index once
echo "  Building index..."
find "$OUTDIR" -maxdepth 1 -type f -name "*.xls" -printf "%f\n" > "$TMP_DIR/existing_outputs.txt" 2>/dev/null

# Load into hash table
declare -A existing_outputs
while IFS= read -r filename; do
    existing_outputs["$filename"]=1
done < "$TMP_DIR/existing_outputs.txt"

num_existing=${#existing_outputs[@]}
echo "  ✓ Indexed $num_existing existing outputs"

# Early exit if 100% complete
total_expected=$((${#PEPTIDE_FILES[@]} * $(cat "$ALLELE_CHUNK_DIR"/allele_chunk_*.txt | wc -l)))

if [ "$num_existing" -ge "$total_expected" ]; then
    echo "  ✓ All outputs exist (100% complete)"
    rm -rf "$TMP_DIR"
    echo ""
    echo "✓ All combinations already processed!"
    exit 0
fi

echo "  Completion: $num_existing/$total_expected ($(awk "BEGIN {printf \"%.1f\", ($num_existing/$total_expected)*100}")%)"

# Optimized batch checking
echo "  Checking combinations in batches..."
UNPROCESSED_COMBOS=()
processed_count=0

for pep_file in "${PEPTIDE_FILES[@]}"; do
    pep_name=$(basename "$pep_file" .fasta)
    
    for allele_file in "${ALLELE_CHUNKS[@]}"; do
        # Extract all alleles and check in one pass
        has_missing=false
        
        while read -r short_allele long_allele; do
            [[ -z "$short_allele" || -z "$long_allele" ]] && continue
            
            safe_allele=$(echo "$long_allele" | tr '*:' '__')
            expected_filename="${safe_allele}_${pep_name}.xls"
            
            if [ -z "${existing_outputs[$expected_filename]+x}" ]; then
                has_missing=true
                break  # Early exit on first missing
            fi
        done < "$allele_file"
        
        if [ "$has_missing" = true ]; then
            UNPROCESSED_COMBOS+=("${pep_file}|${allele_file}")
        else
            ((processed_count++))
        fi
    done
done

rm -rf "$TMP_DIR"
check_end=$(date +%s)
check_time=$((check_end - check_start))

echo "  ✓ Processed: $processed_count"
echo "  ✗ Unprocessed: ${#UNPROCESSED_COMBOS[@]}"
echo "  ⏱  Check time: ${check_time}s (~$(awk "BEGIN {printf \"%.1f\", $check_time/60}") min)"
echo ""

############################################################
# STEP 3: CREATE BATCH FILES
############################################################
echo "[3/5] Creating batch files..."

BATCH_FILES=()
batch_num=1

for ((i=0; i<${#UNPROCESSED_COMBOS[@]}; i+=MAX_JOBS_PER_BATCH)); do
    batch_file="$COMBO_DIR/batch_${batch_num}.txt"
    
    # Get up to MAX_JOBS_PER_BATCH combinations
    end=$((i + MAX_JOBS_PER_BATCH))
    if [ $end -gt ${#UNPROCESSED_COMBOS[@]} ]; then
        end=${#UNPROCESSED_COMBOS[@]}
    fi
    
    # Write combinations to batch file
    for ((j=i; j<end; j++)); do
        echo "${UNPROCESSED_COMBOS[$j]}" >> "$batch_file"
    done
    
    BATCH_FILES+=("$batch_file")
    num_in_batch=$((end - i))
    echo "  Batch $batch_num: $num_in_batch combinations"
    ((batch_num++))
    
    # Safety check: don't create too many batches
    if [ $batch_num -gt $MAX_BATCHES_TO_SUBMIT ]; then
        echo ""
        echo "  (Limited to first $MAX_BATCHES_TO_SUBMIT batches for safety)"
        break
    fi
done

TOTAL_BATCHES=${#BATCH_FILES[@]}
echo "  Total batches to submit: $TOTAL_BATCHES"
echo ""

############################################################
# STEP 4: CONFIRM SUBMISSION
############################################################
echo "[4/5] Ready to submit..."
echo ""
echo "Summary:"
echo "  Peptide chunks: ${#PEPTIDE_FILES[@]}"
echo "  Allele chunks: ${#ALLELE_CHUNKS[@]} (~30 alleles each)"
echo "  Total combinations: $TOTAL_COMBINATIONS"
echo "  Already processed: $processed_count ($(awk "BEGIN {printf \"%.1f\", ($processed_count/$TOTAL_COMBINATIONS)*100}")%)"
echo "  Combinations to process: ${#UNPROCESSED_COMBOS[@]}"
echo "  Batches: $TOTAL_BATCHES (max $MAX_JOBS_PER_BATCH jobs per batch)"
echo "  Time per task: ~3-20 minutes"
echo "  Estimated completion: ~1-2 days (with good parallelization)"
echo ""
echo "Performance:"
echo "  Checking existing outputs took: ${check_time}s (~$(awk "BEGIN {printf \"%.1f\", $check_time/60}") min)"
echo ""

if [ $TOTAL_BATCHES -gt $MAX_BATCHES_TO_SUBMIT ]; then
    echo "NOTE: Only submitting first $MAX_BATCHES_TO_SUBMIT batches."
    echo "      Run this script again after completion to submit more."
    echo ""
fi

read -p "Proceed with submission? (yes/no): " confirm
if [[ ! $confirm =~ ^[Yy][Ee][Ss]$ ]]; then
    echo "Aborted"
    rm -rf "$COMBO_DIR"
    exit 0
fi

cd $TMP_DIR_SUB

############################################################
# STEP 5: SUBMIT JOBS WITH DEPENDENCIES
############################################################
echo ""
echo "[5/5] Submitting jobs..."
echo ""

PREV_JOB_ID=""
SUBMITTED_JOBS=()

for ((i=0; i<${#BATCH_FILES[@]}; i++)); do
    batch_num=$((i + 1))
    batch_file="${BATCH_FILES[$i]}"
    num_combos=$(wc -l < "$batch_file")
    
    echo "[$batch_num/$TOTAL_BATCHES] Submitting batch $batch_num ($num_combos combinations)..."
    
    # Create a complete self-contained wrapper script
    wrapper_script="$COMBO_DIR/wrapper_batch_${batch_num}.sh"
    
    cat > "$wrapper_script" << 'EOFWRAPPER'
#!/bin/bash
#SBATCH -J netmhcpan_auto
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-02:00:00
#SBATCH --mem=8G
#SBATCH --constraint=amd
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=htc
#SBATCH --out=slurm-%A.out

# Get the combination for this array task
COMBO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${COMBO_FILE}")

# Split into peptide and allele files
PEPTIDE_FILE=$(echo "$COMBO" | cut -d'|' -f1)
ALLELE_FILE=$(echo "$COMBO" | cut -d'|' -f2)

# Validate inputs
if [ -z "$PEPTIDE_FILE" ] || [ ! -f "$PEPTIDE_FILE" ]; then
    echo "[ERROR] Invalid peptide file: $PEPTIDE_FILE"
    exit 1
fi

if [ -z "$ALLELE_FILE" ] || [ ! -f "$ALLELE_FILE" ]; then
    echo "[ERROR] Invalid allele file: $ALLELE_FILE"
    exit 1
fi

# Output directory
OUTDIR="PLACEHOLDER_OUTDIR"
LOG_DIR="PLACEHOLDER_LOGDIR"
mkdir -p "$OUTDIR"
mkdir -p "$LOG_DIR"

# Setup temp directory
if [ -d "$SLURM_SCRATCH" ] && [ -w "$SLURM_SCRATCH" ]; then
    TMPDIR="$SLURM_SCRATCH/netmhcpan_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
else
    TMPDIR="/tmp/netmhcpan_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
fi
mkdir -p "$TMPDIR"
export TMPDIR

# Get identifiers for logging
PEPTIDE_NAME=$(basename "$PEPTIDE_FILE" .fasta)
ALLELE_NAME=$(basename "$ALLELE_FILE" .txt)

echo "[INFO] =========================================="
echo "[INFO] Batch PLACEHOLDER_BATCH_NUM/PLACEHOLDER_TOTAL_BATCHES - Task $SLURM_ARRAY_TASK_ID"
echo "[INFO] Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Host: $(hostname)"
echo "[INFO] Start: $(date)"
echo "[INFO] =========================================="
echo "[INFO] Peptide: $PEPTIDE_NAME"
echo "[INFO] Allele chunk: $ALLELE_NAME"
echo "[INFO] =========================================="

# Run NetMHCpan for each allele in this chunk
allele_count=0
success_count=0
fail_count=0
skip_count=0

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

# Cleanup
rm -rf "$TMPDIR"

echo ""
echo "[INFO] =========================================="
echo "[INFO] Summary"
echo "[INFO] Total: $allele_count | Success: $success_count | Skip: $skip_count | Fail: $fail_count"
echo "[INFO] End: $(date)"
echo "[INFO] =========================================="

if [ $fail_count -gt 0 ]; then
    exit 1
fi
EOFWRAPPER

    # Replace placeholders
    sed -i "s|PLACEHOLDER_OUTDIR|$OUTDIR|g" "$wrapper_script"
    sed -i "s|PLACEHOLDER_LOGDIR|$LOG_DIR|g" "$wrapper_script"
    sed -i "s|PLACEHOLDER_BATCH_NUM|$batch_num|g" "$wrapper_script"
    sed -i "s|PLACEHOLDER_TOTAL_BATCHES|$TOTAL_BATCHES|g" "$wrapper_script"
    
    chmod +x "$wrapper_script"
    
    # Build sbatch command
    SBATCH_CMD="sbatch"
    SBATCH_CMD="$SBATCH_CMD --array=1-${num_combos}"
    SBATCH_CMD="$SBATCH_CMD --export=ALL,COMBO_FILE=${batch_file},BATCH_NUMBER=${batch_num},TOTAL_BATCHES=${TOTAL_BATCHES}"
    SBATCH_CMD="$SBATCH_CMD --output=$JOB_OUT_DIR/batch_${batch_num}_%A.out"
    
    # Submit batch with dependency chain
    if [ -n "$PREV_JOB_ID" ]; then
        echo "  └─ Depends on job: $PREV_JOB_ID"
        job_output=$(sbatch --dependency=afterok:${PREV_JOB_ID} \
                            --array=1-${num_combos} \
                            --export=ALL,COMBO_FILE=${batch_file},BATCH_NUMBER=${batch_num},TOTAL_BATCHES=${TOTAL_BATCHES} \
                            "$wrapper_script")
    else
        job_output=$(sbatch --array=1-${num_combos} \
                            --export=ALL,COMBO_FILE=${batch_file},BATCH_NUMBER=${batch_num},TOTAL_BATCHES=${TOTAL_BATCHES} \
                            "$wrapper_script")
    fi

# Robust job ID extraction
job_id=$(echo "$job_output" | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) print $i}')
if [ -z "$job_id" ]; then
    echo "  ✗ Failed to extract job ID. sbatch said:"
    echo "    $job_output"
    exit 1
fi

echo "  ✓ Job ID: $job_id"
PREV_JOB_ID="$job_id"

    
    sleep 0.5
done

############################################################
# SUMMARY
############################################################
echo ""
echo "=========================================="
echo "  Submission Complete!"
echo "=========================================="
echo ""
echo "Submitted $TOTAL_BATCHES batches:"
for ((i=0; i<${#SUBMITTED_JOBS[@]}; i++)); do
    batch_num=$((i + 1))
    num_combos=$(wc -l < "${BATCH_FILES[$i]}")
    echo "  Batch $batch_num: Job ${SUBMITTED_JOBS[$i]} ($num_combos combinations)"
done
echo ""
echo "Jobs will run automatically with dependencies."
echo ""
echo "Monitor:"
echo "  ./monitor_netmhcpan_progress.sh"
echo "  watch -n 60 './monitor_netmhcpan_progress.sh'"
echo ""
echo "Check queue:"
echo "  squeue -u $USER"
echo ""
echo "Cancel all:"
echo "  scancel $(IFS=' '; echo "${SUBMITTED_JOBS[*]}")"
echo ""
echo "Batch info saved in: $COMBO_DIR"
echo "=========================================="

# Save job tracking
JOB_TRACKING="$COMBO_DIR/submitted_jobs.txt"
for ((i=0; i<${#SUBMITTED_JOBS[@]}; i++)); do
    batch_num=$((i + 1))
    echo "Batch_${batch_num}: ${SUBMITTED_JOBS[$i]}" >> "$JOB_TRACKING"
done

echo ""
echo "Run again after completion to submit remaining batches (if any)"
echo ""