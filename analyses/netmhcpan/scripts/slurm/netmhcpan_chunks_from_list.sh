#!/bin/bash
#SBATCH -J netmhcpan_A_chunks
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00:00
#SBATCH --mem=8G
#SBATCH --array=1-100
#SBATCH --constraint=amd
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/type1_chunks/A_chunks/netmhcpan_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=smp

# --- Input files ---
ALLELE_FILE=/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_random_sampling/HLA_A_random200.txt
CHUNK_LIST=/ix/djishnu/Priyamvada/virauto/analyses/misc/chunk_lists/
# --- Output directory ---
OUTDIR=/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_A_chunks
mkdir -p "$OUTDIR"

# Get the specific peptide file for this array task
PEP_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHUNK_LIST")

# Check if we got a valid file
if [ -z "$PEP_FILE" ]; then
    echo "[ERROR] Could not read line ${SLURM_ARRAY_TASK_ID} from $CHUNK_LIST"
    exit 1
fi

if [ ! -f "$PEP_FILE" ]; then
    echo "[ERROR] Peptide file does not exist: $PEP_FILE"
    exit 1
fi



############################################################
# 2️⃣  TEMP DIRECTORY HANDLING
############################################################
if [ -d "$SLURM_SCRATCH" ]; then
    TMPDIR="$SLURM_SCRATCH"
    printf "[INFO] Using SLURM_SCRATCH: $TMPDIR\n"
elif [ -d "$SLURM_SCRATCH" ] && [ -w "$SLURM_SCRATCH" ]; then
    TMPDIR="$SLURM_SCRATCH/netmhcpan_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$TMPDIR"
else
    TMPDIR="/ix/djishnu/Priyamvada/virauto/tmp/netmhcpan_${SLURM_JOB_ID}"
    mkdir -p "$TMPDIR"
fi
export TMPDIR

echo "[INFO] =========================================="
echo "[INFO] Job $SLURM_ARRAY_TASK_ID started"
echo "[INFO] Host: $(hostname)"
echo "[INFO] Start time: $(date)"
echo "[INFO] =========================================="
echo "[INFO] Peptide file: $PEP_FILE"
echo "[INFO] Allele list: $ALLELE_FILE"
echo "[INFO] Output dir: $OUTDIR"
echo "[INFO] TMPDIR: $TMPDIR"
echo "[INFO] =========================================="

############################################################
# RUN NETMHCPAN FOR EACH ALLELE
############################################################
allele_count=0
success_count=0
fail_count=0

while read -r short long; do
    [[ -z "$short" || -z "$long" ]] && continue
    
    ((allele_count++))
    
    SAFE_NAME=$(echo "$long" | tr '*' '_' | tr ':' '_')
    CHUNK_NAME=$(basename "$PEP_FILE" .fasta)
    
    OUT_TXT="${OUTDIR}/${SAFE_NAME}_${CHUNK_NAME}.txt"
    OUT_XLS="${OUTDIR}/${SAFE_NAME}_${CHUNK_NAME}.xls"
    
    echo "[RUN] ($allele_count) Allele: $long → $(basename $OUT_XLS)"
    
    netMHCpan -a "$short" -f "$PEP_FILE" -BA -xls -xlsfile "$OUT_XLS" > "$OUT_TXT" 2>&1
    
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        ((success_count++))
        echo "      ✓ Success"
    else
        ((fail_count++))
        echo "      ✗ Failed (exit code: $exit_code)"
    fi

done < "$ALLELE_FILE"

############################################################
# SUMMARY
############################################################
echo ""
echo "[INFO] =========================================="
echo "[INFO] Job $SLURM_ARRAY_TASK_ID Summary"
echo "[INFO] =========================================="
echo "[INFO] Total alleles processed: $allele_count"
echo "[INFO] Successful: $success_count"
echo "[INFO] Failed: $fail_count"
echo "[INFO] End time: $(date)"
echo "[INFO] =========================================="

if [ $fail_count -gt 0 ]; then
    echo "[WARNING] Some alleles failed. Check output files for details."
    exit 1
fi

echo "[INFO] Job completed successfully"

############################################################
# CLEANUP AND LOGGING
############################################################
echo "[INFO] Cleaning up temporary directory..."
rm -rf "$TMPDIR"

echo "[INFO] Job completed successfully at $(date)"