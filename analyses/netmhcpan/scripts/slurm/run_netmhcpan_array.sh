#!/bin/bash
#SBATCH -J netmhcpan_by_allele
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -t 4-00:00
#SBATCH --mem=8G
#SBATCH --array=1-100
#SBATCH --constraint=amd
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/type1_v2/netmhcpan_%A_%a.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/type1_v2/netmhcpan_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

############################################################
# 1. Input paths
############################################################
# Peptides file (human or viral)
PEP_FILE=/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/human_9mers.fasta

# Allele chunk file for this array task
ALLELE_FILE=$(ls /ix/djishnu/Priyamvada/virauto/data/HLA_alleles/type_I_chunks/* | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Optional: file listing alleles to skip (SAFE_NAME or long name acceptable)
EXCLUDE_FILE=/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/alleles_in_both.txt  # leave empty or comment out if not needed

# Output directory
OUTDIR=/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1_human_9_mers
mkdir -p "$OUTDIR"

############################################################
# 2. Scratch handling
############################################################
if [ -d "$SLURM_TMPDIR" ]; then
    TMPDIR="$SLURM_TMPDIR"
elif [ -d "/scratch" ] && [ -w "/scratch" ]; then
    TMPDIR="/scratch/$USER/netmhcpan_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$TMPDIR"
else
    TMPDIR="/ix/djishnu/Priyamvada/virauto/tmp/netmhcpan_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$TMPDIR"
fi
export TMPDIR
echo "[INFO] Using TMPDIR: $TMPDIR"

echo "[INFO] Job $SLURM_ARRAY_TASK_ID started at $(date)"
echo "[INFO] Using allele chunk: $ALLELE_FILE"
echo "[INFO] Peptide file: $PEP_FILE"
echo "[INFO] Output dir: $OUTDIR"
[ -f "$EXCLUDE_FILE" ] && echo "[INFO] Excluding alleles listed in $EXCLUDE_FILE"

############################################################
# 3. Prepare exclusion list
############################################################
declare -A EXCLUDE
if [ -f "$EXCLUDE_FILE" ]; then
    while read -r ex; do
        [[ -z "$ex" ]] && continue
        SAFE_EX=$(echo "$ex" | tr '*' '_' | tr ':' '_')
        EXCLUDE["$SAFE_EX"]=1
    done < "$EXCLUDE_FILE"
fi

############################################################
# 4. Run NetMHCpan for each allele in this chunk
############################################################
while read -r short long; do
    [[ -z "$short" || -z "$long" ]] && continue

    SAFE_NAME=$(echo "$long" | tr '*' '_' | tr ':' '_')

    # Skip excluded alleles
    if [[ ${EXCLUDE["$SAFE_NAME"]+_} ]]; then
        echo "[SKIP] Allele $long is in exclusion list"
        continue
    fi

    OUT_TXT="${OUTDIR}/${SAFE_NAME}_chunk${SLURM_ARRAY_TASK_ID}.txt"
    OUT_XLS="${OUTDIR}/${SAFE_NAME}_chunk${SLURM_ARRAY_TASK_ID}.xls"

    # Skip if output already exists
    if [[ -s "$OUT_XLS" ]]; then
        echo "[SKIP] Output already exists for $long"
        continue
    fi

    echo "[RUN] Processing allele $long → $OUT_XLS"
    netMHCpan -a "$short" -f "$PEP_FILE" -BA -xlsfile "$OUT_XLS" > "$OUT_TXT"
done < "$ALLELE_FILE"

############################################################
# 5. Cleanup and finish
############################################################
echo "[INFO] Job $SLURM_ARRAY_TASK_ID completed at $(date)"
# rm -rf "$TMPDIR"
