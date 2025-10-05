#!/bin/bash
#SBATCH -J netmhcpan_by_allele
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -t 4-00:00
#SBATCH --mem=8G
#SBATCH --array=1-100
#SBATCH --constraint=amd
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/type1/netmhcpan_%A_%a.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/type1/netmhcpan_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

############################################################
# 1. Input paths
############################################################
# Peptides: all viral/human 9-mers together
PEP_FILE=/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/human_9mers.fasta

# Split allele list files (e.g., 100 chunks)
ALLELE_FILE=$(ls /ix/djishnu/Priyamvada/virauto/data/HLA_alleles/type_I_chunks/* | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Output directory
OUTDIR=/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1_human_9_mers_v2
mkdir -p "$OUTDIR"

############################################################
# 2. Scratch handling
############################################################
# Use fast node-local /scratch if available; otherwise fall back to lab path
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

############################################################
# 3. Run NetMHCpan for each allele in this chunk
############################################################
while read -r short long; do
    # Skip blank or malformed lines
    [[ -z "$short" || -z "$long" ]] && continue

    # Sanitize allele name for filenames
    SAFE_NAME=$(echo "$long" | tr '*' '_' | tr ':' '_')

    # Define output file paths
    OUT_TXT="${OUTDIR}/${SAFE_NAME}_chunk${SLURM_ARRAY_TASK_ID}.txt"
    OUT_XLS="${OUTDIR}/${SAFE_NAME}_chunk${SLURM_ARRAY_TASK_ID}.xls"

    echo "[RUN] Processing allele $long → $OUT_XLS"

    # Run NetMHCpan with tabular (Excel) output
    netMHCpan -a "$short" -f "$PEP_FILE" -BA -xlsfile "$OUT_XLS" > "$OUT_TXT"
done < "$ALLELE_FILE"

############################################################
# 4. Cleanup and finish
############################################################
echo "[INFO] Job $SLURM_ARRAY_TASK_ID completed at $(date)"
# optional cleanup: uncomment if scratch gets big
# rm -rf "$TMPDIR"
