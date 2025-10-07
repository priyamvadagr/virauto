#!/bin/bash
#SBATCH -J netmhcpan_C_chunks
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00:00
#SBATCH --mem=8G
#SBATCH --array=1-100
#SBATCH --constraint=amd
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/type1_chunks/C_chunks/netmhcpan_%A_%a.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/type1_chunks/C_chunks/netmhcpan_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
##SBATCH --cluster=smp

# --- Input files ---
ALLELE_FILE=/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_random_sampling/HLA_C_random200.txt
PEP_FILE=$(ls /ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks/matched_pairs_chunk_*.fasta | sed -n "${SLURM_ARRAY_TASK_ID}p")

# --- Output directory ---
OUTDIR=/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_C_chunks
mkdir -p "$OUTDIR"

############################################################
# 2️⃣  TEMP DIRECTORY HANDLING
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

echo "[INFO] Job $SLURM_ARRAY_TASK_ID started on $(hostname) at $(date)"
echo "[INFO] Peptide file: $PEP_FILE"
echo "[INFO] Allele list: $ALLELE_FILE"
echo "[INFO] Output dir: $OUTDIR"
echo "[INFO] TMPDIR: $TMPDIR"

############################################################
# 3️⃣  RUN NETMHCPAN FOR EACH ALLELE
############################################################
while read -r short long; do
    [[ -z "$short" || -z "$long" ]] && continue

    SAFE_NAME=$(echo "$long" | tr '*' '_' | tr ':' '_')
    CHUNK_NAME=$(basename "$PEP_FILE" .fasta)

    OUT_TXT="${OUTDIR}/${SAFE_NAME}_${CHUNK_NAME}.txt"
    OUT_XLS="${OUTDIR}/${SAFE_NAME}_${CHUNK_NAME}.xls"

    echo "[RUN] Allele: $long → $OUT_XLS"
    netMHCpan -a "$short" -f "$PEP_FILE" -BA -xls -xlsfile "$OUT_XLS" > "$OUT_TXT"

done < "$ALLELE_FILE"

echo "[INFO] Job $SLURM_ARRAY_TASK_ID completed successfully at $(date)"
