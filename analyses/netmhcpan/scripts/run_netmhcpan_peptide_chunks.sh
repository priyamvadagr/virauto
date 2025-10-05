#!/bin/bash
#SBATCH -J netmhcpan_text
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -t 06:00:00
#SBATCH --mem=8G
#SBATCH --array=1-1
#SBATCH --constraint=amd
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/logs/netmhcpan_txt_%A_%a.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/logs/netmhcpan_txt_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

# ======== CONFIG ========
module load gcc/15.1.0
export NETMHC=/ix/djishnu/Priyamvada/virauto/software/NetMHCpan/netMHCpan-4.2
export PATH=$NETMHC:$PATH

ALLELE_FILE=$(ls /ix/djishnu/Priyamvada/virauto/data/HLA_alleles/type_I_chunks/* | sed -n "${SLURM_ARRAY_TASK_ID}p")
PEP_FILE=$(ls /ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks/peptides/* | sed -n "${SLURM_ARRAY_TASK_ID}p")

OUTDIR=/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1_pep_chunks
mkdir -p $OUTDIR

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

echo "[INFO] Starting job $SLURM_ARRAY_TASK_ID"
echo "[INFO] Peptide file: $PEP_FILE"
echo "[INFO] Allele file: $ALLELE_FILE"

# Loop through alleles in this chunk
while read -r short long; do
  OUTFILE=${OUTDIR}/${long//\*/_}_chunk${SLURM_ARRAY_TASK_ID}.txt
  echo "[RUN] Allele $long -> $OUTFILE"
  netMHCpan -a $short -f $PEP_FILE -BA > $OUTFILE
done < $ALLELE_FILE

echo "[DONE] Job $SLURM_ARRAY_TASK_ID completed successfully."
