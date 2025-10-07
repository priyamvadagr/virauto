#!/bin/bash
#SBATCH -J netmhcpan_test
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -t 04:00:00
#SBATCH --mem=8G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_runtime/netmhcpan_test_%j.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_runtime/netmhcpan_test_%j.err


# --- Pick one peptide file and one allele list ---
ALLELE_FILE=/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_random_sampling/HLA_A_random200.txt
PEP_FILE=/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks/matched_pairs_chunk_001.fasta

OUTDIR=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_runtime
mkdir -p "$OUTDIR"

# --- Optional temporary directory handling ---
if [ -d "$SLURM_TMPDIR" ]; then
    TMPDIR="$SLURM_TMPDIR"
else
    TMPDIR="/ix/djishnu/Priyamvada/virauto/tmp/netmhcpan_test_${SLURM_JOB_ID}"
    mkdir -p "$TMPDIR"
fi
export TMPDIR

echo "[INFO] Running test job on $(hostname) at $(date)"
echo "[INFO] Allele file: $ALLELE_FILE"
echo "[INFO] Peptide file: $PEP_FILE"
echo "[INFO] TMPDIR: $TMPDIR"

#############################################################
#   LOG JOB START
############################################################
echo "[INFO] Starting NetMHCpan test"
echo "[INFO] Peptide file: $PEP_FILE"
echo "[INFO] Allele list: $ALLELE_FILE"
echo "[INFO] Output directory: $OUTDIR"
echo "[INFO] Host: $(hostname)"
echo "[INFO] Start time: $(date)"

############################################################
#  RUN NETMHCPAN FOR ALL ALLELES
############################################################
start_time=$(date +%s)

count=0
while read -r short long; do
    [[ -z "$short" || -z "$long" ]] && continue
    SAFE_NAME=$(echo "$long" | tr '*' '_' | tr ':' '_')

    OUT_TXT="${OUTDIR}/${SAFE_NAME}_chunk_test.txt"
    OUT_XLS="${OUTDIR}/${SAFE_NAME}_chunk_test.xls"

    echo "[RUN] ($((++count))) Running allele: $long"
    /usr/bin/time -v netMHCpan -a "$short" -f "$PEP_FILE" -BA -xls -xlsfile "$OUT_XLS" > "$OUT_TXT"

done < "$ALLELE_FILE"

end_time=$(date +%s)
runtime=$((end_time - start_time))

############################################################
#  SUMMARY
############################################################
echo "[INFO] Finished all alleles in $runtime seconds ($(echo "$runtime/60" | bc) min)"
echo "[INFO] Results written to: $OUTDIR"
echo "[INFO] End time: $(date)"