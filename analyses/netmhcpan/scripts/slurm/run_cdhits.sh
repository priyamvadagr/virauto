#!/usr/bin/env bash
#SBATCH -J get_NR_seq
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/type1_extract_NR_alleles
# ==============================================================================
# Script: extract_NR_HLA_alleles.sh
# Description:
#   Extracts pseudo-sequences for all alleles in HLA_allele_list from
#   NetMHCpan's MHC_pseudo.dat, collapses identical sequences (100% identity)
#   using CD-HIT, and outputs:
#       1. Representative FASTA (non-redundant)
#       2. Representative allele list (+ readable names)
#       3. Mapping file: representative → equivalent alleles (internal + readable)
#   The intermediate full FASTA is removed at the end.
# ==============================================================================

# ---- CONFIG ----
DATA_DIR="/ix/djishnu/Priyamvada/auto_immune/NetMHCpan/NetMHCpan/netMHCpan-4.2/Linux_x86_64/data"
ALLELE_LIST="${DATA_DIR}/HLA_allele_list"
PSEUDO_FILE="${DATA_DIR}/MHC_pseudo.dat"

OUT_DIR="/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_NR_alleles"
LOG_DIR="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs"
TMP_FASTA="${OUT_DIR}/tmp_all_pseudo.fasta"
OUT_FASTA="${OUT_DIR}/HLA_all_NR.fasta"
OUT_CLUSTER_INFO="${OUT_FASTA}.clstr"
OUT_ALLELES_WITH_NAMES="${OUT_DIR}/HLA_all_NR_alleles_with_names.txt"
OUT_MAPPING="${OUT_DIR}/HLA_pseudo_equivalence_groups.tsv"
OUT_LOG="${LOG_DIR}/cdhit_unique_log_$(date +%Y%m%d_%H%M%S).txt"

#Load cd-hit module 
module load cd-hit/4.8.1

# ---- 1. Build temporary FASTA ----
echo "[INFO] Extracting pseudo-sequences..."
> "${TMP_FASTA}"

cut -d ' ' -f1 "${ALLELE_LIST}" | while read -r allele; do
    [[ -z "$allele" ]] && continue
    seq=$(grep -w "${allele}" "${PSEUDO_FILE}" | awk '{print $2}')
    if [[ -n "$seq" ]]; then
        echo ">${allele}" >> "${TMP_FASTA}"
        echo "${seq}" >> "${TMP_FASTA}"
    fi
done

# ---- 2. Cluster identical sequences ----
echo "[INFO] Running CD-HIT (100% identity)..."
cd-hit -i "${TMP_FASTA}" -o "${OUT_FASTA}" -c 1.0 -n 5 -d 0 -M 16000 -T 4 > "${OUT_LOG}" 2>&1

# ---- 3. Parse clusters and create output files with Python ----
echo "[INFO] Parsing clusters and creating output files..."

# Path to the Python script (adjust if needed)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/parse_cdhit_clusters.py"

python3 "${PYTHON_SCRIPT}" \
    --cluster-file "${OUT_CLUSTER_INFO}" \
    --allele-list "${ALLELE_LIST}" \
    --out-mapping "${OUT_MAPPING}" \
    --out-alleles "${OUT_ALLELES_WITH_NAMES}"

# ---- 4. Clean up ----
rm -f "${TMP_FASTA}"

# ---- 5. Summary ----
echo "[INFO] Completed successfully!"
echo "-----------------------------------------------"
echo "Representative FASTA:     ${OUT_FASTA}"
echo "Allele list (with names): ${OUT_ALLELES_WITH_NAMES}"
echo "Mapping file:             ${OUT_MAPPING}"
echo "Cluster info:             ${OUT_CLUSTER_INFO}"
echo "CD-HIT log:               ${OUT_LOG}"
echo "-----------------------------------------------"
echo "[DONE] Retained $(wc -l < ${OUT_ALLELES_WITH_NAMES}) representative pseudo-sequences"