#!/bin/bash
#
#SBATCH -N 1                              # Use one node
#SBATCH -t 0-06:00                        # Runtime in D-HH:MM
#SBATCH -J run_ld                        # Job name
#SBATCH --output=ld_out/run_ld_chr%a.out        # Output log
#SBATCH --mail-type=END,FAIL             # Email notifications
#SBATCH --cluster=smp                    # Cluster type
#SBATCH --mail-user=prg65@pitt.edu       # Your email

# Load PLINK 1.9
module purge
module load plink/1.90b6.7

# ------------------------------------------------------------------------------
# Script: 05_run_ld.sh
# Purpose: Run pairwise LD calculation for GWAS SNPs on one chromosome
# Usage: bash scripts/05_run_ld.sh <snp_list> <chr_num> <out_dir>
# Example: bash scripts/05_run_ld.sh gwas_snps/MS/MS_0.001.txt 10 results/MS
# ------------------------------------------------------------------------------

set -euo pipefail

# Input arguments
snp_list=$1   # e.g., gwas_snps/MS/MS_0.001.txt
chr=$2        # e.g., 10
out_dir=$3    # e.g., results/MS

# Create output directory
mkdir -p ${out_dir}

# Strip file name for output prefix
base_name=$(basename "${snp_list}" .txt)
out_prefix="${out_dir}/${base_name}_chr${chr}"

# Reference prefix
ref_prefix="/ix/djishnu/Priyamvada/auto_immune/1000G_LD/data/chr_files/1000G.EUR.hg37.${chr}"

echo "📌 Running LD for ${base_name} on chromosome ${chr}"
echo "📁 Output directory: ${out_dir}"

plink --bfile ${ref_prefix} \
      --extract ${snp_list} \
      --r2 \
      --ld-window-kb 1000 \
      --ld-window 99999 \
      --ld-window-r2 0.8 \
      --out ${out_prefix}

echo "✅ LD results saved to: ${out_prefix}.ld" 
