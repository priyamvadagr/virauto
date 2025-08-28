#!/bin/bash
#
#SBATCH -N 1                                     # Run on one node
#SBATCH -t 1-00:00                               # Runtime (1 day)
#SBATCH -J split_by_chr                          # Job name
#SBATCH --output=split_by_chr.out                # Output log file
#SBATCH --mail-type=END,FAIL                     # Notification settings
##SBATCH --cluster=smp                            # Cluster type (if applicable)
#SBATCH --mail-user=prg65@pitt.edu               # Email address

# Load PLINK 1.9 (needed for --chr)
module purge
module load plink/1.90b6.7

# ------------------------------------------------------------------------------
# Script: 04_split_by_chr.sh
# Purpose: Split PLINK reference files into one per chromosome for faster LD
# Input: data/EUR.hg37.id.ref.alt.{bed,bim,fam}
# Output: data/1000G.EUR.hg37.{CHR}.{bed,bim,fam}
# ------------------------------------------------------------------------------

set -euo pipefail
cd /ix/djishnu/Priyamvada/auto_immune/1000G_LD/data
mkdir chr_files

echo "📂 Splitting reference into per-chromosome files..."

for chr in {1..22}; do
    echo "🧬 Chromosome $chr"
    plink --bfile EUR.hg37.id.ref.alt \
          --chr $chr \
          --make-bed \
          --out chr_files/1000G.EUR.hg37.${chr}
done

echo "✅ Done. Output written to data/chr_files/1000G.EUR.hg37.{chr}.{bed,bim,fam}"

