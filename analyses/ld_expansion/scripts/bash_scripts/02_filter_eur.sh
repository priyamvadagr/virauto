#!/bin/bash
#
#SBATCH -N 1                            # Ensure all cores are on one machine
#SBATCH -t 0-12:00                      # Runtime in D-HH:MM
#SBATCH -J filter_eur_samples           # Job name
#SBATCH --output=filter_eur.out         # Output log file
#SBATCH --mail-type=END,FAIL           # Email on job end or fail                
#SBATCH --mail-user=prg65@pitt.edu     # Email address for notifications
#SBATCH --mem=100000
module purge
module load plink/2.00_20230109

# ------------------------------------------------------------------------------
# Script: 02_filter_eur.sh
# Purpose: Filter EUR (non-FIN) individuals and convert PLINK2 → PLINK1 format
# Input: data/all_phase3.pgen, .pvar, .psam
# Output: data/raw.{bed,bim,fam}
# ------------------------------------------------------------------------------

set -euo pipefail

cd /ix/djishnu/Priyamvada/auto_immune/1000G_LD/data

echo "📋 Filtering EUR (non-FIN) individuals..."

# Create keep list: EUR samples excluding FINs
awk '($5=="EUR" && $6!="FIN") { print 0, $1 }' phase3_corrected.psam > eur.keep

# Exclude nothing (placeholder file)
echo "." > exclude.snps

# Convert to PLINK1 binary format and apply filters
plink2 --pgen all_phase3.pgen \
       --pvar all_phase3.pvar \
       --psam phase3_corrected.psam \
       --keep eur.keep \
       --exclude exclude.snps \
       --maf 0.01 \
       --snps-only just-acgt \
       --max-alleles 2 \
       --rm-dup exclude-all \
       --make-bed \
       --out raw

echo "✅ Output written to: data/raw.bed/.bim/.fam"

