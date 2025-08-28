#!/bin/bash
#
#SBATCH -N 1                              # Use one node
#SBATCH -t 0-02:00                        # Runtime in D-HH:MM
#SBATCH -J rename_snp_ids                # Job name
#SBATCH --output=rename_variants.out     # Output log file
#SBATCH --mail-type=END,FAIL             # Notifications
##SBATCH --cluster=smp                    # Cluster type (if relevant)
#SBATCH --mail-user=prg65@pitt.edu       # Email address for notifications

# Load plink2
module purge
module load plink/2.00_20230109

# ------------------------------------------------------------------------------
# Script: 03_rename_variants.sh
# Purpose: Rename SNP IDs to CHR:POS:REF:ALT format in the filtered 1000G dataset
# Input: data/raw.bed/.bim/.fam
# Output: data/EUR.hg37.id.ref.alt.{bed,bim,fam}
# ------------------------------------------------------------------------------

set -euo pipefail
cd /ix/djishnu/Priyamvada/auto_immune/1000G_LD/data

echo "🔄 Renaming SNP IDs to CHR:POS:REF:ALT format..."

# Create new BIM file with updated SNP IDs
awk '{ $2 = $1 ":" $4 ":" $5 ":" $6; print $0 }' raw.bim > EUR.hg37.new.id.bim

# Copy other files to match the new BIM
cp raw.bed EUR.hg37.new.id.bed
cp raw.fam EUR.hg37.new.id.fam

# Rebuild PLINK binary dataset to remove duplicates again (if any remain)
plink2 --bfile EUR.hg37.new.id \
       --rm-dup exclude-all \
       --make-bed \
       --out EUR.hg37.id.ref.alt

echo "✅ Finished renaming. Output written to: EUR.hg37.id.ref.alt.bed/.bim/.fam"

