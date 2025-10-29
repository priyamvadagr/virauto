#!/bin/bash
# ===========================================================================
# Script to submit Python mimicry analysis jobs via SLURM
# ===========================================================================
#SBATCH -J summary
#SBATCH --cluster htc
#SBATCH --constraint=genoa         # Use AMD EPYC 9374F nodes (fastest CPU)
#SBATCH -N 1
#SBATCH -c 16                      # 64 CPU cores on AMD node
#SBATCH -t 0-12:00
#SBATCH --mem=100G                 # Use high memory for large batch processing
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH -o /ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/logs/summary_%j.out

source activate /ix/djishnu/Priyamvada/envs/virauto

python /ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/python/peptide_pair_summary.py
