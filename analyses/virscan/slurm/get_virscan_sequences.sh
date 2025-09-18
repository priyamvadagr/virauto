#!/bin/bash
#SBATCH -J fasta_map
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/virscan/logs/fasta_map_%x_%j.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/virscan/logs/fasta_map_%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

echo "[info] Job started on $(hostname) at $(date)"

# ----- ensure logs and results directories exist -----
mkdir -p /ix/djishnu/Priyamvada/virauto/analyses/virscan/logs
mkdir -p /ix/djishnu/Priyamvada/virauto/results/virscan

# ----- activate your environment -----
# If your cluster uses conda init:
source activate /ix/djishnu/Priyamvada/envs/virauto/

# ----- anchor to repo root -----
cd /ix/djishnu/Priyamvada/virauto || exit 1
echo "[info] Working directory: $(pwd)"

# ----- run script -----
python analyses/virscan/python/get_virscan_sequences.py

echo "[done] Job finished at $(date)"
