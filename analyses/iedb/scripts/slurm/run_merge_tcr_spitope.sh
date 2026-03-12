#!/bin/bash
#SBATCH -J query_iedb
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH --mem=100G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/iedb/logs/merge_iedb_epitopes_with_tcr.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

source activate /ix/djishnu/Priyamvada/envs/virauto

python /ix/djishnu/Priyamvada/virauto/analyses/iedb/scripts/python/merge_iedb_epitopes_with_tcr.py

