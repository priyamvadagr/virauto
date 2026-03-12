#!/bin/bash
#SBATCH -J filter_iedb
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH --mem=100G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/iedb/logs/filter_iedb_blast_hits.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

source activate /ix/djishnu/Priyamvada/envs/virauto

python /ix/djishnu/Priyamvada/virauto/analyses/iedb/scripts/python/filter_iedb_blast_hits.py

