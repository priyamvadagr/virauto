#!/bin/bash
#SBATCH -J blastdb
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/alignment/logs/get_8mers.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

source activate /ix/djishnu/Priyamvada/envs/virauto/

python /ix/djishnu/Priyamvada/virauto/analyses/alignment/scripts/python/build_kmer_index.py \
    --db-fasta /ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta \
    --k 8 \
    --out-index /ix/djishnu/Priyamvada/virauto/data/refs/prot_kmers/uniprot_human_all_8mers.pkl