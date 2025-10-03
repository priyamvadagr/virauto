#!/bin/bash
#SBATCH -J blastdb
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/blast/logs/makkeblastdb.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

module load blast-plus/2.14.1

makeblastdb -in /ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta -dbtype prot -out /ix/djishnu/Priyamvada/virauto/data/refs/blastdb/uniprot_human_all/uniprot_human_all_db
