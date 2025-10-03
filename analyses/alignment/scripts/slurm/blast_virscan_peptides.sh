#!/bin/bash
#SBATCH -J blastdb
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/alignment/logs/blast_virscan_peptides.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

## BLAST significant VirScan peptides against human proteome to check for potential similarities 
module load blast-plus/2.14.1
blastp -query /ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/VirScan_peptides_api.fasta \
        -db /ix/djishnu/Priyamvada/virauto/data/refs/blastdb/uniprot_human_all/uniprot_human_all_db \
       -out  /ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/all_peptides_blast_out.tsv \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen slen qstart qend sstart send evalue bitscore" \
       -evalue 10 \
       -word_size 2 \
       -num_threads 8


