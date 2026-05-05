#!/bin/bash
#SBATCH -J run_virauto_python
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/logs/run_protein_enrichment_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=smp

source activate /ix/djishnu/Priyamvada/envs/virauto


python /ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/python/compute_mimetopes_per_protein.py \
        --strong-mimicry /ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_i/swissprot/parsed/iedb_mhci_strong_mimicry.csv.gz \
        --pair-map /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/netmhcpan/mhc_i/swissprot/pair_id_mapping.csv.gz \
        --kmer-counts /ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_kmer_counts_filtered.tsv \
        --outfile /ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/iedb/mhc_i/swissprot/ORA_protein_enrichment/mimetopes_per_protein.tsv