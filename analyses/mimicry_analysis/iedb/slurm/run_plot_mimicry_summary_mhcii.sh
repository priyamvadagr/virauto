#!/bin/bash
#SBATCH -J run_virauto_python
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/logs/mhcii_mimicry_summary_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
##SBATCH --cluster=smp

source activate /ix/djishnu/Priyamvada/envs/virauto


python /ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/python/plot_mimicry_summary.py \
    --mhc-class II \
    --input /ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_ii/parsed/iedb_mhcii_strong_mimicry.csv.gz \
    --blast-input /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/mhc_ii/iedb_mhc_ii_pairs_4digit_hla.csv.gz \
    --pair-id-map /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/netmhcpan/mhc_ii/pair_id_mapping.csv.gz \
    --hla-risk /ix/djishnu/Priyamvada/pepSource/data/Ogishi_biorxiv_2019/HLA_autoimmunity_classification_standard.tsv \
    --iedb-fasta /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/fasta/iedb_epitopes_mhc_ii.fasta \
    --fig-dir /ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/iedb/mhc_ii/summary_figures