#!/bin/bash
#SBATCH -J run_virauto_python
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/logs/run_domain_enrichment_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=prg65@pitt.edu

source activate /ix/djishnu/Priyamvada/envs/virauto

BASE_DIR=/ix/djishnu/Priyamvada/virauto/

python /ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/python/compute_domain_enrichment.py \
        --pair-map $BASE_DIR/data/epitopes/iedb/netmhcpan/mhc_i/swissprot/pair_id_mapping.csv.gz\
        --strong-mimicry $BASE_DIR/results/netmhcpan/iedb/mhc_i/swissprot/parsed/iedb_mhci_strong_mimicry.csv.gz \
        --protein2ipr $BASE_DIR/data/refs/interpro/protein2ipr_human.tsv.gz \
        --proteome $BASE_DIR/data/refs/uniprot/uniprot_human_all.fasta \
        --outdir $BASE_DIR/results/mimicry_analysis/iedb/mhc_i/swissprot/ORA_domain_enrichment_full_overlap \
        --fdr 0.05



