#!/bin/bash
#SBATCH -J run_virauto_python
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/logs/run_plot_mimicry_summary_disease_risk_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=smp

source activate /ix/djishnu/Priyamvada/envs/virauto


python /ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/python/plot_mimicry_summary_disease_risk.py \
        --strong-mimicry  /ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_i/swissprot/parsed/iedb_mhci_strong_mimicry.csv.gz \
        --blast-input /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/mhc_i/iedb_mhci_pairs_4digit_hla_swissprot.csv.gz \
        --pair-map /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/netmhcpan/mhc_i/swissprot/pair_id_mapping.csv.gz \
        --hla-risk "/ix/djishnu/Priyamvada/pepSource/data/Ogishi_biorxiv_2019/HLA_autoimmunity_classification.txt" \
        --pssm-matrix /ix/djishnu/Priyamvada/pepSource/hla_domain_atlas/motif_atlas/results/mhc_i_motifs/9_mers/pssm_matrix_9mer.csv  \
        --corr-threshold 1 \
        --top-organisms 20 \
        --out-dir /ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/iedb/mhc_i/swissprot/disease_risk_figures_v2

