#!/bin/bash
#SBATCH -J parse_netmhcpan
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/iedb/logs/parse_netmhciipan_iedb_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

source activate /ix/djishnu/Priyamvada/envs/virauto

python /ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/iedb/scripts/python/parse_netmhciipan_iedb.py \
        --mhc-class II \
        --result-dir /ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_ii/ \
        --pairs-file /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/mhc_ii/iedb_mhc_ii_filtered_blast_hits.csv.gz \
        --pair-map /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/netmhcpan/mhc_ii/pair_id_mapping.csv.gz \
        --out-dir /ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_ii/parsed/ \
        --rank-threshold 5.0

