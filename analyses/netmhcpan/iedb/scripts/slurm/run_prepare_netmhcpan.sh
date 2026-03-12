#!/bin/bash
#SBATCH -J ddebug_omain_enrichment
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/iedb/logs/prepare_netmhcpan_iedb_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

source activate /ix/djishnu/Priyamvada/envs/virauto

python /ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/iedb/scripts/python/prepare_netmhcpan_iedb.py