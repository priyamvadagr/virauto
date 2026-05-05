#!/bin/bash
#SBATCH -J parse_netmhcpan
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/iedb/logs/parse_netmhcpan_iedb_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

source activate /ix/djishnu/Priyamvada/envs/virauto

python /ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/iedb/scripts/python/parse_netmhcpan_iedb_results.py