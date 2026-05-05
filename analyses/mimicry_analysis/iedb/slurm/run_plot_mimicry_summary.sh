#!/bin/bash
#SBATCH -J run_virauto_python
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/logs/run_plot_mimicry_summary_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=smp

source activate /ix/djishnu/Priyamvada/envs/virauto


python /ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/python/plot_mimicry_summary.py 