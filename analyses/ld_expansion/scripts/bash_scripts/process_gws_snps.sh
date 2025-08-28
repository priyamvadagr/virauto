#!/bin/bash
#SBATCH --job-name=process_gws
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH --output=process_gws_snps.out
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=htc

module load gcc/12.2.0
module load r/4.4.0

Rscript --no-save --no-restore --verbose /ix/djishnu/Priyamvada/auto_immune/1000G_LD/scripts/R_scripts/process_gwas_snps.R > /ix/djishnu/Priyamvada/auto_immune/1000G_LD/scripts/R_scripts/process_gwas_snps.Rout