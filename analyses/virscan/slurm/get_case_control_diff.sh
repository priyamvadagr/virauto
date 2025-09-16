#!/bin/bash
#SBATCH -J case_control
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/virscan/logs/get_case_control_differences.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

module load r/4.5.0

Rscript --no-save --no-restore --verbose /ix/djishnu/Priyamvada/virauto/analyses/virscan/R/get_case_control_differences.R  > /ix/djishnu/Priyamvada/virauto/analyses/virscan/R/get_case_control_differences.Rout
