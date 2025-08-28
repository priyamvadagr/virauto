#!/bin/bash
#SBATCH -N 1
#SBATCH --array=1-22
#SBATCH -J batch_ld
#SBATCH -t 1-00:00
#SBATCH --output=ld_out/MS/ld_chr%a-%A_0.001.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

GWAS_FILE=/ix/djishnu/Priyamvada/auto_immune/1000G_LD/gwas_snps/MS/MS_0.001.txt  # Pass like: gwas_snps/MS/MS_0.001.txt
OUTPUT_DIR=/ix/djishnu/Priyamvada/auto_immune/1000G_LD/LD_snps/MS
bash 05_run_ld.sh $GWAS_FILE $SLURM_ARRAY_TASK_ID $OUTPUT_DIR

