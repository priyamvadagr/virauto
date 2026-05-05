#!/bin/bash
#SBATCH -J run_virauto_python
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/logs/run_protein_enrichment_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=smp

source activate /ix/djishnu/Priyamvada/envs/virauto


python /ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/iedb/python/compute_protein_enrichment.py \
        --mimetopes /ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/iedb/mhc_i/swissprot/ORA_protein_enrichment/mimetopes_per_protein.tsv \
        --outdir /ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/iedb/mhc_i/swissprot/ORA_protein_enrichment/ \
        --fdr 0.05