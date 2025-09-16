#!/bin/bash
#SBATCH -J ld_extract
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/ld_expansion/logs/ld_extract_%x_%j.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/ld_expansion/logs/ld_extract_%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

# Activate env
source activate /ix/djishnu/Priyamvada/envs/virauto/
export PYTHONNOUSERSITE=1
export JAVA_HOME="$CONDA_PREFIX"
export PATH="$JAVA_HOME/bin:$PATH"
export PYSPARK_PYTHON="$(which python)"
unset _JAVA_OPTIONS JAVA_TOOL_OPTIONS


cd /ix/djishnu/Priyamvada/virauto/analyses/ld_expansion/scripts/python
 python panukbb_ld_extract.py \
	    /ix/djishnu/Priyamvada/virauto/data/gwas_snps/p_val_filtered/MS/p_val_p1e-2/ms_p1e-2.chr10.txt.part00.txt \
		test.out
