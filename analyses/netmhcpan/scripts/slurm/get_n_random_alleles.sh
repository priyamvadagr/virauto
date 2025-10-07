#!/bin/bash
cd /ix/djishnu/Priyamvada/auto_immune/NetMHCpan/NetMHCpan/netMHCpan-4.2/data
# Make output dir
mkdir -p HLA_random200

# Randomly sample 200 alleles from each class
grep "^HLA-A" HLA_allele_list | shuf -n 200 > /ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_random_sampling/HLA_A_random200.txt
grep "^HLA-B" HLA_allele_list | shuf -n 200 > /ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_random_sampling/HLA_B_random200.txt
grep "^HLA-C" HLA_allele_list | shuf -n 200 > /ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_random_sampling/HLA_C_random200.txt

