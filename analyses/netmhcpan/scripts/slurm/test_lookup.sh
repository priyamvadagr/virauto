#!/bin/bash
#SBATCH -J netmhcpan_test
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-01:00:00
#SBATCH --mem=8G
#SBATCH --constraint=amd
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_lookup.err
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_lookup.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
##SBATCH --cluster=smp

OUTDIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/Type1_NR/All_types"
PEPTIDE_FILES=(/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks/*.fasta)
ALLELE_CHUNKS=(/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_NR_alleles/Type1_chunks/*.txt)

declare -A OUTPUT_MAP
processed_count=0
UNPROCESSED_COMBOS=()

start_time=$(date +%s)
echo "[INFO] Indexing output files..."

while IFS= read -r file; do
    base=$(basename "$file" .xls)
    OUTPUT_MAP["$base"]=1
done < <(find "$OUTDIR" -type f -name "*.xls")

echo "[INFO] Checking peptide × allele-chunk combinations..."
for pep_file in "${PEPTIDE_FILES[@]}"; do
    pep_name=$(basename "$pep_file" .fasta)

    for allele_chunk in "${ALLELE_CHUNKS[@]}"; do
        chunk_name=$(basename "$allele_chunk" .txt)
        allele_count=$(wc -l < "$allele_chunk")

        existing=0
        while read -r short long; do
            [[ -z "$short" || -z "$long" ]] && continue
            safe=$(echo "$long" | tr '*' '_' | tr ':' '_')
            key="${safe}_${pep_name}"
            if [[ -v OUTPUT_MAP[$key] ]]; then ((existing++)); fi
        done < "$allele_chunk"

        if (( existing < allele_count )); then
            UNPROCESSED_COMBOS+=("${pep_file}|${allele_chunk}")
        else
            ((processed_count++))
        fi
    done
done

echo "Already processed: $processed_count combinations"
echo "Remaining: ${#UNPROCESSED_COMBOS[@]} combinations"
if (( ${#UNPROCESSED_COMBOS[@]} == 0 )); then
    echo "✓ All combinations already processed!"
    exit 0
fi
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "[INFO] Completed in ${runtime} seconds (~$(echo "scale=2; $runtime/60" | bc) minutes)"
