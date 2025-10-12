#!/bin/bash
#SBATCH -J netmhcpan_test
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-01:00:00
#SBATCH --mem=8G
#SBATCH --constraint=amd
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_lookup_clause.err
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/test_lookup_claude.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=htc

OUTDIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/Type1_NR/All_types"
PEPTIDE_FILES=(/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks/*.fasta)
echo "[INFO] First few FASTA files detected:"
ls -1 "${PEPTIDE_FILES[@]:0:3}"
ALLELE_CHUNKS=(/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_NR_alleles/Type1_chunks/*.txt)

start_time=$(date +%s)
echo "[INFO] Indexing output files... at $start_time"

UNPROCESSED_COMBOS=()
processed_count=0

# For each combination, check if all expected outputs exist
for pep_file in "${PEPTIDE_FILES[@]}"; do
    pep_name=$(basename "$pep_file" .fasta)
    echo $pep_name
    
    for allele_chunk in "${ALLELE_CHUNKS[@]}"; do
        chunk_name=$(basename "$allele_chunk" .txt)
        allele_count=$(wc -l < "$allele_chunk")
        
        # Count how many outputs exist for this combination
        existing_outputs=$(find "$OUTDIR" -name "*_${pep_name}.xls" -type f 2>/dev/null | \
                          xargs grep -l "$(head -1 "$allele_chunk" | awk '{print $1}')" 2>/dev/null | wc -l)
        
        # Simple heuristic: if any output exists for this combo, check more carefully
        # For speed, we'll just check if the expected number of files exist
        outputs_for_combo=$(find "$OUTDIR" -name "*_${pep_name}.xls" -type f | wc -l)
        
        if [ $outputs_for_combo -lt $allele_count ]; then
            UNPROCESSED_COMBOS+=("${pep_file}|${allele_chunk}")
        else
            ((processed_count++))
        fi
    done
done

echo "  Already processed: $processed_count combinations"
echo "  Remaining: ${#UNPROCESSED_COMBOS[@]} combinations"
echo ""

if [ ${#UNPROCESSED_COMBOS[@]} -eq 0 ]; then
    echo "✓ All combinations already processed!"
    exit 0
fi
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "[INFO] Completed in ${runtime} seconds (~$(echo "scale=2; $runtime/60" | bc) minutes)"