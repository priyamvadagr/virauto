#!/bin/bash
#SBATCH -J netmhcpan_check
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 0-00:30:00
#SBATCH --mem=8G
#SBATCH --constraint=amd
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/robust_fast_lookup.err
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/robust_fast_lookup.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=htc

############################################################
# ROBUST FAST VERSION
# Correctly handles allele chunk files by checking actual
# allele names from within the files
############################################################

OUTDIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/Type1_NR/All_types"
PEPTIDE_DIR="/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks"
ALLELE_CHUNK_DIR="/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_NR_alleles/Type1_chunks"

echo "=========================================="
echo "  Robust Fast Checker"
echo "=========================================="
echo ""

start_time=$(date +%s)

# Count inputs
num_peptides=$(ls -1 "$PEPTIDE_DIR"/*.fasta 2>/dev/null | wc -l)
num_allele_chunks=$(ls -1 "$ALLELE_CHUNK_DIR"/*.txt 2>/dev/null | wc -l)
total_combos=$((num_peptides * num_allele_chunks))

echo "[INFO] Configuration:"
echo "  Peptide chunks: $num_peptides"
echo "  Allele chunks: $num_allele_chunks"
echo "  Total combinations: $total_combos"
echo ""

############################################################
# STEP 1: BUILD SET OF EXISTING OUTPUTS
############################################################
echo "[1/2] Indexing existing outputs..."
index_start=$(date +%s)

TMP_DIR="/tmp/netmhcpan_check_$$"
mkdir -p "$TMP_DIR"

# Get all existing .xls filenames (basename only)
# This is much faster than repeated find operations
ls -1 "$OUTDIR"/*.xls 2>/dev/null | xargs -n 1 basename > "$TMP_DIR/existing_outputs.txt"

# Convert to associative array for O(1) lookup
declare -A existing_outputs
while read -r filename; do
    existing_outputs["$filename"]=1
done < "$TMP_DIR/existing_outputs.txt"

num_existing=${#existing_outputs[@]}

index_end=$(date +%s)
index_time=$((index_end - index_start))

echo "  ✓ Found $num_existing existing output files in ${index_time}s"
echo ""

############################################################
# STEP 2: CHECK COMBINATIONS
############################################################
echo "[2/2] Checking combinations..."
check_start=$(date +%s)

UNPROCESSED_FILE="$TMP_DIR/unprocessed.txt"
> "$UNPROCESSED_FILE"

processed_count=0
unprocessed_count=0

combo_num=0

# For each peptide file
for pep_file in "$PEPTIDE_DIR"/*.fasta; do
    [ -f "$pep_file" ] || continue
    
    pep_name=$(basename "$pep_file" .fasta)
    
    # Progress indicator
    ((combo_num++))
    if (( combo_num % 50 == 0 )); then
        echo "  Progress: $combo_num/$num_peptides peptides..."
    fi
    
    # For each allele chunk
    for allele_file in "$ALLELE_CHUNK_DIR"/*.txt; do
        [ -f "$allele_file" ] || continue
        
        # Check if ALL alleles in this chunk have outputs for this peptide
        all_complete=true
        missing_count=0
        checked_count=0
        
        # Read each allele from the chunk file and check if output exists
        while read -r short_allele long_allele; do
            [ -z "$short_allele" ] || [ -z "$long_allele" ] && continue
            
            ((checked_count++))
            
            # Convert allele name to safe filename (same as NetMHCpan output)
            safe_allele=$(echo "$long_allele" | tr '*' '_' | tr ':' '_')
            expected_filename="${safe_allele}_${pep_name}.xls"
            
            # O(1) lookup in hash table
            if [ -z "${existing_outputs[$expected_filename]}" ]; then
                all_complete=false
                ((missing_count++))
                # Don't break - we want to count all missing for logging if needed
            fi
            
            # Early exit optimization: if we find any missing, mark as unprocessed
            if [ "$all_complete" = false ] && [ $missing_count -ge 1 ]; then
                break
            fi
        done < "$allele_file"
        
        # Record result
        if [ "$all_complete" = false ]; then
            echo "${pep_file}|${allele_file}" >> "$UNPROCESSED_FILE"
            ((unprocessed_count++))
        else
            ((processed_count++))
        fi
    done
done

check_end=$(date +%s)
check_time=$((check_end - check_start))

echo "  ✓ Checked $total_combos combinations in ${check_time}s"
echo ""

############################################################
# SUMMARY
############################################################
end_time=$(date +%s)
total_runtime=$((end_time - start_time))

echo "=========================================="
echo "  Results"
echo "=========================================="
echo "  Total combinations: $total_combos"
echo "  Already processed: $processed_count"
echo "  Unprocessed: $unprocessed_count"
echo "  Completion: $(awk "BEGIN {printf \"%.1f\", ($processed_count/$total_combos)*100}")%"
echo ""
echo "=========================================="
echo "  Performance"
echo "=========================================="
echo "  Indexing: ${index_time}s"
echo "  Checking: ${check_time}s"
echo "  Total: ${total_runtime}s (~$(awk "BEGIN {printf \"%.1f\", $total_runtime/60}") min)"
echo "  Speedup vs original: ~$(awk "BEGIN {printf \"%.0f\", 1500/$total_runtime}")x"
echo ""

if [ $unprocessed_count -eq 0 ]; then
    echo "✓ All combinations already processed!"
    rm -rf "$TMP_DIR"
    exit 0
fi

echo "=========================================="
echo "  Next Steps"
echo "=========================================="
echo "Unprocessed combinations saved to:"
echo "  $UNPROCESSED_FILE"
echo ""
echo "Number of lines: $(wc -l < "$UNPROCESSED_FILE")"
echo ""
echo "First 5 unprocessed combinations:"
head -5 "$UNPROCESSED_FILE"
echo ""
echo "Use this file for batch submission!"
echo "=========================================="