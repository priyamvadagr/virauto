#!/bin/bash

# Script to submit netMHCpan jobs for chunk batches

CHUNK_LIST_DIR="/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunk_lists"
SLURM_SCRIPT="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/scripts/slurm/netmhcpan_chunks_from_list.sh"

# Check if chunk lists exist
if [ ! -d "$CHUNK_LIST_DIR" ]; then
    echo "ERROR: Chunk list directory not found: $CHUNK_LIST_DIR"
    echo "Please run find_unprocessed_chunks.sh first"
    exit 1
fi

# Find all batch files
BATCH_FILES=($(ls "$CHUNK_LIST_DIR"/chunk_batch_*.txt 2>/dev/null | sort -V))

if [ ${#BATCH_FILES[@]} -eq 0 ]; then
    echo "ERROR: No chunk batch files found in $CHUNK_LIST_DIR"
    echo "Please run find_unprocessed_chunks.sh first"
    exit 1
fi

echo "=========================================="
echo "  netMHCpan Chunk Processing Submission"
echo "=========================================="
echo ""
echo "Found ${#BATCH_FILES[@]} batch file(s):"

for batch_file in "${BATCH_FILES[@]}"; do
    batch_name=$(basename "$batch_file")
    num_chunks=$(wc -l < "$batch_file")
    echo "  - $batch_name ($num_chunks chunks)"
done

echo ""
read -p "Which batch do you want to submit? (1-${#BATCH_FILES[@]}, or 'all'): " choice

submit_batch() {
    local batch_file=$1
    local batch_name=$(basename "$batch_file" .txt)
    local num_chunks=$(wc -l < "$batch_file")
    
    echo ""
    echo "Preparing to submit: $batch_name"
    echo "  Chunks to process: $num_chunks"
    
    # Create a temporary SLURM script with the correct parameters
    TEMP_SCRIPT="/tmp/netmhcpan_${batch_name}_$$.sh"
    
    # Copy the template and update the array size and chunk list path
    sed "s|#SBATCH --array=.*|#SBATCH --array=1-${num_chunks}|" "$SLURM_SCRIPT" | \
    sed "s|CHUNK_LIST=.*|CHUNK_LIST=${batch_file}|" > "$TEMP_SCRIPT"
    
    echo "  Array range: 1-${num_chunks}"
    echo "  Chunk list: $batch_file"
    echo ""
    
    read -p "Submit this batch? (y/n): " confirm
    if [[ $confirm =~ ^[Yy]$ ]]; then
        job_id=$(sbatch "$TEMP_SCRIPT" | grep -oP '\d+')
        if [ $? -eq 0 ]; then
            echo "✓ Job submitted: $job_id"
            echo "  Monitor with: squeue -j $job_id"
            echo "  Cancel with: scancel $job_id"
        else
            echo "✗ Failed to submit job"
        fi
    else
        echo "Skipped"
    fi
    
    rm -f "$TEMP_SCRIPT"
}

if [[ $choice == "all" ]]; then
    for batch_file in "${BATCH_FILES[@]}"; do
        submit_batch "$batch_file"
    done
elif [[ $choice =~ ^[0-9]+$ ]] && [ $choice -ge 1 ] && [ $choice -le ${#BATCH_FILES[@]} ]; then
    batch_idx=$((choice - 1))
    submit_batch "${BATCH_FILES[$batch_idx]}"
else
    echo "Invalid choice"
    exit 1
fi

echo ""
echo "=========================================="
echo "Submission complete!"
echo "Check job status with: squeue -u $USER"
echo "=========================================="