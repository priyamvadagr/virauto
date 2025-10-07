#!/bin/bash

# Directories
PROCESSED_DIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_A_chunks"
CHUNK_DIR="/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks"
OUTPUT_DIR="/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunk_lists"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Extract unique chunk numbers that have been processed
echo "Finding processed chunks..."
PROCESSED_CHUNKS=$(ls "$PROCESSED_DIR" | grep -oP 'chunk_\K\d+' | sort -u)

# Save processed chunks to a file
echo "$PROCESSED_CHUNKS" > "$OUTPUT_DIR/processed_chunks.txt"
NUM_PROCESSED=$(echo "$PROCESSED_CHUNKS" | wc -l)
echo "Found $NUM_PROCESSED processed chunks"

# Get all available chunk files
echo "Finding all available chunk files..."
ALL_CHUNKS=$(ls "$CHUNK_DIR"/matched_pairs_chunk_*.fasta | grep -oP 'chunk_\K\d+' | sort -u)

# Save all chunks to a file
echo "$ALL_CHUNKS" > "$OUTPUT_DIR/all_chunks.txt"
NUM_ALL=$(echo "$ALL_CHUNKS" | wc -l)
echo "Found $NUM_ALL total chunks"

# Find unprocessed chunks (chunks in ALL but not in PROCESSED)
echo "Identifying unprocessed chunks..."
UNPROCESSED_CHUNKS=$(comm -23 \
    <(echo "$ALL_CHUNKS" | sort -n) \
    <(echo "$PROCESSED_CHUNKS" | sort -n))

# Save unprocessed chunks to a file
echo "$UNPROCESSED_CHUNKS" > "$OUTPUT_DIR/unprocessed_chunks.txt"
NUM_UNPROCESSED=$(echo "$UNPROCESSED_CHUNKS" | wc -l)
echo "Found $NUM_UNPROCESSED unprocessed chunks"

# Create file lists with full paths, 100 chunks per file
echo "Creating file lists (100 chunks per file)..."

BATCH_SIZE=100
batch_num=1
count=0
current_file="$OUTPUT_DIR/chunk_batch_${batch_num}.txt"

# Clear/create first batch file
> "$current_file"

for chunk in $UNPROCESSED_CHUNKS; do
    # Add the full path to the chunk file
    echo "$CHUNK_DIR/matched_pairs_chunk_${chunk}.fasta" >> "$current_file"
    
    ((count++))
    
    # If we've hit the batch size, start a new file
    if [ $count -eq $BATCH_SIZE ]; then
        echo "Created $current_file (100 chunks)"
        ((batch_num++))
        current_file="$OUTPUT_DIR/chunk_batch_${batch_num}.txt"
        > "$current_file"
        count=0
    fi
done

# Report on the last batch if it has any entries
if [ $count -gt 0 ]; then
    echo "Created $current_file ($count chunks)"
fi

# Summary
echo ""
echo "=========================================="
echo "              SUMMARY"
echo "=========================================="
echo "Total chunks available: $NUM_ALL"
echo "Processed chunks: $NUM_PROCESSED"
echo "Unprocessed chunks: $NUM_UNPROCESSED"
echo "Batch files created: $batch_num"
echo ""
echo "Output files saved in: $OUTPUT_DIR"
echo "  - all_chunks.txt: all available chunks"
echo "  - processed_chunks.txt: already processed"
echo "  - unprocessed_chunks.txt: need to process"
echo "  - chunk_batch_*.txt: batched file lists"
echo "=========================================="