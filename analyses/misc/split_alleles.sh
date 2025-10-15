#!/bin/bash

############################################################
# SPLIT ALLELE FILE INTO CHUNKS
# Divides 3,000 alleles into smaller chunks for better
# parallelization
############################################################

# --- Configuration ---
# UPDATE THIS PATH to your actual 3000-allele file
INPUT_ALLELE_FILE="/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_3000_alleles.txt"
OUTPUT_DIR="/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_chunks"
ALLELES_PER_CHUNK=30  # Adjust this: 30 = 100 chunks, 50 = 60 chunks, 100 = 30 chunks

echo "=========================================="
echo "  Allele File Splitter"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Input file: $INPUT_ALLELE_FILE"
echo "  Output dir: $OUTPUT_DIR"
echo "  Alleles per chunk: $ALLELES_PER_CHUNK"
echo ""

# Check input file exists
if [ ! -f "$INPUT_ALLELE_FILE" ]; then
    echo "ERROR: Allele file not found!"
    echo "  Looking for: $INPUT_ALLELE_FILE"
    echo ""
    echo "Please update INPUT_ALLELE_FILE in this script to point to your 3000-allele file."
    echo ""
    echo "Looking for possible allele files in the area..."
    find /ix/djishnu/Priyamvada/virauto/data/HLA_alleles/ -name "*.txt" -type f 2>/dev/null | head -10
    echo ""
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Count total alleles
total_alleles=$(wc -l < "$INPUT_ALLELE_FILE")
echo "Total alleles: $total_alleles"
echo "Alleles per chunk: $ALLELES_PER_CHUNK"

# Calculate number of chunks
num_chunks=$(( (total_alleles + ALLELES_PER_CHUNK - 1) / ALLELES_PER_CHUNK ))
echo "Number of chunks: $num_chunks"
echo ""

# Confirm
read -p "Split into $num_chunks chunks? (yes/no): " confirm
if [[ ! $confirm =~ ^[Yy][Ee][Ss]$ ]]; then
    echo "Aborted"
    exit 0
fi

echo ""
echo "Splitting..."

# Split the file
split -l $ALLELES_PER_CHUNK -d -a 3 "$INPUT_ALLELE_FILE" "$OUTPUT_DIR/allele_chunk_"

# Rename to add .txt extension
for file in "$OUTPUT_DIR"/allele_chunk_*; do
    if [[ ! "$file" =~ \.txt$ ]]; then
        mv "$file" "${file}.txt"
    fi
done

# Verify
actual_chunks=$(ls "$OUTPUT_DIR"/allele_chunk_*.txt 2>/dev/null | wc -l)

echo ""
echo "=========================================="
echo "  Complete!"
echo "=========================================="
echo ""
echo "Created $actual_chunks chunk files in:"
echo "  $OUTPUT_DIR"
echo ""
echo "Chunk files:"
ls "$OUTPUT_DIR"/allele_chunk_*.txt | head -10
if [ $actual_chunks -gt 10 ]; then
    echo "  ... and $((actual_chunks - 10)) more"
fi
echo ""

# Show sample from first chunk
echo "Sample from first chunk:"
head -3 "$OUTPUT_DIR"/allele_chunk_000.txt | sed 's/^/  /'
echo ""

# Create index file for reference
INDEX_FILE="$OUTPUT_DIR/allele_chunks_index.txt"
echo "# Allele chunk index" > "$INDEX_FILE"
echo "# Created: $(date)" >> "$INDEX_FILE"
echo "# Total alleles: $total_alleles" >> "$INDEX_FILE"
echo "# Alleles per chunk: $ALLELES_PER_CHUNK" >> "$INDEX_FILE"
echo "# Total chunks: $actual_chunks" >> "$INDEX_FILE"
echo "" >> "$INDEX_FILE"

chunk_num=0
for chunk_file in "$OUTPUT_DIR"/allele_chunk_*.txt; do
    chunk_alleles=$(wc -l < "$chunk_file")
    echo "chunk_${chunk_num}: $chunk_file ($chunk_alleles alleles)" >> "$INDEX_FILE"
    ((chunk_num++))
done

echo "Index file created: $INDEX_FILE"
echo ""