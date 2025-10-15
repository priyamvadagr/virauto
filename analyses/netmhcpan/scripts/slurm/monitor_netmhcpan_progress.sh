#!/bin/bash

############################################################
# NETMHCPAN PROGRESS MONITOR (Updated for Type1_NR)
# Tracks completion status of automated batch jobs
############################################################

# === PATHS (update if you change your submission script) ===
OUTDIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/Type1_NR/All_types_chunks/chunk_2"
PEPTIDE_DIR="/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/split_chunks/chunk_2"
ALLELE_CHUNK_DIR="/ix/djishnu/Priyamvada/virauto/data/HLA_alleles/Type1_NR_alleles/Type1_chunks"
LOG_DIR="/ix/djishnu/Priyamvada/virauto/analyses/netmhcpan/logs/Type1_NR_chunks/tmp"

# === Colors ===
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

############################################################
# COLLECT STATISTICS
############################################################

echo "=========================================="
echo "  NetMHCpan Progress Monitor (Type1_NR)"
echo "=========================================="
echo ""

# Count peptide and allele chunks
total_peptide_files=$(find "$PEPTIDE_DIR" -name "*.fasta" -type f | wc -l)
total_allele_chunks=$(find "$ALLELE_CHUNK_DIR" -name "allele_chunk_*.txt" -type f | wc -l)

# Estimate total expected .xls files
# Average alleles per chunk (approx)
if [ $total_allele_chunks -gt 0 ]; then
    avg_alleles_per_chunk=$(awk '{ total += NF } END { print (NR ? total/NR : 0) }' "$ALLELE_CHUNK_DIR"/allele_chunk_*.txt 2>/dev/null)
else
    avg_alleles_per_chunk=0
fi

total_expected_estimate=$(awk "BEGIN {printf \"%d\", $total_peptide_files * $total_allele_chunks * $avg_alleles_per_chunk / $total_allele_chunks}")

echo "Expected outputs (approx):"
echo "  Peptide chunks: $total_peptide_files"
echo "  Allele chunks:  $total_allele_chunks (~30 alleles each)"
echo "  Estimated total .xls files: ≈ $((total_peptide_files * total_allele_chunks * 30))"
echo ""

# Count completed outputs
if [ -d "$OUTDIR" ]; then
    completed_xls=$(find "$OUTDIR" -name "*.xls" -type f | wc -l)
    completed_txt=$(find "$OUTDIR" -name "*.txt" -type f | wc -l)
else
    completed_xls=0
    completed_txt=0
fi

# Estimate percentage
total_expected=$((total_peptide_files * total_allele_chunks * 30))
if [ $total_expected -gt 0 ]; then
    percent_complete=$(awk "BEGIN {printf \"%.2f\", ($completed_xls/$total_expected)*100}")
else
    percent_complete=0
fi

echo "Completed outputs:"
echo "  .xls files: $completed_xls / ~${total_expected} (${percent_complete}%)"
echo "  .txt files: $completed_txt"
echo ""

# Progress bar
bar_length=50
filled_length=$(awk "BEGIN {printf \"%.0f\", ($percent_complete/100)*$bar_length}")
bar=$(printf "%${filled_length}s" | tr ' ' '█')
empty=$(printf "%$((bar_length - filled_length))s" | tr ' ' '░')
echo -e "Progress: [${GREEN}${bar}${NC}${empty}] ${percent_complete}%"
echo ""

############################################################
# JOB STATUS
############################################################

echo "=========================================="
echo "  Job Status"
echo "=========================================="
echo ""
USER=prg65
job_count=$(squeue -u "$USER" -h -t PENDING,RUNNING --name=netmhcpan_auto -M htc | wc -l)

if [ $job_count -eq 0 ]; then
    echo -e "${YELLOW}No jobs currently running or pending${NC}"
else
    echo "Active netMHCpan jobs:"
    squeue -u "$USER" -M htc -t PENDING,RUNNING -o "  %A %t %D %C %m %l %S %j" --name=netmhcpan_auto 2>/dev/null | head -20
    if [ $job_count -gt 20 ]; then
        echo "  ... and $((job_count - 20)) more"
    fi
fi
echo ""

############################################################
# RECENT ERRORS
############################################################

echo "=========================================="
echo "  Recent Errors (last 10)"
echo "=========================================="
echo ""

if [ -d "$LOG_DIR" ]; then
    recent_errors=$(find "$LOG_DIR" -name "*.err" -type f -mtime -2 ! -empty 2>/dev/null | head -10)
    if [ -z "$recent_errors" ]; then
        echo -e "${GREEN}✓ No recent errors found${NC}"
    else
        echo -e "${RED}Found error logs:${NC}"
        echo "$recent_errors" | while read -r err_file; do
            echo "  - $(basename "$err_file")"
            echo "    Last line: $(tail -1 "$err_file")"
        done
    fi
else
    echo "Log directory not found: $LOG_DIR"
fi
echo ""

############################################################
# PER-CHUNK COMPLETION STATUS
############################################################

if [[ "$1" == "--detailed" || "$1" == "-d" ]]; then
    echo "=========================================="
    echo "  Per-Chunk Completion Status"
    echo "=========================================="
    echo ""

    find "$PEPTIDE_DIR" -name "*.fasta" -type f | sort | while read -r pep_file; do
        chunk_name=$(basename "$pep_file" .fasta)
        chunk_outputs=$(find "$OUTDIR" -name "*_${chunk_name}.xls" 2>/dev/null | wc -l)

        # Roughly assume ~30 alleles per chunk
        expected_chunk_outputs=30

        if [ $chunk_outputs -ge $expected_chunk_outputs ]; then
            status="${GREEN}✓${NC}"
        elif [ $chunk_outputs -gt 0 ]; then
            status="${YELLOW}◐${NC}"
        else
            status="${RED}○${NC}"
        fi

        echo -e "  $status $chunk_name: $chunk_outputs/${expected_chunk_outputs}"
    done | head -50

    echo ""
    echo "Legend: ${GREEN}✓${NC}=Complete  ${YELLOW}◐${NC}=Partial  ${RED}○${NC}=Pending"
    echo "(Showing first 50 chunks; use full find for complete listing)"
    echo ""
fi

############################################################
# SUMMARY
############################################################

echo "=========================================="
echo "  Summary"
echo "=========================================="
echo ""

if [ $completed_xls -gt 0 ] && [ $completed_xls -ge $((total_expected * 9 / 10)) ]; then
    echo -e "${GREEN}✓ Nearly all or all predictions complete!${NC}"
    echo ""
    echo "Next steps:"
    echo "  - Verify results integrity"
    echo "  - Merge .xls outputs if needed"
elif [ $job_count -eq 0 ] && [ $completed_xls -lt $((total_expected * 9 / 10)) ]; then
    echo -e "${RED}⚠ Jobs finished but outputs incomplete${NC}"
    echo ""
    echo "Possible issues:"
    echo "  - Some batches failed"
    echo "  - Check logs: $LOG_DIR/*.err"
    echo "  - Re-run auto_submit.sh (it skips completed chunks)"
else
    echo -e "${BLUE}Jobs running...${NC}"
    echo ""
    echo "Monitor live with:"
    echo "  watch -n 60 '$0'"
    echo ""
    echo "Check queue:"
    echo "  squeue -u $USER --name=netmhcpan_auto"
fi

echo ""
echo "=========================================="
