#!/bin/bash
#SBATCH -J copy_xls_files
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -t 0-04:00:00
#SBATCH --mem=4G
#SBATCH --cluster=htc
#SBATCH --output=copy_files_%j.log

echo "Starting file copy at $(date)"
echo "Host: $(hostname)"

SRC="/ix/djishnu/Tracy/AutoimmuneInfectious/results/netmhcpan/virscan/9_mers/Class_I/HLA-A/chunk_10"
DEST="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/Class_I/HLA-A/chunk_10"

mkdir -p "$DEST"

# Count files
total=$(find "$SRC" -maxdepth 1 -name '*.xls' -type f | wc -l)
echo "Total files to copy: $total"

# Use rsync for NFS reliability
rsync -av --progress \
    --include="*.xls" \
    --exclude="*" \
    "$SRC/" "$DEST/"

echo "Completed at $(date)"
echo "Files in destination: $(ls "$DEST"/*.xls 2>/dev/null | wc -l)"