#!/usr/bin/env bash
# Split GWAS SNP lists by chromosome into per-threshold directories
# and chunk large chromosome files to manageable sizes.
#
# Input format per line: CHR:POS:REF:ALT (GRCh37)
#
# Usage:
#   bash scripts/misc/split_by_chr_and_chunk.sh [--delete] [--chunk-size N] [--base PATH]
#
# Defaults:
#   --base        /ix/djishnu/Priyamvada/virauto/data/gwas_snps/p_val_filtered
#   --chunk-size  3000    (max lines per chunk)
#
# Example:
#   bash scripts/misc/split_by_chr_and_chunk.sh --delete --chunk-size 2500

set -euo pipefail

# -------------------- defaults --------------------
BASE="/ix/djishnu/Priyamvada/virauto/data/gwas_snps/p_val_filtered"
CHUNK_SIZE=3000
DELETE=false

# -------------------- args ------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base)       BASE="$2"; shift 2 ;;
    --chunk-size) CHUNK_SIZE="$2"; shift 2 ;;
    --delete)     DELETE=true; shift ;;
    -h|--help)
      grep '^# ' "$0" | cut -c3-
      exit 0
      ;;
    *)
      echo "[error] Unknown arg: $1" >&2; exit 1 ;;
  esac
done

echo "[info] BASE       = $BASE"
echo "[info] CHUNK_SIZE = $CHUNK_SIZE"
echo "[info] DELETE     = $DELETE"

# -------------------- helpers ---------------------
chunk_file () {
  local file="$1"
  local n
  n=$(wc -l < "$file")
  if [[ "$n" -le "$CHUNK_SIZE" ]]; then
    echo "    [chunk] OK: $file ($n lines) <= $CHUNK_SIZE"
    return 0
  fi

  echo "    [chunk] Splitting $file ($n lines) into ~${CHUNK_SIZE}-line parts"
  # split writes file.part00.txt, file.part01.txt, ...
  split -l "$CHUNK_SIZE" -d --additional-suffix=.txt "$file" "${file}.part"
  # Optional: remove the large original chr file after chunking
  rm "$file"
}

# -------------------- main ------------------------
# Expect trait subdirs, e.g., MS/, SLE/, T1D/
for trait_dir in "$BASE"/*; do
  [[ -d "$trait_dir" ]] || continue
  trait=$(basename "$trait_dir")
  echo "[info] Trait: $trait"

  # For each thresholded list like sle_p1e-3.all.txt
  shopt -s nullglob
  for infile in "$trait_dir"/*.all.txt; do
    [[ -f "$infile" ]] || continue

    fname=$(basename "$infile" .all.txt)  # e.g., sle_p1e-3
    # threshold string is everything after the first underscore, e.g., p1e-3
    thresh="${fname#*_}"                  # p1e-3
    outdir="$trait_dir/p_val_${thresh}"   # p_val_p1e-3

    mkdir -p "$outdir"
    echo "  [split] $infile -> $outdir/${fname}.chr*.txt"

    # Split by chromosome based on field before the first colon
    awk -F: -v base="$outdir/${fname}" '{ print > base ".chr" $1 ".txt" }' "$infile"

    # Optionally delete original *.all.txt
    if $DELETE; then
      echo "  [delete] $infile"
      rm "$infile"
    fi

    # Chunk any large per-chromosome files in this threshold dir
    for chrfile in "$outdir/${fname}.chr"*.txt; do
      [[ -f "$chrfile" ]] || continue
      chunk_file "$chrfile"
    done

  done
  shopt -u nullglob
done

echo "[done] Split by chromosome and chunked large files."

