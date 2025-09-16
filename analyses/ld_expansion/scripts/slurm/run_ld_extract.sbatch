#!/bin/bash
#SBATCH -J ld_extract
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/ld_expansion/logs/ld_extract_%x_%j.out
#SBATCH --error=/ix/djishnu/Priyamvada/virauto/analyses/ld_expansion/logs/ld_extract_%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

# Usage: sbatch run_ld_extract.sbatch <DISEASE>
# Example: sbatch run_ld_extract.sbatch MS

if [ $# -lt 1 ]; then
  echo "Usage: sbatch run_ld_extract.sbatch <DISEASE>"
  exit 1
fi

DISEASE="$1"
echo "[info] Processing disease: $DISEASE"

# Activate env
source activate /ix/djishnu/Priyamvada/envs/virauto/

# ----- anchor to repo root so relative paths are correct -----
cd /ix/djishnu/Priyamvada/virauto
echo "[info] Working directory: $(pwd)"

# Input and output roots (relative to repo root)
IN_DIR="data/gwas_snps/p_val_filtered/${DISEASE}"
OUT_BASE="results/ld_expansion/${DISEASE}"
mkdir -p "$OUT_BASE"

# For each p-value directory (e.g., p_val_p1e-3)
for pdir in "$IN_DIR"/p_val_*; do
  [ -d "$pdir" ] || continue
  pname="$(basename "$pdir")"            # e.g., p_val_p1e-3
  OUT_DIR="${OUT_BASE}/${pname}"
  mkdir -p "$OUT_DIR"
  echo "[info] Output directory: $OUT_DIR"

  # Run the extractor for every SNP list in this p-val directory
  for f in "$pdir"/*.txt; do
    [ -f "$f" ] || continue
    base="$(basename "$f" .txt)"
    echo "[run] $f -> ${OUT_DIR}/${base}"
    python analyses/ld_expansion/scripts/python/panukbb_ld_extract.py \
      "$f" \
      "${OUT_DIR}/${base}"
  done
done

echo "[done] Finished $DISEASE"

