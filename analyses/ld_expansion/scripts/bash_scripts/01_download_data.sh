#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 0-05:00
#SBATCH -J 1000_genomes_process
#SBATCH --output=download_data.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

# Load PLINK2 (okay if this is a normal bash run; sbatch will actually use it)
module purge 2>/dev/null || true
module load plink2/2.0.0a.6.9 2>/dev/null || true

set -euo pipefail

REPO_ROOT="/ix/djishnu/Priyamvada/virauto"

OUT_DIR="${REPO_ROOT}/data/refs/1000G"
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "[info] Repo root: ${REPO_ROOT}"
echo "[info] Saving 1000G files under: ${OUT_DIR}"

# ---------------------------------------------------------
# Download sources (replace URLs if you have a different host)
# ---------------------------------------------------------
PGEN_ZST_URL="https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst"
PVAR_ZST_URL="https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst"
PSAM_URL="https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam"

echo "[info] Downloading (resume safe)..."
wget -c -O all_phase3.pgen.zst "${PGEN_ZST_URL}"
wget -c -O all_phase3.pvar.zst "${PVAR_ZST_URL}"
wget -c -O phase3_corrected.psam "${PSAM_URL}"
echo "[ok] Downloads complete."

# ---------------------------------------------------------
# Decompress ZSTD files using PLINK2
# ---------------------------------------------------------
if [[ ! -s all_phase3.pgen ]]; then
  echo "[info] Decompressing all_phase3.pgen.zst -> all_phase3.pgen"
  plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
else
  echo "[skip] all_phase3.pgen already exists."
fi

if [[ ! -s all_phase3.pvar ]]; then
  echo "[info] Decompressing all_phase3.pvar.zst -> all_phase3.pvar"
  plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar
else
  echo "[skip] all_phase3.pvar already exists."
fi

# ---------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------
for f in all_phase3.pgen all_phase3.pvar phase3_corrected.psam; do
  [[ -s "$f" ]] || { echo "[error] Missing or empty: $f" >&2; exit 1; }
done

echo "[done] 1000G files ready in ${OUT_DIR}"

