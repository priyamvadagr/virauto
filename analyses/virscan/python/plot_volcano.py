#!/usr/bin/env python3
#!/usr/bin/env python3
"""
===============================================================================
Script: plot_mimicry_volcano_by_allele_types.py
Author: Priyamvada Guha Roy
===============================================================================
Description:
    This script generates a volcano-style scatter plot summarizing viral–
    human peptide mimicry bias across HLA alleles. Each allele is represented
    as a point positioned by:
        - X-axis: % Human-dominant peptide pairs (tolerogenic mimicry bias)
        - Y-axis: % Viral-dominant peptide pairs (infection activation bias)
    Points are color-coded based on their classification:
        - Green:  Human-dominant (tolerance mimicry)
        - Red:    Viral-dominant (infection mimicry)
        - Purple: Dual-high (both dominant)
        - Gray:   Neutral (below thresholds)
    
    The script automatically detects available HLA allele classes (e.g., A, B, C)
    from the filenames of the input CSVs (named as:
        allele_mimicry_summary_{type_}_{class_}.csv).
    Quantile-based thresholds (90th percentile) are computed dynamically
    across all included alleles to highlight dominant subsets.

Enhancements:
    • Automatically identifies and merges all available HLA allele types.
    • Dynamically adjusts the plot title and output filename to reflect
      the specific HLA allele types present (e.g., A, B, C or A & C).
    • Applies adaptive point scaling based on the log10 of the total
      number of peptide pairs.
    • Annotates top alleles above either threshold using text adjustment
      to minimize overlap.
    • Saves a high-resolution figure (PNG) with consistent formatting.

Inputs:
    - allele_mimicry_summary_{type_}_{class_}.csv
      (one or more per HLA allele type)

Outputs:
    - Volcano plot (PNG) with title dynamically reflecting allele types
      e.g., volcano_type1_HLA_A_C.png
    - Console summary of included files, thresholds, and allele counts

Usage Example:
    python plot_mimicry_volcano_by_allele_types.py

Example Output Filenames:
    volcano_type1_HLA_A.png
    volcano_type1_HLA_A_C.png
    volcano_type1_HLA_A_B_C.png

===============================================================================
"""


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import glob

# =============================================================================
# 1. Input files and setup
# =============================================================================
data_dir = "/ix/djishnu/Priyamvada/virauto/results/virscan/mimicry_analysis/viral_strong_binders/9_mers/type1/data"
out_dir = "/ix/djishnu/Priyamvada/virauto/results/virscan/mimicry_analysis/viral_strong_binders/9_mers/type1/plots"
os.makedirs(out_dir, exist_ok=True)

# Identify all allele summary files (e.g., allele_mimicry_summary_type1_A.csv)
files = glob.glob(os.path.join(data_dir, "allele_mimicry_summary_*_*.csv"))
if not files:
    raise FileNotFoundError(f"No allele_mimicry_summary_*_*.csv files found in {data_dir}")

# Extract unique allele types (A, B, C, etc.)
allele_types = sorted({os.path.basename(f).split("_")[-1].replace(".csv", "") for f in files})
print(f"[INFO] Found allele types: {', '.join(allele_types)}")

# Read and merge all allele summaries
dfs = []
for f in files:
    df = pd.read_csv(f)
    if "Allele" not in df.columns:
        df.rename(columns={df.columns[0]: "Allele"}, inplace=True)
    df["HLA_type"] = df["Allele"].str.extract(r"(HLA-[A-Z])")
    dfs.append(df)

allele_summary = pd.concat(dfs, ignore_index=True)
print(f"[INFO] Combined {len(allele_summary)} alleles from {len(files)} files")

# =============================================================================
# 2. Quantile-based thresholds
# =============================================================================
q_hi_human = allele_summary["Human_dominant_pct"].quantile(0.9)
q_hi_viral = allele_summary["Viral_dominant_pct"].quantile(0.9)
print(f"[INFO] Quantile thresholds → Human ≥ {q_hi_human:.2f}%, Viral ≥ {q_hi_viral:.2f}%")

# ===================================================
# 3. Annotate Alleles with Dominance Classification
# ===================================================

def classify_bias(row):
    """Assign categorical label based on quantile thresholds."""
    if row["Human_dominant_pct"] >= q_hi_human and row["Viral_dominant_pct"] >= q_hi_viral:
        return "Dual_high"
    elif row["Human_dominant_pct"] >= q_hi_human:
        return "Human_dominant"
    elif row["Viral_dominant_pct"] >= q_hi_viral:
        return "Viral_dominant"
    else:
        return "Neutral"

# Apply classification
allele_summary["Dominance_category"] = allele_summary.apply(classify_bias, axis=1)

# ====================================================================================================
# 4. Assign Colors Consistent with Classification &
#    Compute log-scaled point sizes for visualization
# =====================================================================================================
def classify_color(row):
    if row["Dominance_category"] == "Dual_high":
        return "#8e44ad"  # purple
    elif row["Dominance_category"] == "Human_dominant":
        return "#27ae60"  # green
    elif row["Dominance_category"] == "Viral_dominant":
        return "#c0392b"  # red
    else:
        return "#95a5a6"  # gray

allele_summary["color"] = allele_summary.apply(classify_color, axis=1)
allele_summary["log_pairs"] = np.log10(allele_summary["Total"] + 1)
allele_summary["size"] = 40 + 60 * (allele_summary["log_pairs"] / allele_summary["log_pairs"].max())

# ===================================================
# 5. Save Annotated Summary (without log/size columns)
# ===================================================
cols_to_save = [c for c in allele_summary.columns if c not in ["log_pairs", "size", "color"]]
allele_summary[cols_to_save].to_csv(
    os.path.join(data_dir, "allele_summary_with_candidates.csv"), index=False
)

# =============================================================================
# 6. Volcano Plot
# =============================================================================
plt.figure(figsize=(10, 8))
sns.set(style="whitegrid")

plt.scatter(
    allele_summary["Human_dominant_pct"],
    allele_summary["Viral_dominant_pct"],
    s=allele_summary["size"],
    c=allele_summary["color"],
    alpha=0.85,
    edgecolors="black",
    linewidth=0.3
)

# Threshold lines
plt.axvline(q_hi_human, color="#27ae60", linestyle="--", lw=1.2, label=f"Human ≥ {q_hi_human:.1f}%")
plt.axhline(q_hi_viral, color="#c0392b", linestyle="--", lw=1.2, label=f"Viral ≥ {q_hi_viral:.1f}%")

# Label top alleles with adjustText
top_alleles = allele_summary[
    (allele_summary["Human_dominant_pct"] >= q_hi_human) |
    (allele_summary["Viral_dominant_pct"] >= q_hi_viral)
]

texts = []
for _, row in top_alleles.iterrows():
    texts.append(
        plt.text(
            row["Human_dominant_pct"],
            row["Viral_dominant_pct"],
            row["Allele"],
            fontsize=7,
            weight="bold",
            color=row["color"]
        )
    )

adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
    expand_points=(1.5, 1.5),
    expand_text=(1.2, 1.2),
    force_points=(0.5, 0.5),
    force_text=(0.5, 0.5)
)

# =============================================================================
# 5. Dynamic title and filename
# =============================================================================
# Make title context-aware based on number of allele types
if len(allele_types) == 1:
    type_label = f"HLA-{allele_types[0]}"
elif len(allele_types) == 2:
    type_label = f"HLA-{allele_types[0]} & HLA-{allele_types[1]}"
else:
    type_label = f"HLA-{', '.join(allele_types)}"

plt.xlabel("% Human-dominant peptide pairs (tolerogenic mimicry bias)", fontsize=11)
plt.ylabel("% Viral-dominant peptide pairs (infection activation bias)", fontsize=11)
plt.title(f"Mimicry Volcano Plot ({type_label})", fontsize=14, weight="bold")
plt.legend(frameon=False, fontsize=8, loc="upper right")
plt.grid(alpha=0.25)
plt.tight_layout()

# Construct dynamic filename
allele_str = "_".join(allele_types)
plot_filename = f"volcano_type1_HLA_{allele_str}.png"
plot_path = os.path.join(out_dir, plot_filename)

# Save
plt.savefig(plot_path, dpi=300)
plt.show()
print(f"✅ Volcano plot saved → {plot_path}")
