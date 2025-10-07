#!/usr/bin/env python3
"""
===============================================================================
Volcano Plot of Mimicry Bias Across HLA Alleles (Type I: A + C)
===============================================================================
Each point = one HLA allele
X-axis: % Human-dominant (tolerance bias)
Y-axis: % Viral-dominant (infection mimicry bias)
Color = which side is enriched
Quantile thresholds (e.g. 90th percentile) used to highlight top alleles
===============================================================================
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

# =============================================================================
# 1. Input files and setup
# =============================================================================
base_dir = "/ix/djishnu/Priyamvada/virauto/results/virscan/mimicry_analysis/9_mers/data"
out_dir = "/ix/djishnu/Priyamvada/virauto/results/virscan/mimicry_analysis/9_mers/plots"

# Read multiple HLA allele summary files (A + C)
files = [
    os.path.join(base_dir, "allele_mimicry_summary_type1_A.csv"),
    os.path.join(base_dir, "allele_mimicry_summary_type1_B.csv"),
    os.path.join(base_dir, "allele_mimicry_summary_type1_C.csv"),
]

dfs = []
for f in files:
    if os.path.exists(f):
        df = pd.read_csv(f)
        # keep allele name from index or first column
        if "Allele" not in df.columns:
            df.rename(columns={df.columns[0]: "Allele"}, inplace=True)
        df["HLA_type"] = df["Allele"].str.extract(r"(HLA-[A-Z])")
        dfs.append(df)
    else:
        print(f"[WARN] Missing file: {f}")

# Merge datasets
allele_summary = pd.concat(dfs, ignore_index=True)
print(f"[INFO] Combined {len(allele_summary)} HLA alleles from {len(files)} files")

# =============================================================================
# 2. Quantile-based thresholds
# =============================================================================
q_hi_human = allele_summary["Human_dominant_pct"].quantile(0.9)
q_hi_viral = allele_summary["Viral_dominant_pct"].quantile(0.9)
print(f"[INFO] Quantile thresholds → Human ≥ {q_hi_human:.2f}%, Viral ≥ {q_hi_viral:.2f}%")

# Classification colors
def classify_color(row):
    if row["Human_dominant_pct"] >= q_hi_human and row["Viral_dominant_pct"] >= q_hi_viral:
        return "#8e44ad"  # dual high
    elif row["Human_dominant_pct"] >= q_hi_human:
        return "#27ae60"  # human-dominant (tolerance mimicry)
    elif row["Viral_dominant_pct"] >= q_hi_viral:
        return "#c0392b"  # viral-dominant (infection mimicry)
    else:
        return "#95a5a6"  # neutral
    
allele_summary_copy = allele_summary.copy()
allele_summary_copy["candidate"] = allele_summary_copy["Human_dominant_pct"] >= q_hi_human

# Filter to only candidate alleles and save
candidate_alleles = allele_summary_copy[allele_summary_copy["candidate"]]
candidate_alleles.to_csv(os.path.join(base_dir, "human_dominant_candidate_alleles.csv"), index=True)

allele_summary["color"] = allele_summary.apply(classify_color, axis=1)
allele_summary["log_pairs"] = np.log10(allele_summary["Total"] + 1)
allele_summary["size"] = 40 + 60 * (allele_summary["log_pairs"] / allele_summary["log_pairs"].max())

# =============================================================================
# 3. Volcano Plot
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

# Label top alleles with adjustText to avoid overlap
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

# Adjust text positions to avoid overlap
adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
    expand_points=(1.5, 1.5),
    expand_text=(1.2, 1.2),
    force_points=(0.5, 0.5),
    force_text=(0.5, 0.5)
)

# Axes & legend
plt.xlabel("% Human-dominant peptide pairs (tolerogenic mimicry bias)", fontsize=11)
plt.ylabel("% Viral-dominant peptide pairs (infection activation bias)", fontsize=11)
plt.title("Per-HLA Mimicry Volcano Plot (Type I: A + C)", fontsize=14, weight="bold")
plt.legend(frameon=False, fontsize=8, loc="upper right")
plt.grid(alpha=0.25)
plt.tight_layout()

# Save
os.makedirs(out_dir, exist_ok=True)
plot_path = os.path.join(out_dir, "volcano_HLA_A_B_C_combined.png")
plt.savefig(plot_path, dpi=300)
plt.show()
print(f"✅ Volcano plot saved → {plot_path}")