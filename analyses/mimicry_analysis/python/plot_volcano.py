#!/usr/bin/env python3
"""
===============================================================================
Script: plot_mimicry_volcano_from_parquet.py
Author: Priyamvada Guha Roy
===============================================================================
Description:
    Generates a volcano-style scatter plot summarizing viral–human peptide
    mimicry bias across HLA alleles, using the long-format Parquet dataset
    partitioned by allele.

Inputs:
    Parquet dataset directory containing:
        allele=HLA-A_0101/part-*.parquet
        allele=HLA-B_0702/part-*.parquet
        etc.

Outputs:
    volcano_HLA_[A|B|C|...].png — Volcano-style mimicry plot
===============================================================================
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import pyarrow.parquet as pq
import re

# =============================================================================
# 1. Configuration
# =============================================================================
BASE_PATH = "/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/virscan/Type1_NR_all/delta_binding/9_mers/"
DATASET_PATH = os.path.join(BASE_PATH, "data")
OUT_DATA_DIR = os.path.join(BASE_PATH, "summaries")
OUT_DIR = os.path.join(BASE_PATH, "plots")
os.makedirs(OUT_DIR, exist_ok=True)

# =============================================================================
# 2. Collect all allele partitions
# =============================================================================
partitions = [d for d in os.listdir(DATASET_PATH) if d.startswith("allele_safe=")]
if not partitions:
    raise FileNotFoundError(f"No allele partitions found in {DATASET_PATH}")

print(f"[INFO] Found {len(partitions)} allele partitions.")

summary_records = []

for part in partitions:
    allele = part.split("=", 1)[1]
    part_path = os.path.join(DATASET_PATH, part)

    # Read all part files for this allele
    df = pq.ParquetDataset(part_path).read_pandas().to_pandas()

    # Compute counts
    total = len(df)
    if total == 0:
        continue

    human_dom = (df["Binding_Class"] == "Human_dominant").sum()
    viral_dom = (df["Binding_Class"] == "Viral_dominant").sum()
    dual = (df["Binding_Class"] == "Equivalent_binding").sum()
    non_bind = (df["Binding_Class"] == "Non_binder").sum()
    allele = df["allele"].iloc[10]  

    human_pct = 100 * human_dom / total
    viral_pct = 100 * viral_dom / total

    # Extract HLA class (A/B/C)
    match = re.search(r"HLA[-_: ]?([A-G])", allele, flags=re.IGNORECASE)
    hla_class = match.group(1).upper() if match else "Unknown"

    summary_records.append({
        "Allele": allele,
        "HLA_type": f"HLA-{hla_class}",
        "Total": total,
        "Human_dominant_pct": human_pct,
        "Viral_dominant_pct": viral_pct,
        "Dual_candidates": dual,
        "Non_binders": non_bind
    })

# =============================================================================
# 3. Create summary DataFrame
# =============================================================================
allele_summary = pd.DataFrame(summary_records)
if allele_summary.empty:
    raise ValueError("No allele summary could be generated — dataset empty or inconsistent.")

allele_types = sorted(allele_summary["HLA_type"].unique())
print(f"[INFO] Included allele types: {', '.join(allele_types)}")


# =============================================================================
# 4. Quantile thresholds
# =============================================================================
q_hi_human = allele_summary["Human_dominant_pct"].quantile(0.95)
q_hi_viral = allele_summary["Viral_dominant_pct"].quantile(0.95)
print(f"[INFO] Quantile thresholds → Human ≥ {q_hi_human:.2f}%, Viral ≥ {q_hi_viral:.2f}%")


# =============================================================================
# 5. Classification & plotting setup
# =============================================================================
def classify_bias(row):
    if row["Human_dominant_pct"] >= q_hi_human and row["Viral_dominant_pct"] >= q_hi_viral:
        return "Dual_high"
    elif row["Human_dominant_pct"] >= q_hi_human:
        return "Human_dominant"
    elif row["Viral_dominant_pct"] >= q_hi_viral:
        return "Viral_dominant"
    else:
        return "Neutral"

allele_summary["Dominance_category"] = allele_summary.apply(classify_bias, axis=1)

allele_summary["Dominance_category"].unique()

color_map = {
    "Dual_high": "#8e44ad",      # purple
    "Human_dominant": "#27ae60", # green
    "Viral_dominant": "#95a5a6", # red
    "Neutral": "#95a5a6"         # gray
}
allele_summary["color"] = allele_summary["Dominance_category"].map(color_map)

allele_summary["log_pairs"] = np.log10(allele_summary["Total"] + 1)
allele_summary["size"] = 40 + 60 * (allele_summary["log_pairs"] / allele_summary["log_pairs"].max())

# =============================================================================
# 6. Save annotated summary
# =============================================================================
allele_summary.to_csv(os.path.join(OUT_DATA_DIR, "allele_summary_from_parquet.csv"), index=False)
print(f"[INFO] Summary saved → allele_summary_from_parquet.csv")

# =============================================================================
# 7. Volcano plot
# =============================================================================
plt.figure(figsize=(15, 10))
sns.set(style="whitegrid")

plt.scatter(
    allele_summary["Human_dominant_pct"],
    allele_summary["Viral_dominant_pct"],
    s=allele_summary["size"],
    c=allele_summary["color"],
    alpha=0.6,
    edgecolors="black",
    linewidth=0.3
)


# Threshold lines
plt.axvline(q_hi_human, color="#27ae60", linestyle="--", lw=1.2, label=f"Human ≥ {q_hi_human:.1f}%")
plt.axhline(q_hi_viral, color="#c0392b", linestyle="--", lw=1.2, label=f"Viral ≥ {q_hi_viral:.1f}%")

# Annotate top alleles
dual_high_df = allele_summary[allele_summary["Dominance_category"] == "Dual_high"].copy()
human_dom_df = allele_summary[allele_summary["Dominance_category"] == "Human_dominant"].copy()

# Sort and select top 20 in each
top_dual = (
    dual_high_df.sort_values("Viral_dominant_pct", ascending=False)
    .head(40)
)
top_human = (
    human_dom_df.sort_values("Human_dominant_pct", ascending=False)
    .head(20)
)

# Combine labeled set (no duplicates)
top_labels = pd.concat([top_dual, top_human]).drop_duplicates("Allele")
label_colors = {
    "Human_dominant": "#138a36",  # bolder green
    "Dual_high": "#6f1d99"        # bolder purple
}
texts = []
for _, row in top_labels.iterrows():
    cat = row["Dominance_category"]
    label_color = label_colors.get(cat, "#333333")  # fallback gray for others
    texts.append(
        plt.text(
            row["Human_dominant_pct"],
            row["Viral_dominant_pct"],
            row["Allele"],
            fontsize=12,
            fontweight="bold",
            color=label_color
        )
    )


plt.scatter(
    top_labels["Human_dominant_pct"],
    top_labels["Viral_dominant_pct"],
    s=top_labels["size"] * 1.3,
    edgecolors="black",
    linewidth=0.6,
    facecolor=top_labels["color"],
    alpha=0.9,
    zorder=3
)


adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
    expand_points=(1.5, 1.5),
    expand_text=(1.2, 1.2),
    force_points=(0.5, 0.5),
    force_text=(0.5, 0.5)
)

# Title and labels
if len(allele_types) == 1:
    type_label = f"{allele_types[0]}"
elif len(allele_types) == 2:
    type_label = f"{allele_types[0]} & {allele_types[1]}"
else:
    type_label = ", ".join(allele_types)

plt.xlabel("% Human-dominant peptide pairs", fontsize=11)
plt.ylabel("% Viral-dominant peptide pairs", fontsize=11)
plt.title(f"Mimicry Volcano Plot ({type_label})", fontsize=14, weight="bold")
plt.legend(frameon=False, fontsize=8, loc="lower right")
plt.grid(alpha=0.25)
plt.tight_layout()

# Save figure
allele_str = "_".join(a.split("-")[1] for a in allele_types)
plot_filename = f"volcano_HLA_{allele_str}.png"
plot_path = os.path.join(OUT_DIR, plot_filename)

plt.savefig(plot_path, dpi=300)
plt.show()

print(f"✅ Volcano plot saved → {plot_path}")


