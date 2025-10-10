"""
================================================================================
Molecular Mimicry Risk Classification
================================================================================
Goal:
    Identify viral peptides that mimic self (human) peptides based on MHC
    binding affinity and sequence similarity, representing potential triggers
    for autoimmune responses.

Binding Behavior and Interpretation:
    - Mimicry_candidate:
        Both viral and human peptides are strong MHC binders (Rank ≤ 2)
        and have similar binding affinities (|ΔBA_score| ≤ 0.5).
        → High autoimmune risk due to cross-reactive T-cell recognition.

    - Viral_dominant:
        Both peptides bind well, but the viral peptide binds significantly
        stronger (ΔBA_score ≤ -0.5).
        → Possible infection-triggered activation that may later cross-react
          with self; moderate autoimmune risk.

    - Human_dominant:
        Both peptides bind well, but the human peptide binds significantly
        stronger (ΔBA_score ≥ +0.5).
        → Tolerogenic; likely central tolerance to self epitope; low risk.

    - Non_binder:
        Either peptide fails to meet binding threshold (Rank > 2).
        → Immunologically silent; not expected to contribute to mimicry.

Quantitative thresholds:
    • Binding threshold: BA_Rank ≤ 2  (NetMHCpan standard for binders)
    • Affinity similarity: |ΔBA_score| ≤ 0.5  (log-scale window)
================================================================================
"""

# ====================================================
# Input/output paths and packages
# ====================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

out_dir = "/ix/djishnu/Priyamvada/virauto/results/virscan/mimicry_analysis/9_mers/type1/"
hla_type = "A"
netmhc_pan_file = f"/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_processed/type1_{hla_type}_all_predictions_processed.txt.gz"


os.makedirs(os.path.join(out_dir, "data"), exist_ok=True)
os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)

# ===================================================
# Annotate and Classify Viral-Human Peptide Pairs 
# ===================================================
merged_diff = pd.read_csv(netmhc_pan_file, sep="\t", compression = 'gzip')

# --- Define thresholds ---
BA_rank_thresh = 2.0
BA_score_delta = 0.5

# --- Compute ΔBA_score and binding classification ---
merged_diff["BA_score_diff"] = merged_diff["BA_score_viral"] - merged_diff["BA_score_human"]

both_bind = (merged_diff["BA_Rank_viral"] <= BA_rank_thresh) & (
    merged_diff["BA_Rank_human"] <= BA_rank_thresh
)

merged_diff["Mimicry_class"] = "Non_binder"
merged_diff.loc[both_bind & (np.abs(merged_diff["BA_score_diff"]) <= BA_score_delta), "Mimicry_class"] = "Equivalent_binding"
merged_diff.loc[both_bind & (merged_diff["BA_score_diff"] < -BA_score_delta), "Mimicry_class"] = "Viral_dominant"
merged_diff.loc[both_bind & (merged_diff["BA_score_diff"] > BA_score_delta), "Mimicry_class"] = "Human_dominant"

# ===================================================
# Summarize mimicry per HLA allele
# ===================================================
allele_summary = (
    merged_diff.groupby(["Allele", "Mimicry_class"])
    .size()
    .reset_index(name="count")
    .pivot(index="Allele", columns="Mimicry_class", values="count")
    .fillna(0)
)

allele_summary["Total"] = allele_summary.sum(axis=1)
for col in ["Equivalent_binding", "Viral_dominant", "Human_dominant", "Non_binder"]:
    if col in allele_summary.columns:
        allele_summary[f"{col}_pct"] = 100 * allele_summary[col] / allele_summary["Total"]

allele_summary = allele_summary.sort_values("Equivalent_binding_pct", ascending=False)
allele_summary.to_csv(os.path.join(out_dir, f"data/allele_mimicry_summary_type1_{hla_type}.csv"))
print("✅ Saved allele-wise mimicry summary")

# ===================================================
# Barplot: per-allele mimicry class distribution
# ===================================================
plot_df = allele_summary.reset_index().melt(
    id_vars="Allele",
    value_vars=[col for col in allele_summary.columns if col.endswith("_pct")],
    var_name="Class",
    value_name="Percent"
)

plt.figure(figsize=(12, 6))
for cls, color in zip(
    ["Equivalent_binding_pct", "Viral_dominant_pct", "Human_dominant_pct"],
    ["#777978", "#339ade", "#c0392b"],
):
    subset = plot_df[plot_df["Class"] == cls]
    plt.bar(subset["Allele"], subset["Percent"], label=cls.replace("_pct", ""), color=color)

plt.xticks(rotation=90)
plt.ylabel("% of peptide pairs")
plt.title("Per-allele distribution of mimicry patterns")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, f"plots/allele_mimicry_summary_type1_{hla_type}.png"), dpi=300)
plt.close()
print("✅ Saved allele-wise mimicry summary plot")

# ===================================================
# Volcano-style plot: allele bias in viral vs self binding
# ===================================================
allele_stats = (
    merged_diff.groupby("Allele")
    .agg(
        mean_diff=("BA_score_diff", "mean"),
        median_diff=("BA_score_diff", "median"),
        n_pairs=("BA_score_diff", "size"),
    )
    .reset_index()
)

plt.figure(figsize=(8, 6))
plt.scatter(
    allele_stats["mean_diff"],
    np.log10(allele_stats["n_pairs"]),
    c=np.where(allele_stats["mean_diff"] > 0, "#27ae60", "#c0392b"),
    alpha=0.7,
    s=60,
)
plt.axvline(0, color="gray", linestyle="--", lw=1)
plt.axvline(BA_score_delta, color="#27ae60", linestyle="--", lw=1, label="Human dominant (> +0.5)")
plt.axvline(-BA_score_delta, color="#c0392b", linestyle="--", lw=1, label="Viral dominant (< -0.5)")
plt.xlabel("Mean ΔBA_score (Viral − Human)")
plt.ylabel("log₁₀(Number of peptide pairs)")
plt.title(f"Allele bias in viral vs self peptide binding (HLA-{hla_type})")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, f"plots/allele_binding_bias_type1_{hla_type}.png"), dpi=300)
plt.close()
print("✅ Saved allele binding bias (volcano-style) plot")
