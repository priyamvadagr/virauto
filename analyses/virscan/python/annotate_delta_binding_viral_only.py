"""
================================================================================
Molecular Mimicry Risk Classification
================================================================================
Goal:
    Identify viral peptides that mimic self (human) peptides based on MHC
    binding affinity and sequence similarity, representing potential triggers
    for autoimmune responses. Pick pairs where only viral peptide is a strong binder.

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
# Config 
# ====================================================
k_mer = 9
type_ = 'type1'
class_ = 'A'
base_out_dir = f"/ix/djishnu/Priyamvada/virauto/results/virscan/mimicry_analysis/viral_strong_binders/"
out_dir = f"{k_mer}_mers/{type_}/data"
netmhc_pan_file = f"/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/{k_mer}_mers/{type_}/{type_}_processed/{type_}_{class_}_all_predictions_processed.txt"

# ====================================================
# Packages
# ====================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ===================================================
# Annotate and Classify Viral-Human Peptide Pairs 
# ===================================================
merged_diff = pd.read_csv(netmhc_pan_file, sep="\t")

# --- Define thresholds ---
BA_rank_thresh = 2.0
BA_score_delta = 0.5

# --- Compute ΔBA_score and binding classification ---
merged_diff["BA_score_diff"] = merged_diff["BA_score_viral"] - merged_diff["BA_score_human"]

both_bind = (merged_diff["BA_Rank_viral"] <= BA_rank_thresh)

merged_diff["Mimicry_class"] = "Non_binder"
merged_diff.loc[both_bind & (np.abs(merged_diff["BA_score_diff"]) <= BA_score_delta), "Mimicry_class"] = "Equivalent_binding"
merged_diff.loc[both_bind & (merged_diff["BA_score_diff"] < -BA_score_delta), "Mimicry_class"] = "Viral_dominant"
merged_diff.loc[both_bind & (merged_diff["BA_score_diff"] > BA_score_delta), "Mimicry_class"] = "Human_dominant"

# ====================================================
# Save annotated file 
# ====================================================
os.makedirs(os.path.join(out_dir, "data"), exist_ok=True)
merged_diff.to_csv(os.path.join(out_dir, f'{type_}_{class_}_mimicry_class_annotate.csv'))
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
allele_summary.to_csv(os.path.join(out_dir, f"allele_mimicry_summary_{type_}_{class_}.csv"))
print("✅ Saved allele-wise mimicry summary")
