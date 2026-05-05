#!/usr/bin/env python3
"""
======================================================================
Script: plot_mimicry_summary.py
Description:
    Generate summary plots for strong molecular mimicry candidates.
    Enrichment plots are adjusted for IEDB ascertainment bias.
    Statistical tests: binomial (HLA/organism), Fisher's exact (HLA×organism).

Dependencies:
    pandas, matplotlib, seaborn, scipy, statsmodels

Usage:
    python plot_mimicry_summary.py \
        --mhc-class I \
        --input /path/to/strong_mimicry.csv.gz \
        --blast-input /path/to/blast_pairs.csv.gz \
        --pair-id-map /path/to/pair_id_mapping.csv.gz \
        --hla-risk /path/to/HLA_autoimmunity_classification_standard.tsv \
        --iedb-fasta /path/to/iedb_epitopes.fasta \
        --fig-dir /path/to/output_figures
======================================================================
"""

import argparse
import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from scipy.stats import binom, fisher_exact
from statsmodels.stats.multitest import multipletests
from collections import Counter

# ====================================================
# CLI
# ====================================================
parser = argparse.ArgumentParser(description="Plot mimicry summary figures")
parser.add_argument("--mhc-class", choices=["I", "II"], required=True,
                    help="MHC class (I or II)")
parser.add_argument("--input", required=True,
                    help="Path to strong mimicry CSV")
parser.add_argument("--blast-input", required=True,
                    help="Path to BLAST filtered pairs CSV")
parser.add_argument("--pair-id-map", required=True,
                    help="Path to pair ID mapping CSV")
parser.add_argument("--hla-risk", required=True,
                    help="Path to HLA autoimmunity risk annotation TSV (standard nomenclature)")
parser.add_argument("--iedb-fasta", required=True,
                    help="Path to IEDB epitope FASTA (pre-BLAST)")
parser.add_argument("--fig-dir", required=True,
                    help="Output directory for figures")
args = parser.parse_args()

MHC_CLASS = args.mhc_class
INPUT_FILE = args.input
BLAST_INPUT_FILE = args.blast_input
PAIR_ID_MAP_FILE = args.pair_id_map
HLA_RISK_FILE = args.hla_risk
IEDB_FASTA_FILE = args.iedb_fasta
FIG_DIR = args.fig_dir
RISK_CLASS_KEY = "ClassI" if MHC_CLASS == "I" else "ClassII"

os.makedirs(FIG_DIR, exist_ok=True)
print(f"=== MHC Class {MHC_CLASS} ===\n")

# Plot style
sns.set_style("whitegrid")
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
})

# ====================================================
# Helpers
# ====================================================

def binomial_enrichment_test(k_values, n_values, global_rate):
    pvals = []
    for k, n in zip(k_values, n_values):
        if n == 0:
            pvals.append(1.0)
        else:
            pval = 1.0 - binom.cdf(k - 1, n, global_rate)
            pvals.append(min(pval, 1.0))
    pvals = np.array(pvals)
    if len(pvals) > 0 and (pvals < 1.0).any():
        _, qvals, _, _ = multipletests(pvals, method="fdr_bh")
    else:
        qvals = pvals.copy()
    return pvals, qvals


def qval_to_stars(q):
    if q < 0.001:
        return "***"
    elif q < 0.01:
        return "**"
    elif q < 0.05:
        return "*"
    return ""


def netmhciipan_to_standard(allele):
    """Convert NetMHCIIpan allele format to standard nomenclature.

    DRB1_0101              -> HLA-DRB1*01:01
    DRB3_0301              -> HLA-DRB3*03:01
    HLA-DQA10102-DQB10602  -> HLA-DQA1*01:02/DQB1*06:02
    HLA-DPA10103-DPB10401  -> HLA-DPA1*01:03/DPB1*04:01
    Already standard       -> returned as-is
    """
    # already standard format
    if "*" in allele and ":" in allele:
        return allele

    # paired DQ/DP: HLA-DQA10102-DQB10602
    m = re.match(r"HLA-(D[QP]A1)(\d{4})-(D[QP]B1)(\d{4})$", allele)
    if m:
        a_gene, a_digits, b_gene, b_digits = m.groups()
        return (f"HLA-{a_gene}*{a_digits[:2]}:{a_digits[2:]}/"
                f"{b_gene}*{b_digits[:2]}:{b_digits[2:]}")

    # single chain: DRB1_0101, DRB3_0301, DRB4_0101, DRB5_0101
    m = re.match(r"(DRB[1345]|DPA1|DPB1|DQA1|DQB1)_(\d{4})$", allele)
    if m:
        gene, digits = m.groups()
        return f"HLA-{gene}*{digits[:2]}:{digits[2:]}"

    return allele


def load_hla_risk(filepath, class_key):
    """Load HLA risk annotations, filtered to the specified MHC class.

    Expects the standardized TSV (tab-delimited, HLA column already in
    standard nomenclature e.g. HLA-DRB1*03:01).
    """
    risk_df = pd.read_csv(filepath, sep="\t")

    # filter to requested class
    risk_cls = risk_df[risk_df["HLA_Class"] == class_key].copy()

    predisposing = set(risk_cls.loc[risk_cls["Association"] == "Predisposing", "HLA"])
    protective = set(risk_cls.loc[risk_cls["Association"] == "Protective", "HLA"])

    predisposing_diseases = {}
    for _, row in risk_cls[risk_cls["Association"] == "Predisposing"].iterrows():
        predisposing_diseases.setdefault(row["HLA"], []).append(row["Disease"])

    return risk_cls, predisposing, protective, predisposing_diseases


def classify_allele_risk(allele, predisposing, protective):
    """Classify an allele as predisposing/protective/not annotated.

    Converts NetMHCIIpan format to standard before matching.
    For Class II paired alleles, checks if any annotated chain
    appears as a substring.
    """
    std = netmhciipan_to_standard(allele)

    if std in predisposing:
        return "Predisposing"
    if std in protective:
        return "Protective"
    # substring match for paired alleles
    for p in predisposing:
        short = p.replace("HLA-", "")
        if short in std:
            return "Predisposing"
    for p in protective:
        short = p.replace("HLA-", "")
        if short in std:
            return "Protective"
    return "Not annotated"


def get_allele_diseases(allele, predisposing_diseases):
    """Get diseases for an allele, including substring matches for pairs."""
    std = netmhciipan_to_standard(allele)
    diseases = predisposing_diseases.get(std, [])
    if not diseases:
        for p, d_list in predisposing_diseases.items():
            short = p.replace("HLA-", "")
            if short in std:
                diseases.extend(d_list)
    return diseases


# ====================================================
# Load data
# ====================================================
print("Loading strong mimicry candidates...")
df = pd.read_csv(INPUT_FILE, low_memory=False)
print(f"  Total records: {len(df):,}")
print(f"  Unique pair_ids: {df['pair_id'].nunique():,}")
print(f"  Unique viral peptides: {df['viral_peptide'].nunique():,}")
print(f"  Unique human peptides: {df['human_peptide'].nunique():,}")

df["source_organism_short"] = (
    df["source_organism"]
    .str.replace(r"\s*\(.*?\)", "", regex=True)
    .str.strip()
)
df["hla_locus"] = df["mhc_allele"].str.extract(r"(HLA-[ABC])")

# ====================================================
# Load IEDB epitope FASTA (organism denominator)
# ====================================================
print("\nLoading IEDB epitope FASTA (pre-BLAST)...")
fasta_org_per_sid = {}
with open(IEDB_FASTA_FILE) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        parts = line[1:].strip().split("|")
        if len(parts) >= 3:
            sid = parts[0]
            org = parts[2].replace("_", " ")
            fasta_org_per_sid[sid] = org

print(f"  Total unique structure_ids: {len(fasta_org_per_sid):,}")

org_sid_counts = Counter()
for sid, org in fasta_org_per_sid.items():
    org_sid_counts[org] += 1
epitopes_per_org_fasta = pd.Series(org_sid_counts, name="fasta_epitopes")
print(f"  Unique organisms: {len(epitopes_per_org_fasta)}")

org_short_map_fasta = {
    org: pd.Series([org]).str.replace(r"\s*\(.*?\)", "", regex=True).str.strip().iloc[0]
    for org in epitopes_per_org_fasta.index
}

# ====================================================
# Load BLAST input + pair IDs (HLA denominator)
# ====================================================
print("\nLoading BLAST input file + pair IDs...")
blast_df = pd.read_csv(BLAST_INPUT_FILE, low_memory=False)
pair_map = pd.read_csv(PAIR_ID_MAP_FILE, low_memory=False)

blast_df["pair_key"] = (
    blast_df["structure_id"].astype(str) + "_" +
    blast_df["hu_prot_id"].astype(str) + "_" +
    blast_df["viral_sequence"].astype(str) + "_" +
    blast_df["human_mimic_sequence"].astype(str)
)
key_to_pair_id = pair_map.drop_duplicates("pair_key").set_index("pair_key")["pair_id"]
blast_df["pair_id"] = blast_df["pair_key"].map(key_to_pair_id)
blast_df = blast_df.dropna(subset=["pair_id"])
print(f"  BLAST pairs with pair_id: {len(blast_df):,}")

pairs_per_hla_input = (
    blast_df.groupby("mhc_allele")["pair_id"].nunique().rename("input_pairs")
)
print(f"  HLA alleles in input: {len(pairs_per_hla_input)}")

# ====================================================
# Load HLA-autoimmunity risk annotations
# ====================================================
print("\nLoading HLA risk annotations...")
hla_risk_df, predisposing_alleles, protective_alleles, predisposing_diseases = load_hla_risk(
    HLA_RISK_FILE, RISK_CLASS_KEY
)
print(f"  Risk alleles loaded ({RISK_CLASS_KEY}): "
      f"{len(predisposing_alleles)} predisposing, {len(protective_alleles)} protective")

df["hla_risk"] = df["mhc_allele"].apply(
    lambda x: classify_allele_risk(x, predisposing_alleles, protective_alleles)
)
print(f"\n  Mimicry pairs by HLA risk status:")
for status, count in df["hla_risk"].value_counts().items():
    print(f"    {status}: {count:,}")

# Build expanded sets that include matched paired alleles for coloring
all_pred_matched = set()
all_prot_matched = set()
for allele in df["mhc_allele"].unique():
    risk = classify_allele_risk(allele, predisposing_alleles, protective_alleles)
    if risk == "Predisposing":
        all_pred_matched.add(allele)
    elif risk == "Protective":
        all_prot_matched.add(allele)

disease_alleles = all_pred_matched | all_prot_matched

# ====================================================
# Figure 1: HLA allele conversion rate
# ====================================================
print("\nPlotting HLA allele conversion rate (Figure 1)...")

strong_per_hla = df.groupby("mhc_allele")["pair_id"].nunique().rename("strong_pairs")
hla_conversion = pairs_per_hla_input.to_frame().join(strong_per_hla, how="left").fillna(0)
hla_conversion["strong_pairs"] = hla_conversion["strong_pairs"].astype(int)
hla_conversion["fraction"] = hla_conversion["strong_pairs"] / hla_conversion["input_pairs"]

global_rate_hla = hla_conversion["strong_pairs"].sum() / hla_conversion["input_pairs"].sum()
print(f"  Global HLA conversion rate: {global_rate_hla:.4f}")

hla_pvals, hla_qvals = binomial_enrichment_test(
    hla_conversion["strong_pairs"].values,
    hla_conversion["input_pairs"].values,
    global_rate_hla,
)
hla_conversion["qval"] = hla_qvals
hla_conversion["stars"] = [qval_to_stars(q) for q in hla_qvals]

n_sig_hla = (hla_conversion["qval"] < 0.05).sum()
print(f"  Significant alleles (FDR < 0.05): {n_sig_hla}")

hla_conv_plot = (
    hla_conversion[hla_conversion["strong_pairs"] > 0]
    .sort_values("input_pairs", ascending=False)
    .head(20)
)

fig, ax = plt.subplots(figsize=(11, 7))
bar_colors = [
    "#D32F2F" if a in all_pred_matched else
    "#1976D2" if a in all_prot_matched else "#BDBDBD"
    for a in hla_conv_plot.index
]

ax.barh(range(len(hla_conv_plot)), hla_conv_plot["fraction"].values[::-1],
        color=bar_colors[::-1], edgecolor="white")
ax.set_yticks(range(len(hla_conv_plot)))
ax.set_yticklabels(hla_conv_plot.index[::-1], fontsize=10)

allele_list_rev = list(hla_conv_plot.index[::-1])
for i, label in enumerate(ax.get_yticklabels()):
    a = allele_list_rev[i]
    if a in all_pred_matched:
        label.set_fontweight("bold")
        label.set_color("#D32F2F")
    elif a in all_prot_matched:
        label.set_color("#1976D2")

for i, (allele, row) in enumerate(hla_conv_plot.iloc[::-1].iterrows()):
    stars = row.get("stars", "")
    if stars:
        ax.text(row["fraction"] + 0.003, i, stars,
                va="center", fontsize=14, fontweight="bold", color="#333")
    diseases = get_allele_diseases(allele, predisposing_diseases)
    if diseases:
        disease_str = ", ".join(sorted(set(diseases)))
        ax.text(row["fraction"] + 0.003, i - 0.3,
                disease_str, va="center", fontsize=9,
                color="#D32F2F", style="italic")

ax.set_xlabel("Conversion rate (strong mimicry pairs / input pairs)")
ax.set_title(f"HLA allele conversion rate to strong mimicry (MHC Class {MHC_CLASS})\n"
             "(adjusted for number of candidate pairs per allele)")
legend_elements = [
    Patch(facecolor="#D32F2F", label="Predisposing"),
    Patch(facecolor="#1976D2", label="Protective"),
    Patch(facecolor="#BDBDBD", label="Not annotated"),
]
ax.legend(handles=legend_elements, fontsize=11, loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "01_hla_conversion_rate.png"))
plt.close()
print("  ✅ 01_hla_conversion_rate.png")

# ====================================================
# Figure 1b: Disease-associated HLA alleles only
# ====================================================
print("Plotting disease-associated HLA conversion rate (Figure 1b)...")

hla_conv_disease = (
    hla_conversion[
        (hla_conversion.index.isin(disease_alleles)) &
        (hla_conversion["strong_pairs"] > 0)
    ].sort_values("input_pairs", ascending=False)
)

if not hla_conv_disease.empty:
    fig, ax = plt.subplots(figsize=(11, max(7, len(hla_conv_disease) * 0.4)))
    bar_colors = [
        "#D32F2F" if a in all_pred_matched else "#1976D2"
        for a in hla_conv_disease.index
    ]

    ax.barh(range(len(hla_conv_disease)), hla_conv_disease["fraction"].values[::-1],
            color=bar_colors[::-1], edgecolor="white")
    ax.set_yticks(range(len(hla_conv_disease)))
    ax.set_yticklabels(hla_conv_disease.index[::-1], fontsize=10)

    allele_list_rev = list(hla_conv_disease.index[::-1])
    for i, label in enumerate(ax.get_yticklabels()):
        a = allele_list_rev[i]
        if a in all_pred_matched:
            label.set_fontweight("bold")
            label.set_color("#D32F2F")
        else:
            label.set_color("#1976D2")

    for i, (allele, row) in enumerate(hla_conv_disease.iloc[::-1].iterrows()):
        stars = row.get("stars", "")
        if stars:
            ax.text(row["fraction"] + 0.003, i, stars,
                    va="center", fontsize=14, fontweight="bold", color="#333")
        diseases = get_allele_diseases(allele, predisposing_diseases)
        if diseases:
            disease_str = ", ".join(sorted(set(diseases)))
            ax.text(row["fraction"] + 0.003, i - 0.3,
                    disease_str, va="center", fontsize=9,
                    color="#D32F2F", style="italic")

    ax.set_xlabel("Conversion rate (strong mimicry pairs / input pairs)")
    ax.set_title(f"Disease-associated HLA alleles: conversion rate (MHC Class {MHC_CLASS})\n"
                 "(red = predisposing, blue = protective)")
    legend_elements = [
        Patch(facecolor="#D32F2F", label="Predisposing"),
        Patch(facecolor="#1976D2", label="Protective"),
    ]
    ax.legend(handles=legend_elements, fontsize=11, loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "01b_hla_disease_conversion_rate.png"))
    plt.close()
    print("  ✅ 01b_hla_disease_conversion_rate.png")
else:
    print("  ⚠️  No disease-associated alleles with strong mimicry pairs — skipping 1b")

# ====================================================
# Figure 2: Source organism conversion rate
# ====================================================
print("Plotting source organism conversion rate (Figure 2)...")

strong_epi_per_org = df.groupby("source_organism")["structure_id"].nunique().rename("strong_epitopes")
org_conversion = epitopes_per_org_fasta.to_frame().join(strong_epi_per_org, how="left").fillna(0)
org_conversion["strong_epitopes"] = org_conversion["strong_epitopes"].astype(int)
org_conversion["fraction"] = org_conversion["strong_epitopes"] / org_conversion["fasta_epitopes"]

global_rate_org = org_conversion["strong_epitopes"].sum() / org_conversion["fasta_epitopes"].sum()
print(f"  Global organism conversion rate: {global_rate_org:.4f}")

org_pvals, org_qvals = binomial_enrichment_test(
    org_conversion["strong_epitopes"].values,
    org_conversion["fasta_epitopes"].values,
    global_rate_org,
)
org_conversion["qval"] = org_qvals
org_conversion["stars"] = [qval_to_stars(q) for q in org_qvals]

n_sig_org = (org_conversion["qval"] < 0.05).sum()
print(f"  Significant organisms (FDR < 0.05): {n_sig_org}")

org_conversion["org_short"] = org_conversion.index.map(lambda x: org_short_map_fasta.get(x, x))

org_conv_plot = (
    org_conversion[org_conversion["strong_epitopes"] > 0]
    .sort_values("fasta_epitopes", ascending=False)
    .head(15)
)

fig, ax = plt.subplots(figsize=(11, 7))
ax.barh(range(len(org_conv_plot)), org_conv_plot["fraction"].values[::-1],
        color=sns.color_palette("Set2", len(org_conv_plot))[::-1], edgecolor="white")
ax.set_yticks(range(len(org_conv_plot)))
ax.set_yticklabels(org_conv_plot["org_short"].values[::-1], fontsize=10)

for i, (_, row) in enumerate(org_conv_plot.iloc[::-1].iterrows()):
    stars = row.get("stars", "")
    if stars:
        ax.text(row["fraction"] + 0.003, i, stars,
                va="center", fontsize=14, fontweight="bold", color="#333")

ax.set_xlabel("Conversion rate (strong mimicry epitopes / total IEDB epitopes)")
ax.set_title(f"Source organism conversion rate to strong mimicry (MHC Class {MHC_CLASS})\n"
             "(adjusted for total epitopes per organism in IEDB)")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "02_organism_conversion_rate.png"))
plt.close()
print("  ✅ 02_organism_conversion_rate.png")

# ====================================================
# Figure 3: Binding affinity comparison
# ====================================================
print("Plotting binding affinity comparison...")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
ax.scatter(df["viral_score_BA"], df["human_score_BA"],
           alpha=0.4, s=15, c="#2196F3", edgecolors="none")
lims = [0, max(df["viral_score_BA"].max(), df["human_score_BA"].max()) * 1.05]
ax.plot(lims, lims, "k--", alpha=0.5, linewidth=1, label="Equal binding")
ax.set_xlabel("Viral BA_score")
ax.set_ylabel("Human BA_score")
ax.set_title("Binding affinity: viral vs human")
ax.legend(fontsize=10)

ax = axes[1]
ax.hist(df["delta_BA_score"].dropna(), bins=50, color="#FF9800", edgecolor="white", alpha=0.8)
ax.axvline(0, color="black", linestyle="--", linewidth=1)
ax.axvline(0.5, color="red", linestyle=":", linewidth=1, label="+0.5 threshold")
ax.axvline(-0.5, color="blue", linestyle=":", linewidth=1, label="-0.5 threshold")
ax.set_xlabel("ΔBA_score (viral − human)")
ax.set_ylabel("Count")
ax.set_title("Distribution of ΔBA_score")
ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "03_binding_comparison.png"))
plt.close()
print("  ✅ 03_binding_comparison.png")

# ====================================================
# Figure 4: Peptide length and mismatch distributions
# ====================================================
print("Plotting peptide properties...")

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

ax = axes[0]
viral_lens = df["viral_peptide"].str.len()
human_lens = df["human_peptide"].str.len()
all_lens = pd.concat([viral_lens, human_lens])
bins = range(int(all_lens.min()), int(all_lens.max()) + 2)
ax.hist(viral_lens, bins=bins, alpha=0.6, label="Viral", color="#2196F3", edgecolor="white")
ax.hist(human_lens, bins=bins, alpha=0.6, label="Human", color="#FF9800", edgecolor="white")
ax.set_xlabel("Peptide length (aa)")
ax.set_ylabel("Count")
ax.set_title("Peptide length distribution")
ax.legend(fontsize=10)

ax = axes[1]
if "n_mismatches" in df.columns:
    mm_counts = df["n_mismatches"].value_counts().sort_index()
    ax.bar(mm_counts.index, mm_counts.values, color="#4CAF50", edgecolor="white")
    ax.set_xlabel("Number of mismatches")
    ax.set_ylabel("Count")
    ax.set_title("Mismatch distribution")
    ax.set_xticks(mm_counts.index)

ax = axes[2]
if "pident" in df.columns:
    ax.hist(df["pident"].dropna(), bins=30, color="#9C27B0", edgecolor="white", alpha=0.8)
    ax.set_xlabel("Percent identity")
    ax.set_ylabel("Count")
    ax.set_title("Sequence identity distribution")

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "04_peptide_properties.png"))
plt.close()
print("  ✅ 04_peptide_properties.png")

# ====================================================
# Figure 5: HLA × Organism dot plot
#   Fisher's exact test for interaction
# ====================================================
print("Plotting HLA × organism conversion rate dot plot...")

input_counts = blast_df.groupby(["mhc_allele", "source_organism"])["pair_id"].nunique().rename("input_pairs")
strong_counts = df.groupby(["mhc_allele", "source_organism"])["pair_id"].nunique().rename("strong_pairs")

conversion = input_counts.to_frame().join(strong_counts, how="left").fillna(0)
conversion["strong_pairs"] = conversion["strong_pairs"].astype(int)
conversion["fraction"] = conversion["strong_pairs"] / conversion["input_pairs"]

# Fisher's exact test
strong_per_hla_all = df.groupby("mhc_allele")["pair_id"].nunique()
strong_per_org_all = df.groupby("source_organism")["pair_id"].nunique()
total_strong_all = df["pair_id"].nunique()

cell_pvals = []
for (allele, org), row in conversion.iterrows():
    k = int(row["strong_pairs"])
    if k == 0:
        cell_pvals.append(1.0)
        continue
    a = k
    b = max(strong_per_hla_all.get(allele, 0) - a, 0)
    c = max(strong_per_org_all.get(org, 0) - a, 0)
    d = max(total_strong_all - a - b - c, 0)
    _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
    cell_pvals.append(pval)

conversion["pval_interaction"] = cell_pvals

mask_tested = conversion["strong_pairs"] > 0
conversion["qval_interaction"] = 1.0
if mask_tested.sum() > 0:
    tested_pvals = conversion.loc[mask_tested, "pval_interaction"].values
    _, tested_qvals, _, _ = multipletests(tested_pvals, method="fdr_bh")
    conversion.loc[mask_tested, "qval_interaction"] = tested_qvals

conversion["stars_interaction"] = [qval_to_stars(q) for q in conversion["qval_interaction"]]

n_sig_cells = (conversion["qval_interaction"] < 0.05).sum()
print(f"  Fisher's exact: {n_sig_cells} significant interactions (FDR < 0.05)")

blast_df["source_organism_short"] = (
    blast_df["source_organism"].str.replace(r"\s*\(.*?\)", "", regex=True).str.strip()
)
org_short_map = (
    blast_df[["source_organism", "source_organism_short"]]
    .drop_duplicates("source_organism")
    .set_index("source_organism")["source_organism_short"]
)

top_alleles = df["mhc_allele"].value_counts().head(15).index.tolist()
top_orgs_full = df["source_organism"].value_counts().head(10).index.tolist()
top_orgs_short = [org_short_map.get(o, o) for o in top_orgs_full]

heatmap_frac = pd.DataFrame(np.nan, index=top_orgs_short, columns=top_alleles)
heatmap_input = pd.DataFrame(0, index=top_orgs_short, columns=top_alleles)
heatmap_stars = pd.DataFrame("", index=top_orgs_short, columns=top_alleles, dtype=str)

for org_full, org_short in zip(top_orgs_full, top_orgs_short):
    for allele in top_alleles:
        key = (allele, org_full)
        if key in conversion.index:
            row = conversion.loc[key]
            heatmap_frac.loc[org_short, allele] = row["fraction"]
            heatmap_input.loc[org_short, allele] = int(row["input_pairs"])
            heatmap_stars.loc[org_short, allele] = row["stars_interaction"]

heatmap_frac = heatmap_frac.fillna(0)

fig, ax = plt.subplots(figsize=(14, 7))
vmax = min(heatmap_frac.max().max() * 1.1, 1.0)
if vmax == 0:
    vmax = 0.01
norm = plt.Normalize(vmin=0, vmax=vmax)
cmap = plt.cm.YlOrRd

max_input = max(heatmap_input.max().max(), 1)
min_dot, max_dot = 30, 600

for i, org_short in enumerate(top_orgs_short):
    for j, allele in enumerate(top_alleles):
        frac = heatmap_frac.loc[org_short, allele]
        n_input = heatmap_input.loc[org_short, allele]
        if n_input == 0:
            continue
        dot_size = min_dot + (n_input / max_input) * (max_dot - min_dot)
        ax.scatter(j + 0.5, i + 0.5, s=dot_size, c=[cmap(norm(frac))],
                   edgecolors="white", linewidths=0.5, zorder=3)
        stars = heatmap_stars.loc[org_short, allele]
        if stars:
            ax.text(j + 0.5, i + 0.9, stars, ha="center", va="center",
                    fontsize=14, fontweight="bold", color="#333", zorder=4)

for i in range(len(top_orgs_short) + 1):
    ax.axhline(i, color="#E0E0E0", linewidth=0.5, zorder=1)
for j in range(len(top_alleles) + 1):
    ax.axvline(j, color="#E0E0E0", linewidth=0.5, zorder=1)

ax.set_xlim(0, len(top_alleles))
ax.set_ylim(0, len(top_orgs_short))
ax.set_xticks([j + 0.5 for j in range(len(top_alleles))])
ax.set_xticklabels(top_alleles, rotation=45, ha="right", fontsize=10)
ax.set_yticks([i + 0.5 for i in range(len(top_orgs_short))])
ax.set_yticklabels(top_orgs_short, fontsize=10)
ax.invert_yaxis()

for label in ax.get_xticklabels():
    a = label.get_text()
    if a in all_pred_matched:
        label.set_color("#D32F2F")
        label.set_fontweight("bold")
    elif a in all_prot_matched:
        label.set_color("#1976D2")

ax.set_xlabel("HLA allele")
ax.set_ylabel("Source organism")
ax.set_title(f"HLA × organism: conversion rate and interaction (MHC Class {MHC_CLASS})\n"
             "(color = conversion rate; dot size = input pairs; "
             "stars = Fisher's exact FDR)")

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.02)
cbar.set_label("Conversion rate")

legend_sizes = sorted(set(max(1, s) for s in [
    int(max_input * 0.1), int(max_input * 0.5), int(max_input)
]))
legend_dots = [
    ax.scatter([], [], s=min_dot + (s / max_input) * (max_dot - min_dot),
               c="gray", edgecolors="white", linewidths=0.5, label=f"{s:,} pairs")
    for s in legend_sizes
]
ax.legend(handles=legend_dots, title="Input pairs", loc="upper left",
          bbox_to_anchor=(1.15, 1.0), labelspacing=2.5,
          fontsize=10, title_fontsize=11, frameon=True)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "05_hla_organism_dotplot.png"))
plt.close()
print("  ✅ 05_hla_organism_dotplot.png")

# ====================================================
# Figure 5b: Disease-associated HLA × Organism dot plot
# ====================================================
print("Plotting disease-associated HLA × organism dot plot (Figure 5b)...")

all_disease_in_data = sorted(
    [a for a in hla_conversion.index if a in disease_alleles and hla_conversion.loc[a, "strong_pairs"] > 0],
    key=lambda a: hla_conversion.loc[a, "input_pairs"], reverse=True,
)

if all_disease_in_data:
    disease_frac = pd.DataFrame(np.nan, index=top_orgs_short, columns=all_disease_in_data)
    disease_input = pd.DataFrame(0, index=top_orgs_short, columns=all_disease_in_data)
    disease_stars = pd.DataFrame("", index=top_orgs_short, columns=all_disease_in_data, dtype=str)

    for org_full, org_short in zip(top_orgs_full, top_orgs_short):
        for allele in all_disease_in_data:
            key = (allele, org_full)
            if key in conversion.index:
                row = conversion.loc[key]
                disease_frac.loc[org_short, allele] = row["fraction"]
                disease_input.loc[org_short, allele] = int(row["input_pairs"])
                disease_stars.loc[org_short, allele] = row["stars_interaction"]

    disease_frac = disease_frac.fillna(0)

    fig, ax = plt.subplots(figsize=(max(10, len(all_disease_in_data) * 0.8), 7))
    vmax_d = max(min(disease_frac.max().max() * 1.1, 1.0), 0.01)
    norm_d = plt.Normalize(vmin=0, vmax=vmax_d)
    max_input_d = max(disease_input.max().max(), 1)

    for i, org_short in enumerate(top_orgs_short):
        for j, allele in enumerate(all_disease_in_data):
            frac = disease_frac.loc[org_short, allele]
            n_input = disease_input.loc[org_short, allele]
            if n_input == 0:
                continue
            dot_size = min_dot + (n_input / max_input_d) * (max_dot - min_dot)
            ax.scatter(j + 0.5, i + 0.5, s=dot_size, c=[cmap(norm_d(frac))],
                       edgecolors="white", linewidths=0.5, zorder=3)
            stars = disease_stars.loc[org_short, allele]
            if stars:
                ax.text(j + 0.5, i + 0.9, stars, ha="center", va="center",
                        fontsize=14, fontweight="bold", color="#333", zorder=4)

    for i in range(len(top_orgs_short) + 1):
        ax.axhline(i, color="#E0E0E0", linewidth=0.5, zorder=1)
    for j in range(len(all_disease_in_data) + 1):
        ax.axvline(j, color="#E0E0E0", linewidth=0.5, zorder=1)

    ax.set_xlim(0, len(all_disease_in_data))
    ax.set_ylim(0, len(top_orgs_short))
    ax.set_xticks([j + 0.5 for j in range(len(all_disease_in_data))])
    ax.set_xticklabels(all_disease_in_data, rotation=45, ha="right", fontsize=9)
    ax.set_yticks([i + 0.5 for i in range(len(top_orgs_short))])
    ax.set_yticklabels(top_orgs_short, fontsize=10)
    ax.invert_yaxis()

    for label in ax.get_xticklabels():
        a = label.get_text()
        if a in all_pred_matched:
            label.set_color("#D32F2F")
            label.set_fontweight("bold")
        elif a in all_prot_matched:
            label.set_color("#1976D2")
            label.set_fontweight("bold")

    ax.set_xlabel("HLA allele (disease-associated only)")
    ax.set_ylabel("Source organism")
    ax.set_title(f"Disease-associated HLA × organism (MHC Class {MHC_CLASS})\n"
                 "(stars = Fisher's exact FDR; red = predisposing, blue = protective)")

    sm_d = plt.cm.ScalarMappable(cmap=cmap, norm=norm_d)
    sm_d.set_array([])
    cbar_d = plt.colorbar(sm_d, ax=ax, pad=0.02)
    cbar_d.set_label("Conversion rate")

    legend_sizes_d = sorted(set(max(1, s) for s in [
        int(max_input_d * 0.1), int(max_input_d * 0.5), int(max_input_d)
    ]))
    legend_dots_d = [
        ax.scatter([], [], s=min_dot + (s / max_input_d) * (max_dot - min_dot),
                   c="gray", edgecolors="white", linewidths=0.5, label=f"{s:,} pairs")
        for s in legend_sizes_d
    ]
    ax.legend(handles=legend_dots_d, title="Input pairs", loc="upper left",
              bbox_to_anchor=(1.15, 1.0), labelspacing=2.5,
              fontsize=10, title_fontsize=11, frameon=True)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "05b_disease_hla_organism_dotplot.png"))
    plt.close()
    print("  ✅ 05b_disease_hla_organism_dotplot.png")
else:
    print("  ⚠️  No disease-associated alleles with strong mimicry pairs — skipping 5b")

# ====================================================
# Figure 6: Human protein enrichment (raw counts)
# ====================================================
print("Plotting human protein enrichment...")

fig, ax = plt.subplots(figsize=(10, 6))
hu_prot_epitopes = (
    df.groupby("hu_prot_name")["viral_peptide"].nunique()
    .sort_values(ascending=False).head(20)
)
ax.barh(range(len(hu_prot_epitopes)), hu_prot_epitopes.values,
        color=sns.color_palette("coolwarm", len(hu_prot_epitopes)))
ax.set_yticks(range(len(hu_prot_epitopes)))
ax.set_yticklabels(hu_prot_epitopes.index, fontsize=9)
ax.set_xlabel("Number of unique viral epitopes with mimicry")
ax.set_title(f"Top 20 human proteins targeted by molecular mimicry (MHC Class {MHC_CLASS})")
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "06_human_protein_enrichment.png"))
plt.close()
print("  ✅ 06_human_protein_enrichment.png")

# ====================================================
# Figure 7: Autoimmune risk HLA analysis
# ====================================================
print("Plotting autoimmune risk HLA analysis...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

risk_counts = df["hla_risk"].value_counts()
ax = axes[0]
risk_colors = {"Predisposing": "#D32F2F", "Protective": "#1976D2", "Not annotated": "#BDBDBD"}
ax.bar(range(len(risk_counts)), risk_counts.values,
       color=[risk_colors.get(r, "#999") for r in risk_counts.index], edgecolor="white")
ax.set_xticks(range(len(risk_counts)))
ax.set_xticklabels(risk_counts.index, fontsize=11)
ax.set_ylabel("Number of mimicry pairs")
ax.set_title(f"Mimicry pairs by HLA autoimmune risk status (MHC Class {MHC_CLASS})")
for i, (status, count) in enumerate(risk_counts.items()):
    ax.text(i, count + max(risk_counts) * 0.01, str(count),
            ha="center", fontsize=11, fontweight="bold")

predisposing_pairs = df[df["hla_risk"] == "Predisposing"]
if not predisposing_pairs.empty:
    disease_records = []
    for _, row in predisposing_pairs.iterrows():
        allele = row["mhc_allele"]
        diseases = get_allele_diseases(allele, predisposing_diseases)
        if not diseases:
            diseases = ["Unknown"]
        for disease in diseases:
            disease_records.append({"disease": disease, "pair_id": row["pair_id"]})

    disease_df = pd.DataFrame(disease_records)
    disease_counts = disease_df["disease"].value_counts().head(15)

    ax = axes[1]
    ax.barh(range(len(disease_counts)), disease_counts.values,
            color="#D32F2F", alpha=0.7, edgecolor="white")
    ax.set_yticks(range(len(disease_counts)))
    ax.set_yticklabels(disease_counts.index, fontsize=10)
    ax.set_xlabel("Number of mimicry pairs")
    ax.set_title("Autoimmune diseases linked to predisposing HLA alleles")
    ax.invert_yaxis()
else:
    axes[1].text(0.5, 0.5, "No predisposing alleles matched",
                 ha="center", va="center", transform=axes[1].transAxes, fontsize=12)
    axes[1].set_title("Autoimmune diseases linked to predisposing HLA alleles")

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "07_autoimmune_risk_hla.png"))
plt.close()
print("  ✅ 07_autoimmune_risk_hla.png")

# ====================================================
# Summary statistics
# ====================================================
print(f"\n{'=' * 60}")
print(f"Summary Statistics (MHC Class {MHC_CLASS})")
print(f"{'=' * 60}")
print(f"  Total strong mimicry records: {len(df):,}")
print(f"  Unique pair IDs: {df['pair_id'].nunique():,}")
print(f"  Unique viral peptides: {df['viral_peptide'].nunique():,}")
print(f"  Unique human peptides: {df['human_peptide'].nunique():,}")
print(f"  Unique HLA alleles: {df['mhc_allele'].nunique()}")
print(f"  Unique source organisms: {df['source_organism_short'].nunique()}")
print(f"  Unique human proteins: {df['hu_prot_name'].nunique()}")
print(f"\n  Median ΔBA_score: {df['delta_BA_score'].median():.4f}")
print(f"  Mean ΔBA_score: {df['delta_BA_score'].mean():.4f}")
print(f"\n  Figures saved to: {FIG_DIR}")
print(f"\n✅ Done.")