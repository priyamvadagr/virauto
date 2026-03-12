#!/usr/bin/env python3
"""
======================================================================
Script: plot_mimicry_summary.py
Description:
    Generate summary plots for strong molecular mimicry candidates.

Input:
    iedb_mhci_strong_mimicry.csv.gz

Output:
    Figures in results/mimicry_analysis/iedb/mhc_i/figures/

Notes:
    Figures 9 (domain ORA) and 10 (protein ORA) are now handled by
    standalone scripts:
      - compute_domain_enrichment.py  (positional domain ORA)
      - compute_protein_enrichment.py (binomial protein ORA)

Dependencies:
    pandas, matplotlib, seaborn

Usage:
    python plot_mimicry_summary.py
======================================================================
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import os
import numpy as np

# ====================================================
# Config
# ====================================================
INPUT_FILE = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/iedb/mhc_i/parsed/iedb_mhci_strong_mimicry.csv.gz"
BLAST_INPUT_FILE = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla.csv.gz"
FIG_DIR = "/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/iedb/mhc_i/figures"
os.makedirs(FIG_DIR, exist_ok=True)

# Plot style
sns.set_style("whitegrid")
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
})

# ====================================================
# Load data
# ====================================================
print("Loading strong mimicry candidates...")
df = pd.read_csv(INPUT_FILE, low_memory=False)
print(f"  Total records: {len(df):,}")
print(f"  Unique pair_ids: {df['pair_id'].nunique():,}")
print(f"  Unique viral peptides: {df['viral_peptide'].nunique():,}")
print(f"  Unique human peptides: {df['human_peptide'].nunique():,}")

# Clean up source organism names
df["source_organism_short"] = (
    df["source_organism"]
    .str.replace(r"\s*\(.*?\)", "", regex=True)
    .str.strip()
)

# Extract locus from allele name
df["hla_locus"] = df["mhc_allele"].str.extract(r"(HLA-[ABC])")

# ====================================================
# Load HLA-autoimmunity risk annotations
# ====================================================
HLA_RISK_FILE = "/ix/djishnu/Priyamvada/pepSource/data/Ogishi_biorxiv_2019/HLA_autoimmunity_classification.txt"

def load_hla_risk(filepath):
    risk_df = pd.read_csv(filepath, sep=r"\s+", engine="python")
    # Strip any surrounding quotes from column names and string values
    risk_df.columns = risk_df.columns.str.strip('"')
    for col in risk_df.select_dtypes(include="object").columns:
        risk_df[col] = risk_df[col].str.strip('"')
    print(f"  HLA risk annotations loaded: {len(risk_df)} entries")

    def convert_hla_name(name):
        parts = name.split("_")
        if len(parts) == 3:
            locus = parts[1]
            digits = parts[2]
            if len(digits) == 4:
                return f"HLA-{locus}*{digits[:2]}:{digits[2:]}"
            elif len(digits) == 2:
                return f"HLA-{locus}*{digits}"
        return name

    risk_df["hla_standard"] = risk_df["HLA"].apply(convert_hla_name)

    predisposing = set(
        risk_df[risk_df["Association"] == "Predisposing"]["hla_standard"]
    )
    protective = set(
        risk_df[risk_df["Association"] == "Protective"]["hla_standard"]
    )

    predisposing_diseases = {}
    for _, row in risk_df[risk_df["Association"] == "Predisposing"].iterrows():
        allele = row["hla_standard"]
        disease = row["Disease"]
        if allele not in predisposing_diseases:
            predisposing_diseases[allele] = []
        predisposing_diseases[allele].append(disease)

    class_i_predisposing = set(
        risk_df[
            (risk_df["Association"] == "Predisposing") &
            (risk_df["HLA_Class"] == "ClassI")
        ]["hla_standard"]
    )
    class_i_protective = set(
        risk_df[
            (risk_df["Association"] == "Protective") &
            (risk_df["HLA_Class"] == "ClassI")
        ]["hla_standard"]
    )

    print(f"  Class I predisposing alleles: {len(class_i_predisposing)}")
    print(f"  Class I protective alleles: {len(class_i_protective)}")

    return risk_df, class_i_predisposing, class_i_protective, predisposing_diseases

print("\nLoading HLA risk annotations...")
hla_risk_df, predisposing_alleles, protective_alleles, predisposing_diseases = load_hla_risk(HLA_RISK_FILE)

df["hla_risk"] = df["mhc_allele"].apply(
    lambda x: "Predisposing" if x in predisposing_alleles else
              ("Protective" if x in protective_alleles else "Not annotated")
)
print(f"\n  Mimicry pairs by HLA risk status:")
for status, count in df["hla_risk"].value_counts().items():
    print(f"    {status}: {count:,}")

# ====================================================
# Figure 1: HLA allele enrichment
# ====================================================
print("\nPlotting HLA allele enrichment...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

allele_counts = df["mhc_allele"].value_counts().head(20)
ax = axes[0]
bar_colors = []
for allele in allele_counts.index:
    if allele in predisposing_alleles:
        bar_colors.append("#D32F2F")
    elif allele in protective_alleles:
        bar_colors.append("#1976D2")
    else:
        bar_colors.append("#BDBDBD")

bars = ax.barh(range(len(allele_counts)), allele_counts.values, color=bar_colors)
ax.set_yticks(range(len(allele_counts)))
ax.set_yticklabels(allele_counts.index, fontsize=9)

ytick_labels = ax.get_yticklabels()
for i, allele in enumerate(allele_counts.index):
    if allele in predisposing_alleles:
        ytick_labels[i].set_fontweight("bold")
        ytick_labels[i].set_color("#D32F2F")
        diseases = predisposing_diseases.get(allele, [])
        if diseases:
            disease_str = ", ".join(sorted(set(diseases)))
            ax.annotate(
                disease_str,
                xy=(allele_counts.values[i], i),
                xytext=(5, 0), textcoords="offset points",
                fontsize=7, color="#D32F2F", va="center"
            )
    elif allele in protective_alleles:
        ytick_labels[i].set_color("#1976D2")

ax.set_xlabel("Number of mimicry pairs")
ax.set_title("Top 20 HLA alleles in strong mimicry candidates")
ax.invert_yaxis()

legend_elements = [
    Patch(facecolor="#D32F2F", label="Predisposing"),
    Patch(facecolor="#1976D2", label="Protective"),
    Patch(facecolor="#BDBDBD", label="Not annotated"),
]
ax.legend(handles=legend_elements, fontsize=8, loc="lower right")

locus_counts = df["hla_locus"].value_counts()
ax = axes[1]
colors = {"HLA-A": "#2196F3", "HLA-B": "#FF9800", "HLA-C": "#4CAF50"}
ax.pie(
    locus_counts.values,
    labels=locus_counts.index,
    colors=[colors.get(l, "#999") for l in locus_counts.index],
    autopct=lambda pct: f"{pct:.1f}%\n({int(pct/100*sum(locus_counts.values))})",
    startangle=90,
    textprops={"fontsize": 11}
)
ax.set_title("HLA locus distribution")

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "01_hla_enrichment.png"))
plt.close()
print("  ✅ 01_hla_enrichment.png")

# ====================================================
# Figure 2: Source organism enrichment
# ====================================================
print("Plotting source organism enrichment...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

org_counts = df["source_organism_short"].value_counts().head(15)
ax = axes[0]
ax.barh(
    range(len(org_counts)),
    org_counts.values,
    color=sns.color_palette("Set2", len(org_counts))
)
ax.set_yticks(range(len(org_counts)))
ax.set_yticklabels(org_counts.index, fontsize=9)
ax.set_xlabel("Number of mimicry pairs")
ax.set_title("Top 15 source organisms")
ax.invert_yaxis()

epitopes_per_org = (
    df.groupby("source_organism_short")["viral_peptide"]
    .nunique()
    .sort_values(ascending=False)
    .head(15)
)
ax = axes[1]
ax.barh(
    range(len(epitopes_per_org)),
    epitopes_per_org.values,
    color=sns.color_palette("Set2", len(epitopes_per_org))
)
ax.set_yticks(range(len(epitopes_per_org)))
ax.set_yticklabels(epitopes_per_org.index, fontsize=9)
ax.set_xlabel("Number of unique viral epitopes")
ax.set_title("Unique viral epitopes per organism")
ax.invert_yaxis()

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "02_source_organism_enrichment.png"))
plt.close()
print("  ✅ 02_source_organism_enrichment.png")

# ====================================================
# Figure 3: Binding affinity comparison (viral vs human)
# ====================================================
print("Plotting binding affinity comparison...")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

ax = axes[0]
ax.scatter(
    df["viral_score_BA"], df["human_score_BA"],
    alpha=0.4, s=15, c="#2196F3", edgecolors="none"
)
lims = [0, max(df["viral_score_BA"].max(), df["human_score_BA"].max()) * 1.05]
ax.plot(lims, lims, "k--", alpha=0.5, linewidth=1, label="Equal binding")
ax.set_xlabel("Viral BA_score")
ax.set_ylabel("Human BA_score")
ax.set_title("Binding affinity: viral vs human")
ax.legend(fontsize=9)

ax = axes[1]
ax.hist(
    df["delta_BA_score"].dropna(), bins=50,
    color="#FF9800", edgecolor="white", alpha=0.8
)
ax.axvline(0, color="black", linestyle="--", linewidth=1)
ax.axvline(0.5, color="red", linestyle=":", linewidth=1, label="+0.5 threshold")
ax.axvline(-0.5, color="blue", linestyle=":", linewidth=1, label="-0.5 threshold")
ax.set_xlabel("ΔBA_score (viral − human)")
ax.set_ylabel("Count")
ax.set_title("Distribution of ΔBA_score")
ax.legend(fontsize=9)

ax = axes[2]
ax.scatter(
    df["viral_rank_BA"], df["human_rank_BA"],
    alpha=0.4, s=15, c="#4CAF50", edgecolors="none"
)
ax.axhline(2, color="red", linestyle=":", alpha=0.5, label="Binder threshold (Rank=2)")
ax.axvline(2, color="red", linestyle=":", alpha=0.5)
ax.set_xlabel("Viral BA_Rank")
ax.set_ylabel("Human BA_Rank")
ax.set_title("Binding rank: viral vs human")
ax.legend(fontsize=9)

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
bins = range(7, 16)
ax.hist(viral_lens, bins=bins, alpha=0.6, label="Viral", color="#2196F3", edgecolor="white")
ax.hist(human_lens, bins=bins, alpha=0.6, label="Human", color="#FF9800", edgecolor="white")
ax.set_xlabel("Peptide length (aa)")
ax.set_ylabel("Count")
ax.set_title("Peptide length distribution")
ax.legend()

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
# Figure 5: HLA × Organism heatmap (fraction of input
#           epitopes that become strong mimicry candidates)
# ====================================================
print("Plotting HLA × organism conversion rate heatmap...")

if os.path.exists(BLAST_INPUT_FILE):
    input_df = pd.read_csv(BLAST_INPUT_FILE, low_memory=False)
    input_df["source_organism_short"] = (
        input_df["source_organism"]
        .str.replace(r"\s*\(.*?\)", "", regex=True)
        .str.strip()
    )

    # Count unique viral epitopes per (HLA, organism) in input
    input_epitopes = (
        input_df
        .groupby(["mhc_allele", "source_organism_short"])["viral_sequence"]
        .nunique()
        .rename("input_epitopes")
    )

    # Count unique viral epitopes per (HLA, organism) in strong mimicry
    strong_epitopes = (
        df
        .groupby(["mhc_allele", "source_organism_short"])["viral_peptide"]
        .nunique()
        .rename("strong_epitopes")
    )

    # Merge and compute fraction
    conversion = pd.concat([input_epitopes, strong_epitopes], axis=1).fillna(0)
    conversion["fraction"] = conversion["strong_epitopes"] / conversion["input_epitopes"]
    conversion = conversion.replace([np.inf, np.nan], 0)

    # Reshape to matrix for top alleles × organisms
    top_alleles = df["mhc_allele"].value_counts().head(15).index
    top_orgs = df["source_organism_short"].value_counts().head(10).index

    # Build heatmap matrix
    heatmap_frac = pd.DataFrame(0.0, index=top_orgs, columns=top_alleles)
    heatmap_count = pd.DataFrame("", index=top_orgs, columns=top_alleles, dtype=str)

    for org in top_orgs:
        for allele in top_alleles:
            key = (allele, org)
            if key in conversion.index:
                row = conversion.loc[key]
                frac = row["fraction"]
                n_strong = int(row["strong_epitopes"])
                n_input = int(row["input_epitopes"])
                heatmap_frac.loc[org, allele] = frac
                heatmap_count.loc[org, allele] = f"{n_strong}/{n_input}"

    fig, ax = plt.subplots(figsize=(14, 7))

    sns.heatmap(
        heatmap_frac,
        annot=heatmap_count,
        fmt="",
        cmap="YlOrRd",
        linewidths=0.5, linecolor="white",
        ax=ax,
        vmin=0, vmax=min(heatmap_frac.max().max() * 1.1, 1.0),
        cbar_kws={"label": "Fraction of input epitopes → strong mimicry"},
    )
    ax.set_xlabel("HLA allele")
    ax.set_ylabel("Source organism")
    ax.set_title(
        "Conversion rate: unique viral epitopes → strong mimicry candidates\n"
        "(annotations: strong / input epitopes per cell)"
    )

    # Color x-tick labels by risk status
    xtick_labels = ax.get_xticklabels()
    for label in xtick_labels:
        allele = label.get_text()
        if allele in predisposing_alleles:
            label.set_color("#D32F2F")
            label.set_fontweight("bold")
        elif allele in protective_alleles:
            label.set_color("#1976D2")

    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "05_hla_organism_conversion_heatmap.png"))
    plt.close()
    print("  ✅ 05_hla_organism_conversion_heatmap.png")
else:
    print("  ⚠️ BLAST input file not found, skipping heatmap")

# ====================================================
# Figure 6: Human protein enrichment
# ====================================================
print("Plotting human protein enrichment...")

fig, ax = plt.subplots(figsize=(10, 6))

hu_prot_epitopes = (
    df.groupby("hu_prot_name")["viral_peptide"]
    .nunique()
    .sort_values(ascending=False)
    .head(20)
)

ax.barh(
    range(len(hu_prot_epitopes)),
    hu_prot_epitopes.values,
    color=sns.color_palette("coolwarm", len(hu_prot_epitopes))
)
ax.set_yticks(range(len(hu_prot_epitopes)))
ax.set_yticklabels(hu_prot_epitopes.index, fontsize=8)
ax.set_xlabel("Number of unique viral epitopes with mimicry")
ax.set_title("Top 20 human proteins targeted by molecular mimicry")
ax.invert_yaxis()

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "06_human_protein_enrichment.png"))
plt.close()
print("  ✅ 06_human_protein_enrichment.png")

# ====================================================
# Figure 7: Enrichment — fold change over input
# ====================================================
print("Plotting enrichment analysis...")

if os.path.exists(BLAST_INPUT_FILE):
    # input_df already loaded above for Figure 5

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 7a: HLA allele fold enrichment
    input_allele_frac = input_df["mhc_allele"].value_counts(normalize=True)
    mimicry_allele_frac = df["mhc_allele"].value_counts(normalize=True)

    common_alleles = mimicry_allele_frac.head(15).index
    fold_change = []
    for a in common_alleles:
        input_frac = input_allele_frac.get(a, 0)
        mimicry_frac = mimicry_allele_frac.get(a, 0)
        if input_frac > 0:
            fold_change.append(mimicry_frac / input_frac)
        else:
            fold_change.append(np.nan)

    fc_series = pd.Series(fold_change, index=common_alleles).dropna().sort_values(ascending=True)

    ax = axes[0]
    bar_colors = []
    for allele in fc_series.index:
        if allele in predisposing_alleles:
            bar_colors.append("#D32F2F")
        elif allele in protective_alleles:
            bar_colors.append("#1976D2")
        else:
            bar_colors.append("#BDBDBD" if fc_series[allele] < 1 else "#616161")

    ax.barh(range(len(fc_series)), fc_series.values, color=bar_colors, edgecolor="white")
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1, alpha=0.7)
    ax.set_yticks(range(len(fc_series)))

    ytick_labels_text = list(fc_series.index)
    ax.set_yticklabels(ytick_labels_text, fontsize=9)
    for i, label in enumerate(ax.get_yticklabels()):
        allele = ytick_labels_text[i]
        if allele in predisposing_alleles:
            label.set_color("#D32F2F")
            label.set_fontweight("bold")
        elif allele in protective_alleles:
            label.set_color("#1976D2")

    ax.set_xlabel("Fold enrichment (mimicry / input)")
    ax.set_title("HLA allele enrichment in strong mimicry")
    legend_elements = [
        Patch(facecolor="#D32F2F", label="Predisposing"),
        Patch(facecolor="#1976D2", label="Protective"),
        Patch(facecolor="#616161", label="Not annotated"),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc="lower right")

    # 7b: Organism fold enrichment
    input_org_frac = input_df["source_organism_short"].value_counts(normalize=True)
    mimicry_org_frac = df["source_organism_short"].value_counts(normalize=True)

    common_orgs = mimicry_org_frac.head(10).index
    org_fc = []
    for o in common_orgs:
        input_frac = input_org_frac.get(o, 0)
        mimicry_frac = mimicry_org_frac.get(o, 0)
        if input_frac > 0:
            org_fc.append(mimicry_frac / input_frac)
        else:
            org_fc.append(np.nan)

    org_fc_series = pd.Series(org_fc, index=common_orgs).dropna().sort_values(ascending=True)

    ax = axes[1]
    colors = ["#2E7D32" if v >= 1 else "#C8E6C9" for v in org_fc_series.values]
    ax.barh(range(len(org_fc_series)), org_fc_series.values, color=colors, edgecolor="white")
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1, alpha=0.7)
    ax.set_yticks(range(len(org_fc_series)))
    ax.set_yticklabels(org_fc_series.index, fontsize=9)
    ax.set_xlabel("Fold enrichment (mimicry / input)")
    ax.set_title("Organism enrichment in strong mimicry")

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "07_enrichment_comparison.png"))
    plt.close()
    print("  ✅ 07_enrichment_comparison.png")
else:
    print("  ⚠️ Input file not found, skipping enrichment plot")

# ====================================================
# Figure 8: Autoimmune risk HLA analysis
# ====================================================
print("Plotting autoimmune risk HLA analysis...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

risk_counts = df["hla_risk"].value_counts()
ax = axes[0]
risk_colors = {
    "Predisposing": "#D32F2F",
    "Protective": "#1976D2",
    "Not annotated": "#BDBDBD"
}
ax.bar(
    range(len(risk_counts)),
    risk_counts.values,
    color=[risk_colors.get(r, "#999") for r in risk_counts.index],
    edgecolor="white"
)
ax.set_xticks(range(len(risk_counts)))
ax.set_xticklabels(risk_counts.index, fontsize=10)
ax.set_ylabel("Number of mimicry pairs")
ax.set_title("Mimicry pairs by HLA autoimmune risk status")
for i, (status, count) in enumerate(risk_counts.items()):
    ax.text(i, count + max(risk_counts) * 0.01, str(count),
            ha="center", fontsize=10, fontweight="bold")

predisposing_pairs = df[df["hla_risk"] == "Predisposing"]
if not predisposing_pairs.empty:
    disease_records = []
    for _, row in predisposing_pairs.iterrows():
        allele = row["mhc_allele"]
        diseases = predisposing_diseases.get(allele, ["Unknown"])
        for disease in diseases:
            disease_records.append({
                "disease": disease,
                "allele": allele,
                "pair_id": row["pair_id"]
            })

    disease_df = pd.DataFrame(disease_records)
    disease_counts = disease_df["disease"].value_counts().head(15)

    ax = axes[1]
    ax.barh(
        range(len(disease_counts)),
        disease_counts.values,
        color="#D32F2F", alpha=0.7, edgecolor="white"
    )
    ax.set_yticks(range(len(disease_counts)))
    ax.set_yticklabels(disease_counts.index, fontsize=9)
    ax.set_xlabel("Number of mimicry pairs")
    ax.set_title("Autoimmune diseases linked to predisposing HLA alleles")
    ax.invert_yaxis()
else:
    axes[1].text(0.5, 0.5, "No predisposing alleles\nin mimicry candidates",
                 ha="center", va="center", fontsize=12, transform=axes[1].transAxes)
    axes[1].set_title("Autoimmune diseases")

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "08_autoimmune_risk_hla.png"))
plt.close()
print("  ✅ 08_autoimmune_risk_hla.png")

# ====================================================
# Summary statistics
# ====================================================
print(f"\n{'=' * 60}")
print("Summary Statistics")
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
print(f"  Median viral BA_Rank: {df['viral_rank_BA'].median():.4f}")
print(f"  Median human BA_Rank: {df['human_rank_BA'].median():.4f}")

print(f"\n  Figures saved to: {FIG_DIR}")
print(f"\n  NOTE: Domain and protein ORA figures are generated by:")
print(f"    python compute_domain_enrichment.py")
print(f"    python compute_protein_enrichment.py")
print(f"\n✅ Done.")