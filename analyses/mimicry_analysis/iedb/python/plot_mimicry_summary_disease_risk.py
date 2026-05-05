#!/usr/bin/env python3
"""
======================================================================
Script: plot_per_disease_mimicry.py
Description:
    For each autoimmune disease, plot:
      1. HLA conversion rate bar chart (disease alleles + motif-similar)
      2. HLA × organism dot plot (disease alleles + motif-similar)

    Motif-similar alleles are identified from PSSM correlation.

Input:
    --strong-mimicry    : strong mimicry candidates CSV
    --blast-input       : BLAST filtered pairs CSV
    --pair-map          : pair_id_mapping.csv.gz
    --hla-risk          : HLA autoimmunity classification file
    --pssm-matrix       : PSSM matrix CSV (allele × features)
    --out-dir           : output directory for figures

Optional:
    --corr-threshold    : correlation threshold for motif similarity (default: 0.8)
    --top-organisms     : number of top organisms to show (default: 10)

Usage:
    python plot_per_disease_mimicry.py \
        --strong-mimicry iedb_mhci_strong_mimicry.csv.gz \
        --blast-input iedb_mhci_pairs_4digit_hla_swissprot.csv.gz \
        --pair-map pair_id_mapping.csv.gz \
        --hla-risk HLA_autoimmunity_classification.txt \
        --pssm-matrix pssm_matrix_9mer.csv \
        --out-dir results/per_disease_plots/
======================================================================
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from scipy.stats import binom, fisher_exact
from statsmodels.stats.multitest import multipletests

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
# Helpers
# ====================================================

def qval_to_stars(q):
    if q < 0.001:
        return "***"
    elif q < 0.01:
        return "**"
    elif q < 0.05:
        return "*"
    return ""


def convert_hla_name(name):
    """HLA_A_0201 → HLA-A*02:01"""
    parts = name.split("_")
    if len(parts) == 3:
        locus = parts[1]
        digits = parts[2]
        if len(digits) == 4:
            return f"HLA-{locus}*{digits[:2]}:{digits[2:]}"
    return name


def standard_to_pssm_name(allele):
    """HLA-A*02:01 → A0201"""
    m = allele.replace("HLA-", "").replace("*", "").replace(":", "")
    return m


def pssm_to_standard_name(pssm_name):
    """A0201 → HLA-A*02:01"""
    # Parse locus and digits
    # A0201 → A, 0201
    # B0702 → B, 0702
    # C0602 → C, 0602
    for prefix_len in [1, 2, 3]:
        locus = pssm_name[:prefix_len]
        digits = pssm_name[prefix_len:]
        if len(digits) == 4 and digits.isdigit():
            return f"HLA-{locus}*{digits[:2]}:{digits[2:]}"
    return pssm_name


# ====================================================
# Load data
# ====================================================

def load_all_data(args):
    print("=" * 60)
    print("Loading data")
    print("=" * 60)

    # Strong mimicry
    print(f"  Loading strong mimicry: {args.strong_mimicry}")
    strong = pd.read_csv(args.strong_mimicry, low_memory=False)
    strong["source_organism_short"] = (
        strong["source_organism"]
        .str.replace(r"\s*\(.*?\)", "", regex=True)
        .str.strip()
    )
    print(f"    {len(strong):,} records")

    # BLAST input + pair IDs
    print(f"  Loading BLAST input: {args.blast_input}")
    blast = pd.read_csv(args.blast_input, low_memory=False)
    pair_map = pd.read_csv(args.pair_map, low_memory=False)

    blast["pair_key"] = (
        blast["structure_id"].astype(str) + "_" +
        blast["hu_prot_id"].astype(str) + "_" +
        blast["viral_sequence"].astype(str) + "_" +
        blast["human_mimic_sequence"].astype(str)
    )
    key_to_pid = pair_map.drop_duplicates("pair_key").set_index("pair_key")["pair_id"]
    blast["pair_id"] = blast["pair_key"].map(key_to_pid)
    blast = blast.dropna(subset=["pair_id"])
    blast["source_organism_short"] = (
        blast["source_organism"]
        .str.replace(r"\s*\(.*?\)", "", regex=True)
        .str.strip()
    )
    print(f"    {len(blast):,} pairs with pair_ids")

    # HLA risk annotations
    print(f"  Loading HLA risk: {args.hla_risk}")
    risk_df = pd.read_csv(args.hla_risk, sep=r"\s+", engine="python")
    risk_df.columns = risk_df.columns.str.strip('"')
    for col in risk_df.select_dtypes(include="object").columns:
        risk_df[col] = risk_df[col].str.strip('"')
    risk_df["hla_standard"] = risk_df["HLA"].apply(convert_hla_name)

    # Class I only
    risk_ci = risk_df[risk_df["HLA_Class"] == "ClassI"].copy()
    print(f"    Class I entries: {len(risk_ci)}")

    # Build disease → alleles mapping
    diseases = {}
    for _, row in risk_ci[risk_ci["Association"] == "Predisposing"].iterrows():
        d = row["Disease"]
        a = row["hla_standard"]
        if d not in diseases:
            diseases[d] = {"predisposing": set(), "protective": set()}
        diseases[d]["predisposing"].add(a)

    for _, row in risk_ci[risk_ci["Association"] == "Protective"].iterrows():
        d = row["Disease"]
        a = row["hla_standard"]
        if d not in diseases:
            diseases[d] = {"predisposing": set(), "protective": set()}
        diseases[d]["protective"].add(a)

    print(f"    Diseases: {len(diseases)}")

    # PSSM matrix
    print(f"  Loading PSSM matrix: {args.pssm_matrix}")
    pssm = pd.read_csv(args.pssm_matrix, index_col=0)
    print(f"    {pssm.shape[0]} alleles × {pssm.shape[1]} features")

    return strong, blast, diseases, pssm


# ====================================================
# Motif similarity
# ====================================================

def find_similar_alleles(pssm_df, target_alleles, corr_threshold):
    """
    For each target allele, find alleles with Pearson correlation
    above threshold in the PSSM matrix.

    Parameters
    ----------
    pssm_df : DataFrame
        Allele × features (PSSM names like A0201)
    target_alleles : set of str
        Standard HLA names (HLA-A*02:01)
    corr_threshold : float

    Returns
    -------
    similar : dict
        standard_allele → list of (similar_standard_allele, correlation)
    all_similar_alleles : set of str
        All motif-similar alleles (standard names)
    """
    # Convert target alleles to PSSM names
    target_pssm = {}
    for a in target_alleles:
        pn = standard_to_pssm_name(a)
        if pn in pssm_df.index:
            target_pssm[a] = pn

    if not target_pssm:
        return {}, set()

    # Compute correlation matrix
    corr_matrix = pssm_df.T.corr()

    similar = {}
    all_similar = set()

    for std_name, pssm_name in target_pssm.items():
        corrs = corr_matrix.loc[pssm_name].drop(pssm_name)
        above = corrs[corrs >= corr_threshold].sort_values(ascending=False)

        hits = []
        for sim_pssm, r in above.items():
            sim_std = pssm_to_standard_name(sim_pssm)
            hits.append((sim_std, round(r, 3)))
            all_similar.add(sim_std)

        similar[std_name] = hits

    return similar, all_similar


# ====================================================
# Conversion rate computation
# ====================================================

def compute_conversion_rates(strong, blast):
    """Compute conversion rates per HLA and per (HLA, organism)."""

    # Per HLA
    input_per_hla = blast.groupby("mhc_allele")["pair_id"].nunique().rename("input_pairs")
    strong_per_hla = strong.groupby("mhc_allele")["pair_id"].nunique().rename("strong_pairs")

    hla_conv = input_per_hla.to_frame().join(strong_per_hla, how="left").fillna(0)
    hla_conv["strong_pairs"] = hla_conv["strong_pairs"].astype(int)
    hla_conv["fraction"] = hla_conv["strong_pairs"] / hla_conv["input_pairs"]

    # Per (HLA, organism)
    input_per_cell = (
        blast.groupby(["mhc_allele", "source_organism"])["pair_id"]
        .nunique().rename("input_pairs")
    )
    strong_per_cell = (
        strong.groupby(["mhc_allele", "source_organism"])["pair_id"]
        .nunique().rename("strong_pairs")
    )

    cell_conv = input_per_cell.to_frame().join(strong_per_cell, how="left").fillna(0)
    cell_conv["strong_pairs"] = cell_conv["strong_pairs"].astype(int)
    cell_conv["fraction"] = cell_conv["strong_pairs"] / cell_conv["input_pairs"]

    return hla_conv, cell_conv


# ====================================================
# Per-disease plotting
# ====================================================

def plot_disease(disease_name, pred_alleles, prot_alleles,
                 similar_map, all_similar,
                 hla_conv, cell_conv, strong,
                 out_dir, top_n_org=10):
    """Generate bar chart + dot plot for one disease."""

    safe_name = disease_name.replace(" ", "_").replace("/", "_").lower()
    disease_dir = os.path.join(out_dir, safe_name)
    os.makedirs(disease_dir, exist_ok=True)

    # Build allele sets
    # Categories: predisposing, protective, motif-similar
    pred_set = pred_alleles
    prot_set = prot_alleles
    sim_set = all_similar - pred_set - prot_set

    # All alleles to plot (that have data)
    all_alleles = set()
    for a in pred_set | prot_set | sim_set:
        if a in hla_conv.index:
            all_alleles.add(a)

    if not all_alleles:
        print(f"    ⚠️ No alleles with data for {disease_name}")
        return

    # Sort: predisposing first, then similar, then protective
    def sort_key(a):
        if a in pred_set:
            return (0, -hla_conv.loc[a, "input_pairs"])
        elif a in sim_set:
            return (1, -hla_conv.loc[a, "input_pairs"])
        else:
            return (2, -hla_conv.loc[a, "input_pairs"])

    sorted_alleles = sorted(all_alleles, key=sort_key)

    # ---- Figure A: HLA conversion rate bar chart ----
    fig, ax = plt.subplots(figsize=(11, max(5, len(sorted_alleles) * 0.35)))

    # Global rate for binomial test
    global_rate = hla_conv["strong_pairs"].sum() / hla_conv["input_pairs"].sum()

    # Binomial test for plotted alleles
    plot_data = hla_conv.loc[[a for a in sorted_alleles if a in hla_conv.index]].copy()
    pvals = []
    for _, row in plot_data.iterrows():
        if row["input_pairs"] == 0:
            pvals.append(1.0)
        else:
            p = 1.0 - binom.cdf(int(row["strong_pairs"]) - 1,
                                 int(row["input_pairs"]), global_rate)
            pvals.append(min(p, 1.0))
    if pvals:
        _, qvals, _, _ = multipletests(pvals, method="fdr_bh")
    else:
        qvals = []
    plot_data["qval"] = qvals
    plot_data["stars"] = [qval_to_stars(q) for q in qvals]

    # Colors
    bar_colors = []
    for a in sorted_alleles:
        if a in pred_set:
            bar_colors.append("#D32F2F")
        elif a in sim_set:
            bar_colors.append("#FF9800")
        elif a in prot_set:
            bar_colors.append("#1976D2")
        else:
            bar_colors.append("#BDBDBD")

    ax.barh(
        range(len(sorted_alleles)),
        [plot_data.loc[a, "fraction"] if a in plot_data.index else 0
         for a in reversed(sorted_alleles)],
        color=list(reversed(bar_colors)),
        edgecolor="white",
    )

    ax.set_yticks(range(len(sorted_alleles)))
    labels_rev = list(reversed(sorted_alleles))
    ax.set_yticklabels(labels_rev, fontsize=8)

    # Style labels
    for i, label in enumerate(ax.get_yticklabels()):
        a = labels_rev[i]
        if a in pred_set:
            label.set_color("#D32F2F")
            label.set_fontweight("bold")
        elif a in sim_set:
            label.set_color("#FF9800")
        elif a in prot_set:
            label.set_color("#1976D2")

    # Annotate
    for i, a in enumerate(reversed(sorted_alleles)):
        if a in plot_data.index:
            row = plot_data.loc[a]
            stars = row["stars"]

            # Build annotation
            ann = f"{int(row['strong_pairs'])}/{int(row['input_pairs'])} {stars}"

            # Add similarity info for motif-similar alleles
            if a in sim_set:
                # Find which disease allele it's similar to
                for da, sims in similar_map.items():
                    for sim_a, r in sims:
                        if sim_a == a:
                            ann += f" (r={r:.2f} w/ {da})"
                            break

            ax.text(
                row["fraction"] + 0.002, i,
                ann, va="center", fontsize=6, color="#333",
            )

    ax.set_xlabel("Conversion rate (strong mimicry / input pairs)")
    ax.set_title(f"{disease_name}\nHLA conversion rates (predisposing + motif-similar)")

    legend_elements = [
        Patch(facecolor="#D32F2F", label="Predisposing"),
        Patch(facecolor="#FF9800", label="Motif-similar"),
        Patch(facecolor="#1976D2", label="Protective"),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc="lower right")

    plt.tight_layout()
    plt.savefig(os.path.join(disease_dir, "hla_conversion_rate.png"))
    plt.close()

    # ---- Figure B: HLA × organism dot plot ----
    top_orgs_full = strong["source_organism"].value_counts().head(top_n_org).index.tolist()
    top_orgs_short = [
        strong[strong["source_organism"] == o]["source_organism_short"].iloc[0]
        for o in top_orgs_full
    ]

    # Build matrices
    frac_matrix = pd.DataFrame(0.0, index=top_orgs_short, columns=sorted_alleles)
    input_matrix = pd.DataFrame(0, index=top_orgs_short, columns=sorted_alleles)
    stars_matrix = pd.DataFrame("", index=top_orgs_short, columns=sorted_alleles, dtype=str)

    # Fisher's exact test for HLA × organism interaction
    # For each cell:
    #           | this organism | other organisms |
    # this HLA  |      a        |       b         |
    # other HLA |      c        |       d         |
    strong_per_hla_all = strong.groupby("mhc_allele")["pair_id"].nunique()
    strong_per_org_all = strong.groupby("source_organism")["pair_id"].nunique()
    total_strong_all = strong["pair_id"].nunique()

    # Fill matrices and compute interaction p-values
    cell_pvals = []
    cell_indices = []

    for org_full, org_short in zip(top_orgs_full, top_orgs_short):
        for allele in sorted_alleles:
            key = (allele, org_full)
            if key in cell_conv.index:
                row = cell_conv.loc[key]
                frac_matrix.loc[org_short, allele] = row["fraction"]
                input_matrix.loc[org_short, allele] = int(row["input_pairs"])

                k = int(row["strong_pairs"])
                if k == 0:
                    cell_pvals.append(1.0)
                else:
                    a = k
                    b = max(strong_per_hla_all.get(allele, 0) - a, 0)
                    c = max(strong_per_org_all.get(org_full, 0) - a, 0)
                    d = max(total_strong_all - a - b - c, 0)
                    _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
                    cell_pvals.append(pval)
                cell_indices.append((org_short, allele))
            else:
                cell_pvals.append(1.0)
                cell_indices.append((org_short, allele))

    # BH correction
    if cell_pvals:
        _, cell_qvals, _, _ = multipletests(cell_pvals, method="fdr_bh")
        for (org_short, allele), q in zip(cell_indices, cell_qvals):
            stars_matrix.loc[org_short, allele] = qval_to_stars(q)

    # Plot
    fig, ax = plt.subplots(figsize=(max(8, len(sorted_alleles) * 0.7), 7))

    vmax = min(frac_matrix.max().max() * 1.1, 1.0)
    if vmax == 0:
        vmax = 1.0
    norm = plt.Normalize(vmin=0, vmax=vmax)
    cmap = plt.cm.YlOrRd

    max_input = max(input_matrix.max().max(), 1)
    min_dot, max_dot = 30, 500

    for i, org_short in enumerate(top_orgs_short):
        for j, allele in enumerate(sorted_alleles):
            frac = frac_matrix.loc[org_short, allele]
            n_input = input_matrix.loc[org_short, allele]

            if n_input == 0:
                continue

            dot_size = min_dot + (n_input / max_input) * (max_dot - min_dot)

            ax.scatter(
                j + 0.5, i + 0.5,
                s=dot_size,
                c=[cmap(norm(frac))],
                edgecolors="white",
                linewidths=0.5,
                zorder=3,
            )

            stars = stars_matrix.loc[org_short, allele]
            if stars:
                ax.text(
                    j + 0.5, i + 0.9,
                    stars,
                    ha="center", va="center",
                    fontsize=8, fontweight="bold", color="#333",
                    zorder=4,
                )

    # Grid
    for i in range(len(top_orgs_short) + 1):
        ax.axhline(i, color="#E0E0E0", linewidth=0.5, zorder=1)
    for j in range(len(sorted_alleles) + 1):
        ax.axvline(j, color="#E0E0E0", linewidth=0.5, zorder=1)

    ax.set_xlim(0, len(sorted_alleles))
    ax.set_ylim(0, len(top_orgs_short))
    ax.set_xticks([j + 0.5 for j in range(len(sorted_alleles))])
    ax.set_xticklabels(sorted_alleles, rotation=45, ha="right", fontsize=7)
    ax.set_yticks([i + 0.5 for i in range(len(top_orgs_short))])
    ax.set_yticklabels(top_orgs_short, fontsize=9)
    ax.invert_yaxis()

    # Color x labels
    for label in ax.get_xticklabels():
        a = label.get_text()
        if a in pred_set:
            label.set_color("#D32F2F")
            label.set_fontweight("bold")
        elif a in sim_set:
            label.set_color("#FF9800")
        elif a in prot_set:
            label.set_color("#1976D2")

    ax.set_xlabel("HLA allele")
    ax.set_ylabel("Source organism")
    ax.set_title(
        f"{disease_name}\nHLA × organism interaction "
        f"(stars = Fisher's exact FDR: * <.05, ** <.01, *** <.001)"
    )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Conversion rate")

    legend_sizes = sorted(set([
        max(1, int(max_input * 0.1)),
        max(1, int(max_input * 0.5)),
        int(max_input),
    ]))
    legend_dots = [
        ax.scatter([], [], s=min_dot + (s / max_input) * (max_dot - min_dot),
                   c="gray", edgecolors="white", linewidths=0.5,
                   label=f"{s:,} pairs")
        for s in legend_sizes
    ]
    ax.legend(
        handles=legend_dots, title="Input pairs",
        loc="upper left", bbox_to_anchor=(1.15, 1.0),
        labelspacing=1.5,
        fontsize=7, title_fontsize=8, frameon=True,
    )

    plt.tight_layout()
    plt.savefig(os.path.join(disease_dir, "hla_organism_dotplot.png"))
    plt.close()

    print(f"    ✅ {disease_name}: {len(sorted_alleles)} alleles "
          f"({len(pred_set & all_alleles)} pred, {len(sim_set & all_alleles)} sim, "
          f"{len(prot_set & all_alleles)} prot)")


# ====================================================
# Main
# ====================================================

def main():
    parser = argparse.ArgumentParser(
        description="Per-disease mimicry analysis with motif-similar HLA alleles"
    )
    parser.add_argument("--strong-mimicry", required=True)
    parser.add_argument("--blast-input", required=True)
    parser.add_argument("--pair-map", required=True)
    parser.add_argument("--hla-risk", required=True)
    parser.add_argument("--pssm-matrix", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--corr-threshold", type=float, default=0.8,
                        help="PSSM correlation threshold for motif similarity (default: 0.8)")
    parser.add_argument("--top-organisms", type=int, default=10,
                        help="Number of top organisms to show (default: 10)")
    args = parser.parse_args()

    # Load data
    strong, blast, diseases, pssm = load_all_data(args)

    # Compute conversion rates
    print(f"\n{'=' * 60}")
    print("Computing conversion rates")
    print("=" * 60)
    hla_conv, cell_conv = compute_conversion_rates(strong, blast)

    # Process each disease
    print(f"\n{'=' * 60}")
    print(f"Generating per-disease plots ({len(diseases)} diseases)")
    print(f"PSSM correlation threshold: {args.corr_threshold}")
    print("=" * 60)

    for disease_name, allele_sets in sorted(diseases.items()):
        pred = allele_sets["predisposing"]
        prot = allele_sets["protective"]

        if not pred:
            continue

        print(f"\n  {disease_name}:")
        print(f"    Predisposing: {sorted(pred)}")

        # Find motif-similar alleles
        similar_map, all_similar = find_similar_alleles(
            pssm, pred, args.corr_threshold
        )

        n_new = len(all_similar - pred - prot)
        if n_new > 0:
            print(f"    Motif-similar (r ≥ {args.corr_threshold}): {n_new} additional alleles")
            for da, sims in similar_map.items():
                if sims:
                    sim_str = ", ".join(f"{a} (r={r})" for a, r in sims[:5])
                    print(f"      {da}: {sim_str}")

        plot_disease(
            disease_name, pred, prot,
            similar_map, all_similar,
            hla_conv, cell_conv, strong,
            args.out_dir, top_n_org=args.top_organisms,
        )

    print(f"\n✅ Done. Figures in {args.out_dir}/")


if __name__ == "__main__":
    main()