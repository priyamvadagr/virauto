#!/usr/bin/env python3
"""
compute_protein_enrichment.py

Compute per-protein enrichment/depletion of mimicry candidates (mimetopes)
using a Binomial null model, separately for each peptide length k.
Results collapsed across k: max fold enrichment, min q-value.

Model (per k):
    rho_k = mean(Ni_k) / mean(Li_k)        # proteome-wide mimetope rate
    P(Ni_k >= n | Li_k, rho_k) = 1 - Pbinom(n-1 | Li_k, rho_k)   [enrichment]
    log2_odds_k = log2(Ni_k / (rho_k * Li_k))

Input:
    - mimetopes_per_protein.tsv (Li_k and Ni_k per protein,
      from compute_mimetopes_per_protein.py)

Output:
    - protein_enrichment_per_k.tsv (all per-k results)
    - protein_enrichment_collapsed.tsv (collapsed across k)
    - rho_summary.tsv
    - Figures: dot plot + volcano

Usage:
    python compute_protein_enrichment.py \
        --mimetopes mimetopes_per_protein.tsv \
        --outdir results/protein_enrichment \
        --fdr 0.05
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests

from adjustText import adjust_text

sns.set_style("whitegrid")
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
})


# ================================================================
# Parsing
# ================================================================

def parse_proteome_fasta(fasta_path):
    print(f"  Parsing {fasta_path} ...")
    ac2info = {}
    current_ac = None
    current_name = None
    current_len = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_ac is not None:
                    ac2info[current_ac] = {"name": current_name, "seq_len": current_len}
                parts = line[1:].split("|")
                if len(parts) >= 3:
                    current_ac = parts[1].strip()
                    name_field = parts[2].split(" ", 1)
                    current_name = name_field[1].split(" OS=")[0].strip() if len(name_field) > 1 else current_ac
                else:
                    current_ac = None
                    current_name = None
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_ac is not None:
        ac2info[current_ac] = {"name": current_name, "seq_len": current_len}
    print(f"    {len(ac2info):,} proteins parsed")
    return ac2info


# ================================================================
# Binomial test per k
# ================================================================

def test_enrichment_per_k(proteome_df, k_values, fdr_threshold):
    """
    Binomial enrichment test per protein per k.
    """
    print(f"\n  Running binomial enrichment tests ...")

    all_results = []
    rho_summary = []

    for k in k_values:
        li_col = f"Li_{k}"
        ni_col = f"Ni_{k}"

        mask = proteome_df[li_col] > 0
        sub = proteome_df.loc[mask].copy()

        if sub.empty:
            continue

        Li = sub[li_col].values
        Ni = sub[ni_col].values

        # Estimate rho (ratio of means)
        rho = Ni.mean() / Li.mean()

        print(f"    k={k}: {mask.sum():,} proteins, rho={rho:.4e}, "
              f"{(Ni > 0).sum():,} with Ni>0")

        expected = rho * Li

        # Log2 odds
        with np.errstate(divide="ignore", invalid="ignore"):
            log2_odds = np.log2(Ni / expected)
        log2_odds[~np.isfinite(log2_odds) & (Ni == 0)] = np.nan

        # P-values (enrichment only — depletion has no power here)
        pval_enrich = 1.0 - binom.cdf(Ni - 1, Li, rho)
        pval_enrich = np.clip(pval_enrich, 0, 1)

        # BH correction
        _, qval_enrich, _, _ = multipletests(
            pval_enrich, alpha=fdr_threshold, method="fdr_bh"
        )

        # Calls
        call = np.where(qval_enrich < fdr_threshold, "MEP", "NS")

        for i in range(len(sub)):
            if Ni[i] == 0:
                continue  # Only store proteins with mimetopes for collapsed output
            all_results.append({
                "uniprot_ac": sub.iloc[i]["uniprot_ac"],
                "protein_name": sub.iloc[i]["protein_name"],
                "seq_len": int(sub.iloc[i]["seq_len"]),
                "k": k,
                "Li": int(Li[i]),
                "Ni": int(Ni[i]),
                "rho": rho,
                "expected": round(expected[i], 4),
                "log2_odds": round(log2_odds[i], 4) if np.isfinite(log2_odds[i]) else np.nan,
                "fold_enrichment": round(Ni[i] / expected[i], 4) if expected[i] > 0 else np.inf,
                "pval": pval_enrich[i],
                "qval": qval_enrich[i],
                "call": call[i],
            })

        n_mep = (call == "MEP").sum()
        print(f"      {n_mep} MEP (FDR < {fdr_threshold})")

        rho_summary.append({
            "k": k,
            "rho": rho,
            "n_proteins_tested": int(mask.sum()),
            "n_with_mimetopes": int((Ni > 0).sum()),
            "n_MEP": int(n_mep),
        })

    return pd.DataFrame(all_results), pd.DataFrame(rho_summary)


# ================================================================
# Collapse across k
# ================================================================

def collapse_across_k(per_k_df, fdr_threshold):
    """
    Per protein: max fold enrichment, min q-value across k.
    """
    print(f"\n  Collapsing across k-mers per protein ...")

    collapsed = []
    for (ac, name), group in per_k_df.groupby(["uniprot_ac", "protein_name"]):
        best_fe_row = group.loc[group["fold_enrichment"].idxmax()]
        min_qval = group["qval"].min()
        total_ni = group["Ni"].sum()

        collapsed.append({
            "uniprot_ac": ac,
            "protein_name": name,
            "seq_len": int(best_fe_row["seq_len"]),
            "max_fold_enrichment": best_fe_row["fold_enrichment"],
            "best_k_fe": int(best_fe_row["k"]),
            "Ni_at_best_k": int(best_fe_row["Ni"]),
            "total_Ni_all_k": int(total_ni),
            "min_qval": min_qval,
            "min_qval_k": int(group.loc[group["qval"].idxmin(), "k"]),
            "n_k_tested": len(group),
            "n_k_MEP": int((group["call"] == "MEP").sum()),
            "significant": min_qval < fdr_threshold,
        })

    collapsed_df = pd.DataFrame(collapsed).sort_values("min_qval")
    n_sig = collapsed_df["significant"].sum()
    print(f"    {len(collapsed_df)} proteins with mimetopes, {n_sig} MEP")

    return collapsed_df


# ================================================================
# Plotting
# ================================================================

def plot_protein_ora(collapsed_df, outdir, fdr_threshold):
    fig_dir = os.path.join(outdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    sig = collapsed_df[collapsed_df["significant"]].copy()
    insig = collapsed_df[~collapsed_df["significant"]]

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # --- Dot plot: top significant proteins ---
    ax = axes[0]
    top = sig.copy()
    if not top.empty:
        top["_log2_fe"] = np.log2(top["max_fold_enrichment"].clip(lower=1e-10))
        top["_score"] = top["_log2_fe"] * (-np.log10(top["min_qval"].clip(lower=1e-300)))
        top = top.nlargest(25, "_score")

    if top.empty:
        ax.text(0.5, 0.5, f"No significant proteins\nat FDR < {fdr_threshold}",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_title("Protein ORA: mimetope-enriched proteins")
    else:
        neg_log_q = -np.log10(top["min_qval"].clip(lower=1e-300))
        norm = plt.Normalize(vmin=1.3, vmax=10)
        cmap = plt.cm.RdYlBu_r

        # Size proportional to mimetope count
        n_vals = top["total_Ni_all_k"].values.astype(float)
        min_dot = 50
        max_dot = 500
        n_max = n_vals.max()
        dot_sizes = min_dot + (n_vals / n_max) * (max_dot - min_dot)

        scatter = ax.scatter(
            np.log2(top["max_fold_enrichment"].clip(lower=1e-10)),
            range(len(top)),
            s=dot_sizes,
            c=neg_log_q.values,
            cmap=cmap,
            norm=norm,
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )

        labels = top["protein_name"].tolist()
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("log2(max fold enrichment)")
        ax.set_title(f"Top 25 mimetope-enriched proteins (FDR < {fdr_threshold})")

        cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
        cbar.set_label("-log10(min q-value)")

        # Size legend
        legend_sizes = sorted(set([
            max(1, int(n_max * 0.1)),
            max(1, int(n_max * 0.5)),
            int(n_max),
        ]))
        legend_dots = []
        for s in legend_sizes:
            ds = min_dot + (s / n_max) * (max_dot - min_dot)
            legend_dots.append(
                ax.scatter([], [], s=ds, c="gray", edgecolors="white",
                           linewidths=0.5, label=f"{s} mimetopes")
            )
        ax.legend(
            handles=legend_dots,
            title="Mimetope count",
            loc="lower right",
            fontsize=7,
            title_fontsize=8,
            frameon=True,
        )

    # --- Volcano ---
    ax = axes[1]
    ax.scatter(
        np.log2(insig["max_fold_enrichment"].clip(lower=1e-10)),
        -np.log10(insig["min_qval"].clip(lower=1e-300)),
        s=20, alpha=0.3, color="#BDBDBD",
    )
    if not sig.empty:
        ax.scatter(
            np.log2(sig["max_fold_enrichment"].clip(lower=1e-10)),
            -np.log10(sig["min_qval"].clip(lower=1e-300)),
            s=40, alpha=0.8, color="#D32F2F",
            label=f"FDR < {fdr_threshold}",
        )

        sig_score = sig.copy()
        sig_score["_log2_fe"] = np.log2(sig_score["max_fold_enrichment"].clip(lower=1e-10))
        sig_score["score"] = sig_score["_log2_fe"] * (
            -np.log10(sig_score["min_qval"].clip(lower=1e-300))
        )
        top_labels = sig_score.nlargest(12, "score")
        texts = []
        for _, row in top_labels.iterrows():
            x = np.log2(max(row["max_fold_enrichment"], 1e-10))
            y = -np.log10(max(row["min_qval"], 1e-300))
            name = row["protein_name"]
            if len(name) > 45:
                name = name[:42] + "..."
            texts.append(ax.text(x, y, name, fontsize=6, color="#333"))

        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="#999", lw=0.5),
            expand=(1.5, 1.5),
            force_text=(0.5, 1.0),
            force_points=(0.3, 0.5),
        )

    ax.axvline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axhline(-np.log10(fdr_threshold), color="gray", linestyle=":",
               linewidth=0.8, label=f"FDR = {fdr_threshold}")
    ax.set_xlabel("log2(max fold enrichment)")
    ax.set_ylabel("-log10(min q-value)")
    ax.set_title("Protein ORA volcano\n(binomial test, collapsed across k)")
    ax.legend(fontsize=8)

    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "protein_ora_collapsed.png")
    plt.savefig(fig_path)
    plt.close()
    print(f"    Saved {fig_path}")


# ================================================================
# Main
# ================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Binomial enrichment/depletion of mimetopes per protein"
    )
    parser.add_argument("--mimetopes", required=True,
                        help="mimetopes_per_protein.tsv (Li and Ni per protein)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="BH q-value threshold (default: 0.05)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --- Load precomputed data ---
    print("\n" + "=" * 60)
    print("Loading data")
    print("=" * 60)

    print(f"  Loading {args.mimetopes} ...")
    proteome_df = pd.read_csv(args.mimetopes, sep="\t")
    print(f"    {len(proteome_df):,} proteins loaded")

    # Detect k values from Li columns
    k_values = sorted(
        int(c.split("_")[1]) for c in proteome_df.columns if c.startswith("Li_")
    )
    print(f"    Peptide lengths detected: {k_values}")

    # Standardize column names
    col_map = {}
    if "hu_prot_id" in proteome_df.columns and "uniprot_ac" not in proteome_df.columns:
        col_map["hu_prot_id"] = "uniprot_ac"
    if "description" in proteome_df.columns and "protein_name" not in proteome_df.columns:
        col_map["description"] = "protein_name"
    if "length" in proteome_df.columns and "seq_len" not in proteome_df.columns:
        col_map["length"] = "seq_len"
    if col_map:
        proteome_df = proteome_df.rename(columns=col_map)
        print(f"    Renamed columns: {col_map}")

    for k in k_values:
        n_with = (proteome_df[f"Ni_{k}"] > 0).sum()
        print(f"    k={k}: {n_with:,} proteins with Ni > 0")

    # --- Binomial tests ---
    print("\n" + "=" * 60)
    print("Binomial enrichment tests")
    print("=" * 60)

    per_k_df, rho_df = test_enrichment_per_k(proteome_df, k_values, args.fdr)

    if per_k_df.empty:
        print("\n  No results. Check input data.")
        return

    # --- Collapse across k ---
    collapsed_df = collapse_across_k(per_k_df, args.fdr)

    # --- Save ---
    print(f"\n{'=' * 60}")
    print("Saving results")
    print(f"{'=' * 60}")

    per_k_path = os.path.join(args.outdir, "protein_enrichment_per_k.tsv")
    per_k_df.sort_values(["k", "qval"]).to_csv(per_k_path, sep="\t", index=False)
    print(f"  Per-k results: {per_k_path}")

    collapsed_path = os.path.join(args.outdir, "protein_enrichment_collapsed.tsv")
    collapsed_df.to_csv(collapsed_path, sep="\t", index=False)
    print(f"  Collapsed results: {collapsed_path}")

    rho_path = os.path.join(args.outdir, "rho_summary.tsv")
    rho_df.to_csv(rho_path, sep="\t", index=False)
    print(f"  Rho summary: {rho_path}")
    print(rho_df.to_string(index=False))

    # Summary
    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")
    print(f"  Proteins with mimetopes: {len(collapsed_df)}")
    print(f"  MEP (FDR < {args.fdr}): {collapsed_df['significant'].sum()}")
    top5 = collapsed_df[collapsed_df["significant"]].head(5)
    if not top5.empty:
        print(f"\n  Top 5 by q-value:")
        for _, row in top5.iterrows():
            print(f"    {row['uniprot_ac']} — {row['protein_name']}: "
                  f"FE={row['max_fold_enrichment']:.1f} (k={int(row['best_k_fe'])}), "
                  f"n={int(row['total_Ni_all_k'])}, "
                  f"q={row['min_qval']:.2e}")

    # --- Plots ---
    print(f"\n{'=' * 60}")
    print("Generating plots")
    print(f"{'=' * 60}")
    plot_protein_ora(collapsed_df, args.outdir, args.fdr)

    print(f"\n✅ Done.")


if __name__ == "__main__":
    main()