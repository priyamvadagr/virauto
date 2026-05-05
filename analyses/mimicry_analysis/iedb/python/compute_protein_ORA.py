#!/usr/bin/env python3
"""
compute_functional_ora.py

Functional overrepresentation analysis for mimetope-enriched proteins (MEPs)
using gProfiler (g:GOSt). Tests GO:BP, GO:MF, GO:CC, and KEGG pathways.

Results are cached locally for reproducibility.

Input:
    - protein_enrichment_collapsed.tsv (from compute_protein_enrichment.py)
    - Swiss-Prot background list (from proteome FASTA)

Output:
    - functional_ora_results.tsv (all terms tested)
    - functional_ora_significant.tsv (FDR < threshold)
    - Cached API response (JSON)
    - Figures: dot plots per ontology

Dependencies:
    requests, pandas, matplotlib, seaborn

Usage:
    python compute_functional_ora.py \
        --enrichment protein_enrichment_collapsed.tsv \
        --proteome uniprot_human_all.fasta \
        --outdir results/functional_ora \
        --fdr 0.05
"""

import argparse
import json
import os
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
})

# Source name mapping for display
SOURCE_NAMES = {
    "GO:BP": "GO Biological Process",
    "GO:MF": "GO Molecular Function",
    "GO:CC": "GO Cellular Component",
    "KEGG": "KEGG Pathways",
}


# ================================================================
# Parse Swiss-Prot ACs from FASTA
# ================================================================

def get_swissprot_acs(fasta_path):
    print(f"  Parsing Swiss-Prot ACs from {fasta_path} ...")
    sp_acs = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">sp|"):
                ac = line.split("|")[1]
                sp_acs.append(ac)
    print(f"    {len(sp_acs):,} Swiss-Prot ACs")
    return sp_acs


# ================================================================
# Query gProfiler
# ================================================================

def query_gprofiler(query_acs, background_acs, cache_path):
    """
    Query gProfiler g:GOSt API for functional enrichment.
    Uses custom background (Swiss-Prot proteome).
    Caches results locally.
    """

    # Check cache
    if os.path.exists(cache_path):
        print(f"  Loading cached gProfiler results: {cache_path}")
        with open(cache_path) as f:
            return json.load(f)

    print(f"  Querying gProfiler API ...")
    print(f"    Query proteins: {len(query_acs):,}")
    print(f"    Background proteins: {len(background_acs):,}")
    print(f"    Sources: GO:BP, GO:MF, GO:CC, KEGG")

    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

    payload = {
        "organism": "hsapiens",
        "query": list(query_acs),
        "domain_scope": "custom",
        "background": list(background_acs),
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG"],
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "no_evidences": False,
        "all_results": True,
    }

    try:
        response = requests.post(url, json=payload, timeout=120)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        print(f"  ⚠️ gProfiler API error: {e}")
        return None

    # Cache
    with open(cache_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  Cached results: {cache_path}")

    return data


# ================================================================
# Parse gProfiler response
# ================================================================

def parse_gprofiler_results(data):
    """
    Parse gProfiler JSON response into a DataFrame.
    """
    if data is None or "result" not in data:
        print("  ⚠️ No results from gProfiler")
        return pd.DataFrame()

    results = data["result"]
    print(f"  Parsing {len(results):,} terms ...")

    rows = []
    for r in results:
        rows.append({
            "source": r.get("source", ""),
            "term_id": r.get("native", ""),
            "term_name": r.get("name", ""),
            "p_value": r.get("p_value", 1.0),
            "term_size": r.get("term_size", 0),
            "query_size": r.get("query_size", 0),
            "intersection_size": r.get("intersection_size", 0),
            "effective_domain_size": r.get("effective_domain_size", 0),
            "precision": r.get("precision", 0),
            "recall": r.get("recall", 0),
            "intersections": ",".join(
                str(x) for x in (r.get("intersections", []) or [])
                if not isinstance(x, list)
            ) if not any(isinstance(x, list) for x in (r.get("intersections", []) or [])) else
            ",".join(
                item for sublist in (r.get("intersections", []) or [])
                for item in (sublist if isinstance(sublist, list) else [sublist])
            ),
        })

    df = pd.DataFrame(rows)

    if df.empty:
        return df

    # Compute fold enrichment
    # FE = (intersection/query) / (term_size/domain_size)
    df["expected_fraction"] = df["term_size"] / df["effective_domain_size"]
    df["observed_fraction"] = df["intersection_size"] / df["query_size"]
    df["fold_enrichment"] = df["observed_fraction"] / df["expected_fraction"]
    df["fold_enrichment"] = df["fold_enrichment"].replace([np.inf, np.nan], 0)

    # Summary per source
    for source in df["source"].unique():
        sub = df[df["source"] == source]
        n_sig = (sub["p_value"] < 0.05).sum()
        print(f"    {source}: {len(sub)} terms tested, {n_sig} significant (FDR < 0.05)")

    return df


# ================================================================
# Plotting
# ================================================================

def plot_functional_ora(results_df, outdir, fdr_threshold):
    fig_dir = os.path.join(outdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    sources = ["GO:BP", "GO:MF", "GO:CC", "KEGG"]
    available = [s for s in sources if s in results_df["source"].unique()]

    if not available:
        print("  No significant terms to plot")
        return

    # --- Combined dot plot (top terms per source) ---
    fig, axes = plt.subplots(
        1, len(available),
        figsize=(6 * len(available), 8),
        squeeze=False,
    )

    for idx, source in enumerate(available):
        ax = axes[0, idx]
        sub = results_df[
            (results_df["source"] == source) &
            (results_df["p_value"] < fdr_threshold)
        ].copy()

        if sub.empty:
            ax.text(0.5, 0.5, f"No significant terms\n(FDR < {fdr_threshold})",
                    ha="center", va="center", transform=ax.transAxes, fontsize=11)
            ax.set_title(SOURCE_NAMES.get(source, source))
            continue

        # Top 15 by fold enrichment
        top = sub.nlargest(15, "fold_enrichment")

        neg_log_q = -np.log10(top["p_value"].clip(lower=1e-300))
        norm = plt.Normalize(vmin=0, vmax=max(neg_log_q.max(), 1))
        cmap = plt.cm.YlOrRd

        # Size = intersection count
        n_vals = top["intersection_size"].values.astype(float)
        min_dot = 50
        max_dot = 400
        n_max = max(n_vals.max(), 1)
        dot_sizes = min_dot + (n_vals / n_max) * (max_dot - min_dot)

        scatter = ax.scatter(
            top["fold_enrichment"],
            range(len(top)),
            s=dot_sizes,
            c=neg_log_q.values,
            cmap=cmap,
            norm=norm,
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )

        # Truncate long term names
        labels = [
            name[:50] + "..." if len(name) > 50 else name
            for name in top["term_name"]
        ]
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("Fold enrichment")
        ax.set_title(SOURCE_NAMES.get(source, source))

        cbar = plt.colorbar(scatter, ax=ax, pad=0.02, shrink=0.6)
        cbar.set_label("-log10(FDR)", fontsize=8)

        # Size legend (only on first panel)
        if idx == 0:
            legend_sizes = sorted(set([
                max(1, int(n_max * 0.2)),
                max(1, int(n_max * 0.6)),
                int(n_max),
            ]))
            legend_dots = []
            for s in legend_sizes:
                ds = min_dot + (s / n_max) * (max_dot - min_dot)
                legend_dots.append(
                    ax.scatter([], [], s=ds, c="gray", edgecolors="white",
                               linewidths=0.5, label=f"{s} proteins")
                )
            ax.legend(
                handles=legend_dots,
                title="MEPs in term",
                loc="lower right",
                fontsize=6,
                title_fontsize=7,
                frameon=True,
            )

    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "functional_ora_dotplot.png")
    plt.savefig(fig_path)
    plt.close()
    print(f"    Saved {fig_path}")

    # --- Individual volcano per source ---
    for source in available:
        sub = results_df[results_df["source"] == source].copy()
        sig = sub[sub["p_value"] < fdr_threshold]
        insig = sub[sub["p_value"] >= fdr_threshold]

        fig, ax = plt.subplots(figsize=(8, 6))

        ax.scatter(
            insig["fold_enrichment"],
            -np.log10(insig["p_value"].clip(lower=1e-300)),
            s=15, alpha=0.3, color="#BDBDBD",
        )
        if not sig.empty:
            ax.scatter(
                sig["fold_enrichment"],
                -np.log10(sig["p_value"].clip(lower=1e-300)),
                s=30, alpha=0.8, color="#D32F2F",
                label=f"FDR < {fdr_threshold}",
            )

            sig_score = sig.copy()
            sig_score["score"] = sig_score["fold_enrichment"] * (
                -np.log10(sig_score["p_value"].clip(lower=1e-300))
            )
            for _, row in sig_score.nlargest(8, "score").iterrows():
                name = row["term_name"]
                if len(name) > 40:
                    name = name[:40] + "..."
                ax.annotate(
                    name,
                    (row["fold_enrichment"],
                     -np.log10(max(row["p_value"], 1e-300))),
                    fontsize=6, xytext=(4, 2), textcoords="offset points",
                    color="#333",
                )

        ax.axvline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.axhline(-np.log10(fdr_threshold), color="gray", linestyle=":",
                   linewidth=0.8, label=f"FDR = {fdr_threshold}")
        ax.set_xlabel("Fold enrichment")
        ax.set_ylabel("-log10(FDR)")
        ax.set_title(f"{SOURCE_NAMES.get(source, source)} — Functional ORA volcano")
        ax.legend(fontsize=8)

        plt.tight_layout()
        safe_source = source.replace(":", "_")
        fig_path = os.path.join(fig_dir, f"functional_ora_volcano_{safe_source}.png")
        plt.savefig(fig_path)
        plt.close()
        print(f"    Saved {fig_path}")


# ================================================================
# Main
# ================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Functional ORA for mimetope-enriched proteins using gProfiler"
    )
    parser.add_argument("--enrichment", required=True,
                        help="protein_enrichment_collapsed.tsv")
    parser.add_argument("--proteome", required=True,
                        help="Human proteome FASTA (for Swiss-Prot background)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="FDR threshold (default: 0.05)")
    parser.add_argument("--force-requery", action="store_true",
                        help="Force re-query of gProfiler API (ignore cache)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --- Load MEPs ---
    print("\n" + "=" * 60)
    print("Loading data")
    print("=" * 60)

    print(f"  Loading enrichment results: {args.enrichment}")
    enrich_df = pd.read_csv(args.enrichment, sep="\t")
    print(f"    {len(enrich_df):,} proteins total")

    mep = enrich_df[enrich_df["significant"] == True]
    mep_acs = sorted(mep["uniprot_ac"].unique())
    print(f"    MEPs (significant): {len(mep_acs):,}")

    if len(mep_acs) == 0:
        print("  ⚠️ No MEPs found. Nothing to test.")
        return

    # --- Swiss-Prot background ---
    sp_acs = get_swissprot_acs(args.proteome)

    # Ensure MEPs are in background
    mep_in_bg = [ac for ac in mep_acs if ac in set(sp_acs)]
    print(f"    MEPs in Swiss-Prot background: {len(mep_in_bg):,} / {len(mep_acs):,}")

    if len(mep_in_bg) < len(mep_acs):
        missing = set(mep_acs) - set(sp_acs)
        print(f"    ⚠️ {len(missing)} MEPs not in Swiss-Prot (TrEMBL entries?)")
        print(f"       These will still be queried but may not map in gProfiler")

    # --- Query gProfiler ---
    print(f"\n{'=' * 60}")
    print("Querying gProfiler")
    print("=" * 60)

    cache_path = os.path.join(args.outdir, "gprofiler_cache.json")
    if args.force_requery and os.path.exists(cache_path):
        os.remove(cache_path)
        print("  Removed cached results (--force-requery)")

    data = query_gprofiler(mep_acs, sp_acs, cache_path)

    # --- Parse results ---
    print(f"\n{'=' * 60}")
    print("Parsing results")
    print("=" * 60)

    results_df = parse_gprofiler_results(data)

    if results_df.empty:
        print("  No results returned.")
        return

    # --- Save ---
    print(f"\n{'=' * 60}")
    print("Saving results")
    print("=" * 60)

    all_path = os.path.join(args.outdir, "functional_ora_results.tsv")
    results_df.sort_values(["source", "p_value"]).to_csv(
        all_path, sep="\t", index=False
    )
    print(f"  All results: {all_path}")

    sig_df = results_df[results_df["p_value"] < args.fdr]
    sig_path = os.path.join(args.outdir, "functional_ora_significant.tsv")
    sig_df.sort_values(["source", "p_value"]).to_csv(
        sig_path, sep="\t", index=False
    )
    print(f"  Significant (FDR < {args.fdr}): {sig_path}")
    print(f"    {len(sig_df)} significant terms")

    # Summary
    print(f"\n{'=' * 60}")
    print("Summary")
    print("=" * 60)
    print(f"  MEPs queried: {len(mep_acs)}")
    print(f"  Background: {len(sp_acs):,} Swiss-Prot proteins")
    print(f"  Total terms tested: {len(results_df):,}")
    print(f"  Significant terms (FDR < {args.fdr}): {len(sig_df):,}")

    for source in ["GO:BP", "GO:MF", "GO:CC", "KEGG"]:
        sub_sig = sig_df[sig_df["source"] == source]
        if not sub_sig.empty:
            print(f"\n  {SOURCE_NAMES.get(source, source)}: {len(sub_sig)} terms")
            for _, row in sub_sig.head(5).iterrows():
                print(f"    {row['term_name']}: FE={row['fold_enrichment']:.1f}, "
                      f"n={int(row['intersection_size'])}, "
                      f"FDR={row['p_value']:.2e}")

    # --- Plots ---
    print(f"\n{'=' * 60}")
    print("Generating plots")
    print("=" * 60)
    plot_functional_ora(results_df, args.outdir, args.fdr)

    print(f"\n✅ Done.")


if __name__ == "__main__":
    main()