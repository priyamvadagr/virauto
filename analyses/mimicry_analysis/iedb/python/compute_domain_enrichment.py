#!/usr/bin/env python3
"""
compute_domain_enrichment.py

Positional domain overrepresentation analysis for mimetopes.

Question: Are mimetopes enriched within specific protein domains,
relative to what we'd expect if they were randomly distributed
across all possible k-mer windows in the human proteome?

Null model (k-mer window counting):
    For each domain D and peptide length k:
        bg_windows_D = k-mer windows in proteome overlapping D (majority)
        bg_windows_total = total k-mer windows in proteome
        fg_windows_D = mimetopes of length k overlapping D (majority)
        fg_windows_total = total mimetopes of length k

    Fisher's exact test per (domain, k).
    Results collapsed by InterPro ID: max fold enrichment, min q-value.

Overlap criterion: full containment (entire k-mer within domain boundaries).

Input:
    - pair_id_mapping.csv.gz (with hu_uniprot, sstart, send)
    - iedb_mhci_strong_mimicry.csv.gz
    - protein2ipr_human.tsv.gz (from filter_interpro_human.py)
    - uniprot_human_all.fasta

Output:
    - domain_ora_per_k.tsv (all per-k results)
    - domain_ora_collapsed.tsv (collapsed by InterPro ID)
    - Figures: dot plot + volcano

Usage:
    python compute_domain_enrichment.py \
        --pair-map pair_id_mapping.csv.gz \
        --strong-mimicry iedb_mhci_strong_mimicry.csv.gz \
        --protein2ipr protein2ipr_human.tsv.gz \
        --proteome uniprot_human_all.fasta \
        --outdir results/domain_enrichment \
        --fdr 0.05
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from collections import defaultdict

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
    ac2len = {}
    current_ac = None
    current_len = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_ac is not None:
                    ac2len[current_ac] = current_len
                parts = line[1:].split("|")
                current_ac = parts[1].strip() if len(parts) >= 3 else None
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_ac is not None:
        ac2len[current_ac] = current_len
    print(f"    {len(ac2len):,} proteins parsed")
    return ac2len


def load_interpro_domains(filepath, keep_acs):
    """
    Load human-only InterPro annotations from protein2ipr_human.tsv.gz.
    Groups by InterPro ID (ipr_id) to collapse redundant source DB entries.

    Returns:
        ac2domains: dict AC → list of (ipr_id, ipr_name, start, end)
        domain_info: dict ipr_id → {"name": str, "acs": set}
    """
    print(f"  Loading InterPro annotations from {filepath} ...")

    ipr = pd.read_csv(filepath, sep="\t", low_memory=False)
    print(f"    Total rows: {len(ipr):,}")

    ipr = ipr[ipr["uniprot_ac"].isin(keep_acs)]
    print(f"    Rows matching proteome: {len(ipr):,}")

    ac2domains = defaultdict(list)
    domain_info = {}

    for _, row in ipr.iterrows():
        ac = row["uniprot_ac"]
        ipr_id = row["ipr_id"]
        ipr_name = row["ipr_name"]
        start = int(row["start"])
        end = int(row["end"])

        ac2domains[ac].append((ipr_id, ipr_name, start, end))

        if ipr_id not in domain_info:
            domain_info[ipr_id] = {"name": ipr_name, "acs": set()}
        domain_info[ipr_id]["acs"].add(ac)

    print(f"    {len(ac2domains):,} proteins with domains")
    print(f"    {len(domain_info):,} unique InterPro IDs")
    return ac2domains, domain_info


def kmer_overlaps_domain(kmer_start, kmer_end, dom_start, dom_end):
    """Check if the entire k-mer falls within the domain boundaries."""
    return kmer_start >= dom_start and kmer_end <= dom_end


# ================================================================
# Background: count k-mer windows overlapping each domain
# ================================================================

def count_background_windows(ac2len, ac2domains, k_values):
    """
    For each InterPro ID and each k, analytically count k-mer windows
    in the proteome that overlap the domain by majority.
    """
    print(f"\n  Computing background k-mer windows ...")

    bg_total = {}
    for k in k_values:
        bg_total[k] = sum(max(0, slen - k + 1) for slen in ac2len.values())

    # ipr_id → {k: count}
    bg_per_domain = defaultdict(lambda: {k: 0 for k in k_values})

    n_proteins = len(ac2domains)
    for idx, (ac, domains) in enumerate(ac2domains.items()):
        if (idx + 1) % 5000 == 0:
            print(f"    Processed {idx + 1:,} / {n_proteins:,} annotated proteins ...")

        seq_len = ac2len.get(ac, 0)
        if seq_len == 0:
            continue

        for (ipr_id, ipr_name, dom_start, dom_end) in domains:
            for k in k_values:
                n_windows = seq_len - k + 1
                if n_windows <= 0:
                    continue

                # Full containment: k-mer [s, s+k-1] must satisfy
                #   s >= dom_start  AND  s+k-1 <= dom_end
                # So: s >= dom_start  AND  s <= dom_end - k + 1
                s_min = max(1, dom_start)
                s_max = min(n_windows, dom_end - k + 1)

                if s_min <= s_max:
                    bg_per_domain[ipr_id][k] += (s_max - s_min + 1)

    for k in k_values:
        n_with = sum(1 for d in bg_per_domain if bg_per_domain[d][k] > 0)
        print(f"    k={k}: {bg_total[k]:,} total windows, "
              f"{n_with:,} InterPro IDs with ≥1 window")

    return bg_total, bg_per_domain


# ================================================================
# Foreground: count mimetopes overlapping each domain
# ================================================================

def count_foreground_overlaps(mimetopes_df, ac2domains, k_values):
    print(f"\n  Counting foreground mimetope-domain overlaps ...")

    fg_total = {k: 0 for k in k_values}
    # ipr_id → {k: count}
    fg_per_domain = defaultdict(lambda: {k: 0 for k in k_values})
    mimetope_records = []

    for _, row in mimetopes_df.iterrows():
        ac = row["hu_prot_id"]
        s_start = int(row["sstart"])
        s_end = int(row["send"])
        k = int(row["pep_len"])

        if k not in fg_total:
            continue

        fg_total[k] += 1

        domains = ac2domains.get(ac, [])
        overlapping = []
        for (ipr_id, ipr_name, dom_start, dom_end) in domains:
            if kmer_overlaps_domain(s_start, s_end, dom_start, dom_end):
                fg_per_domain[ipr_id][k] += 1
                overlapping.append(f"{ipr_id}:{ipr_name}")

        mimetope_records.append({
            "hu_prot_id": ac,
            "sstart": s_start,
            "send": s_end,
            "pep_len": k,
            "n_domains": len(overlapping),
            "domains": "; ".join(overlapping) if overlapping else "none",
        })

    for k in k_values:
        n_hit = sum(1 for d in fg_per_domain if fg_per_domain[d][k] > 0)
        print(f"    k={k}: {fg_total[k]:,} mimetopes, {n_hit:,} InterPro IDs hit")

    return fg_total, fg_per_domain, mimetope_records


# ================================================================
# Statistical testing
# ================================================================

def run_domain_ora(fg_total, fg_per_domain, bg_total, bg_per_domain,
                   domain_info, k_values, fdr_threshold, min_fg=2):
    """
    Fisher's exact test per (InterPro ID, k).
    Then collapse across k: max fold enrichment, min q-value.
    """
    print(f"\n  Running Fisher's exact tests ...")

    # --- Per-k tests ---
    all_per_k = []

    for k in k_values:
        fg_t = fg_total[k]
        bg_t = bg_total[k]
        if fg_t == 0:
            continue

        test_results = []
        pvals = []

        for ipr_id in fg_per_domain:
            a = fg_per_domain[ipr_id].get(k, 0)
            if a < min_fg:
                continue

            b = fg_t - a
            c = bg_per_domain[ipr_id].get(k, 0)
            d = bg_t - c
            c = max(c, 0)
            d = max(d, 0)

            _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

            fg_frac = a / fg_t
            bg_frac = c / bg_t if bg_t > 0 else 0
            fe = fg_frac / bg_frac if bg_frac > 0 else np.inf

            info = domain_info.get(ipr_id, {"name": ipr_id})

            test_results.append({
                "ipr_id": ipr_id,
                "domain_name": info["name"],
                "k": k,
                "fg_mimetopes": a,
                "fg_total": fg_t,
                "bg_windows": c,
                "bg_total": bg_t,
                "fold_enrichment": round(fe, 4),
                "pval": pval,
            })
            pvals.append(pval)

        if not test_results:
            continue

        _, qvals, _, _ = multipletests(pvals, alpha=fdr_threshold, method="bonferroni")
        for i, res in enumerate(test_results):
            res["qval"] = qvals[i]
            res["significant"] = qvals[i] < fdr_threshold

        all_per_k.extend(test_results)
        n_sig = sum(1 for r in test_results if r["significant"])
        print(f"    k={k}: {len(test_results)} domains tested, {n_sig} significant (Bonferroni)")

    if not all_per_k:
        return pd.DataFrame(), pd.DataFrame()

    per_k_df = pd.DataFrame(all_per_k)

    # --- Collapse across k per InterPro ID ---
    print(f"\n  Collapsing across k-mers per InterPro ID ...")

    collapsed = []
    for ipr_id, group in per_k_df.groupby("ipr_id"):
        # Max fold enrichment row
        best_fe_row = group.loc[group["fold_enrichment"].idxmax()]
        # Min q-value
        min_qval = group["qval"].min()
        # Total mimetopes across all k
        total_fg = group["fg_mimetopes"].sum()

        collapsed.append({
            "ipr_id": ipr_id,
            "domain_name": best_fe_row["domain_name"],
            "max_fold_enrichment": best_fe_row["fold_enrichment"],
            "best_k_fe": int(best_fe_row["k"]),
            "fg_at_best_k": int(best_fe_row["fg_mimetopes"]),
            "total_fg_all_k": int(total_fg),
            "min_qval": min_qval,
            "min_qval_k": int(group.loc[group["qval"].idxmin(), "k"]),
            "n_k_tested": len(group),
            "n_k_significant": int(group["significant"].sum()),
            "significant": min_qval < fdr_threshold,
        })

    collapsed_df = pd.DataFrame(collapsed).sort_values("min_qval")
    n_sig = collapsed_df["significant"].sum()
    print(f"    {len(collapsed_df)} InterPro IDs collapsed, {n_sig} significant")

    return per_k_df, collapsed_df


# ================================================================
# Label repulsion helper
# ================================================================

def _add_repelled_labels(ax, df, x_col, q_col, label_col, n_labels=10):
    """
    Add labels to volcano plot with adjustText repulsion.
    Selects top n_labels by score = log2(FE) × -log10(q).
    """
    from adjustText import adjust_text

    df = df.copy()
    df["_neg_log_q"] = -np.log10(df[q_col].clip(lower=1e-300))
    df["_score"] = df[x_col] * df["_neg_log_q"]
    top = df.nlargest(n_labels, "_score")

    if top.empty:
        return

    texts = []
    for _, row in top.iterrows():
        x = row[x_col]
        y = row["_neg_log_q"]
        name = row[label_col]
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


# ================================================================
# Plotting
# ================================================================

def plot_domain_ora(per_k_df, collapsed_df, outdir, fdr_threshold):
    fig_dir = os.path.join(outdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    # --- Per-k figures ---
    if not per_k_df.empty:
        k_values = sorted(per_k_df["k"].unique())
        for k in k_values:
            sub = per_k_df[per_k_df["k"] == k].copy()
            if sub.empty:
                continue

            sig_k = sub[sub["significant"]].copy()
            insig_k = sub[~sub["significant"]]

            fig, axes = plt.subplots(1, 2, figsize=(18, 8))

            # Dot plot
            ax = axes[0]
            top_k = sig_k.copy()

            if top_k.empty:
                ax.text(0.5, 0.5, f"No significant domains\n(k={k}, FDR < {fdr_threshold})",
                        ha="center", va="center", transform=ax.transAxes, fontsize=12)
                ax.set_title(f"k={k}: Top enriched domains")
            else:
                # Rank by score = log2(FE) × -log10(q)
                top_k = top_k.copy()
                top_k["_log2_fe"] = np.log2(top_k["fold_enrichment"].clip(lower=1e-10))
                top_k["_score"] = top_k["_log2_fe"] * (-np.log10(top_k["qval"].clip(lower=1e-300)))
                top_k = top_k.nlargest(25, "_score")
                neg_log_q = -np.log10(top_k["qval"].clip(lower=1e-300))
                norm = plt.Normalize(vmin=1.3, vmax=10)
                cmap = plt.cm.RdYlBu_r

                n_vals = top_k["fg_mimetopes"].values.astype(float)
                min_dot, max_dot = 50, 500
                n_max = max(n_vals.max(), 1)
                dot_sizes = min_dot + (n_vals / n_max) * (max_dot - min_dot)

                scatter = ax.scatter(
                    np.log2(top_k["fold_enrichment"].clip(lower=1e-10)),
                    range(len(top_k)),
                    s=dot_sizes, c=neg_log_q.values, cmap=cmap, norm=norm,
                    edgecolors="white", linewidths=0.5, zorder=3,
                )

                ax.set_yticks(range(len(top_k)))
                ax.set_yticklabels(top_k["domain_name"].tolist(), fontsize=7)
                ax.invert_yaxis()
                ax.set_xlabel("log2(fold enrichment)")
                ax.set_title(f"k={k}: Top 25 enriched domains (FDR < {fdr_threshold})")

                cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
                cbar.set_label("-log10(q-value)")

                legend_sizes = sorted(set([
                    max(1, int(n_max * 0.1)),
                    max(1, int(n_max * 0.5)),
                    int(n_max),
                ]))
                legend_dots = [
                    ax.scatter([], [], s=min_dot + (s / n_max) * (max_dot - min_dot),
                               c="gray", edgecolors="white", linewidths=0.5,
                               label=f"{s} mimetopes")
                    for s in legend_sizes
                ]
                ax.legend(handles=legend_dots, title="Mimetope count",
                          loc="lower right", fontsize=7, title_fontsize=8, frameon=True)

            # Volcano
            ax = axes[1]
            ax.scatter(
                np.log2(insig_k["fold_enrichment"].clip(lower=1e-10)),
                -np.log10(insig_k["qval"].clip(lower=1e-300)),
                s=15, alpha=0.3, color="#BDBDBD",
            )
            if not sig_k.empty:
                ax.scatter(
                    np.log2(sig_k["fold_enrichment"].clip(lower=1e-10)),
                    -np.log10(sig_k["qval"].clip(lower=1e-300)),
                    s=30, alpha=0.8, color="#D32F2F",
                    label=f"FDR < {fdr_threshold}",
                )
                sig_k_plot = sig_k.copy()
                sig_k_plot["_plot_x"] = np.log2(sig_k_plot["fold_enrichment"].clip(lower=1e-10))
                _add_repelled_labels(ax, sig_k_plot, "_plot_x", "qval",
                                     "domain_name", n_labels=10)

            ax.axvline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
            ax.axhline(-np.log10(fdr_threshold), color="gray", linestyle=":",
                       linewidth=0.8, label=f"FDR = {fdr_threshold}")
            ax.set_xlabel("log2(fold enrichment)")
            ax.set_ylabel("-log10(q-value)")
            ax.set_title(f"k={k}: Domain ORA volcano")
            ax.legend(fontsize=8)

            plt.tight_layout()
            fig_path = os.path.join(fig_dir, f"domain_ora_k{k}.png")
            plt.savefig(fig_path)
            plt.close()
            print(f"    Saved {fig_path}")

    # --- Collapsed figure ---
    sig = collapsed_df[collapsed_df["significant"]].copy()
    insig = collapsed_df[~collapsed_df["significant"]]

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # Dot plot
    ax = axes[0]
    if not sig.empty:
        sig = sig.copy()
        sig["_log2_fe"] = np.log2(sig["max_fold_enrichment"].clip(lower=1e-10))
        sig["_score"] = sig["_log2_fe"] * (-np.log10(sig["min_qval"].clip(lower=1e-300)))
        top = sig.nlargest(25, "_score")
    else:
        top = sig

    if top.empty:
        ax.text(0.5, 0.5, f"No significant domains\nat FDR < {fdr_threshold}",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_title("Domain ORA: top enriched domains")
    else:
        neg_log_q = -np.log10(top["min_qval"].clip(lower=1e-300))
        norm = plt.Normalize(vmin=1.3, vmax=10)
        cmap = plt.cm.RdYlBu_r

        n_vals = top["total_fg_all_k"].values.astype(float)
        min_dot, max_dot = 50, 500
        n_max = max(n_vals.max(), 1)
        dot_sizes = min_dot + (n_vals / n_max) * (max_dot - min_dot)

        scatter = ax.scatter(
            np.log2(top["max_fold_enrichment"].clip(lower=1e-10)),
            range(len(top)),
            s=dot_sizes, c=neg_log_q.values, cmap=cmap, norm=norm,
            edgecolors="white", linewidths=0.5, zorder=3,
        )

        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["domain_name"].tolist(), fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("log2(max fold enrichment)")
        ax.set_title(f"Top 25 enriched domains (FDR < {fdr_threshold})\n(collapsed by InterPro ID)")

        cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
        cbar.set_label("-log10(min q-value)")

        legend_sizes = sorted(set([
            max(1, int(n_max * 0.1)),
            max(1, int(n_max * 0.5)),
            int(n_max),
        ]))
        legend_dots = [
            ax.scatter([], [], s=min_dot + (s / n_max) * (max_dot - min_dot),
                       c="gray", edgecolors="white", linewidths=0.5,
                       label=f"{s} mimetopes")
            for s in legend_sizes
        ]
        ax.legend(handles=legend_dots, title="Mimetope count",
                  loc="lower right", fontsize=7, title_fontsize=8, frameon=True)

    # Volcano
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
        sig_plot = sig.copy()
        sig_plot["_plot_x"] = np.log2(sig_plot["max_fold_enrichment"].clip(lower=1e-10))
        _add_repelled_labels(ax, sig_plot, "_plot_x", "min_qval",
                             "domain_name", n_labels=10)

    ax.axvline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axhline(-np.log10(fdr_threshold), color="gray", linestyle=":",
               linewidth=0.8, label=f"FDR = {fdr_threshold}")
    ax.set_xlabel("log2(max fold enrichment)")
    ax.set_ylabel("-log10(min q-value)")
    ax.set_title("Domain ORA volcano\n(collapsed by InterPro ID, best across k)")
    ax.legend(fontsize=8)

    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "domain_ora_collapsed.png")
    plt.savefig(fig_path)
    plt.close()
    print(f"    Saved {fig_path}")


# ================================================================
# Main
# ================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Positional domain ORA for mimetopes"
    )
    parser.add_argument("--pair-map", required=True,
                        help="pair_id_mapping.csv.gz")
    parser.add_argument("--strong-mimicry", required=True,
                        help="iedb_mhci_strong_mimicry.csv.gz")
    parser.add_argument("--protein2ipr", required=True,
                        help="protein2ipr_human.tsv.gz (from filter_interpro_human.py)")
    parser.add_argument("--proteome", required=True,
                        help="Human proteome FASTA (UniProt)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="FDR threshold (default: 0.05)")
    parser.add_argument("--min-fg", type=int, default=10,
                        help="Minimum foreground mimetopes to test a domain "
                             "(default: 10)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --- Load proteome ---
    print("\n" + "=" * 60)
    print("Loading data")
    print("=" * 60)
    ac2len = parse_proteome_fasta(args.proteome)
    proteome_acs = set(ac2len.keys())

    # --- Load InterPro domains ---
    ac2domains, domain_info = load_interpro_domains(
        args.protein2ipr, keep_acs=proteome_acs
    )

    # --- Load mimetopes ---
    print(f"  Loading pair ID mapping: {args.pair_map}")
    pair_map = pd.read_csv(args.pair_map, low_memory=False)
    print(f"    {len(pair_map):,} rows")

    print(f"  Loading strong mimicry: {args.strong_mimicry}")
    strong = pd.read_csv(args.strong_mimicry, low_memory=False)
    strong_pair_ids = set(strong["pair_id"].unique())
    print(f"    {len(strong_pair_ids):,} unique strong pair_ids")

    # Filter pair_map to strong mimicry pairs
    mimetopes = pair_map[pair_map["pair_id"].isin(strong_pair_ids)].copy()
    print(f"    {len(mimetopes):,} mimetope rows (strong pairs with coordinates)")

    # Standardize human UniProt AC column name
    if "hu_prot_id" not in mimetopes.columns and "hu_uniprot" in mimetopes.columns:
        mimetopes = mimetopes.rename(columns={"hu_uniprot": "hu_prot_id"})

    mimetopes["pep_len"] = mimetopes["human_mimic_sequence"].str.len()

    # Deduplicate
    mimetopes = mimetopes.drop_duplicates(
        subset=["hu_prot_id", "sstart", "send"]
    )
    print(f"    {len(mimetopes):,} unique mimetope positions after dedup")

    k_values = sorted(mimetopes["pep_len"].unique())
    print(f"    Peptide lengths: {k_values}")

    # --- Background ---
    print("\n" + "=" * 60)
    print("Background computation")
    print("=" * 60)
    bg_total, bg_per_domain = count_background_windows(
        ac2len, ac2domains, k_values
    )

    # --- Foreground ---
    print("\n" + "=" * 60)
    print("Foreground computation")
    print("=" * 60)
    fg_total, fg_per_domain, mimetope_records = count_foreground_overlaps(
        mimetopes, ac2domains, k_values
    )

    # Save mimetope annotations
    annot_path = os.path.join(args.outdir, "mimetope_domain_annotations.tsv")
    pd.DataFrame(mimetope_records).to_csv(annot_path, sep="\t", index=False)
    print(f"  Mimetope annotations: {annot_path}")

    # --- Mimetope count distribution per domain ---
    fig_dir = os.path.join(args.outdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    print(f"\n{'=' * 60}")
    print("Mimetope distribution per domain")
    print("=" * 60)

    domain_mimetope_counts = {}
    for ipr_id in fg_per_domain:
        total = sum(fg_per_domain[ipr_id].get(k, 0) for k in k_values)
        if total > 0:
            name = domain_info.get(ipr_id, {"name": ipr_id})["name"]
            domain_mimetope_counts[f"{name} ({ipr_id})"] = total

    count_series = pd.Series(domain_mimetope_counts).sort_values(ascending=False)

    print(f"  Domains with ≥1 mimetope: {len(count_series)}")
    print(f"  Mimetope count distribution:")
    print(f"    Min: {count_series.min()}")
    print(f"    Median: {count_series.median():.0f}")
    print(f"    Mean: {count_series.mean():.1f}")
    print(f"    Max: {count_series.max()}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.hist(count_series.values, bins=50, color="#5C6BC0", edgecolor="white", alpha=0.8)
    ax.set_xlabel("Number of mimetopes per domain")
    ax.set_ylabel("Number of domains")
    ax.set_title("Distribution of mimetope counts per InterPro domain")
    ax.axvline(count_series.median(), color="red", linestyle="--", linewidth=1,
               label=f"Median = {count_series.median():.0f}")
    ax.legend(fontsize=10)

    ax = axes[1]
    bins_log = np.logspace(0, np.log10(count_series.max() + 1), 40)
    ax.hist(count_series.values, bins=bins_log, color="#5C6BC0", edgecolor="white", alpha=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Number of mimetopes per domain (log scale)")
    ax.set_ylabel("Number of domains")
    ax.set_title("Distribution of mimetope counts (log scale)")
    ax.axvline(args.min_fg, color="orange", linestyle="--", linewidth=1,
               label=f"min_fg = {args.min_fg}")
    ax.legend(fontsize=10)

    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "mimetope_count_distribution.png")
    plt.savefig(fig_path)
    plt.close()
    print(f"    Saved {fig_path}")

    # --- Tests ---
    print("\n" + "=" * 60)
    print("Statistical testing")
    print("=" * 60)
    per_k_df, collapsed_df = run_domain_ora(
        fg_total, fg_per_domain, bg_total, bg_per_domain,
        domain_info, k_values, args.fdr, min_fg=args.min_fg
    )

    if collapsed_df.empty:
        print("\n  No domains tested. Check input data.")
        return

    # Save results
    per_k_path = os.path.join(args.outdir, "domain_ora_per_k.tsv")
    per_k_df.sort_values(["k", "qval"]).to_csv(per_k_path, sep="\t", index=False)
    print(f"\n  Per-k results: {per_k_path}")

    collapsed_path = os.path.join(args.outdir, "domain_ora_collapsed.tsv")
    collapsed_df.to_csv(collapsed_path, sep="\t", index=False)
    print(f"  Collapsed results: {collapsed_path}")

    # Summary
    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")
    print(f"  InterPro IDs tested: {len(collapsed_df)}")
    print(f"  Significant (FDR < {args.fdr}): {collapsed_df['significant'].sum()}")
    top5 = collapsed_df[collapsed_df["significant"]].head(5)
    if not top5.empty:
        print(f"\n  Top 5 by q-value:")
        for _, row in top5.iterrows():
            print(f"    {row['ipr_id']} — {row['domain_name']}: "
                  f"FE={row['max_fold_enrichment']:.1f} (k={int(row['best_k_fe'])}), "
                  f"n={int(row['total_fg_all_k'])}, "
                  f"q={row['min_qval']:.2e}")

    # --- Plots ---
    print(f"\n{'=' * 60}")
    print("Generating plots")
    print(f"{'=' * 60}")
    plot_domain_ora(per_k_df, collapsed_df, args.outdir, args.fdr)

    print(f"\n✅ Done.")


if __name__ == "__main__":
    main()