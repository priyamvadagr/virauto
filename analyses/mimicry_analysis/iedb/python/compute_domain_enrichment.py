#!/usr/bin/env python3
"""
compute_domain_enrichment.py

Positional domain overrepresentation analysis for mimetopes.

Question: Are mimetopes enriched within specific protein domains,
relative to what we'd expect if they were randomly distributed
across all possible k-mer windows in the human proteome?

Null model (Null 3 — k-mer window counting):
    For each domain D and peptide length k:
        bg_windows_D = total k-mer windows in the proteome that overlap D
                       by majority (>50% of residues inside D)
        bg_windows_total = total k-mer windows in the proteome
        fg_windows_D = mimetopes of length k that overlap D by majority
        fg_windows_total = total mimetopes of length k

    Fisher's exact test on the 2×2 table:
                        overlaps D      does not overlap D
        mimetopes       fg_D            fg_total - fg_D
        background      bg_D            bg_total - bg_D

    BH correction across all domains tested.

Overlap criterion: majority overlap (>50% of the k-mer's residues
fall within the domain boundaries).

Input:
    - pair_id_mapping.csv.gz (with hu_uniprot, sstart, send, human_mimic_sequence)
    - iedb_mhci_strong_mimicry.csv.gz (for pair_ids in the strong set)
    - protein2ipr.dat.gz (InterPro domain annotations with coordinates)
    - uniprot_human_all.fasta (human proteome for background window counts)

Output:
    - domain_ora_positional_k{k}.tsv per peptide length
    - domain_ora_positional_summary.tsv (all k combined)
    - Figures: domain ORA bar chart + volcano per k

Usage:
    python compute_domain_enrichment.py \
        --pair-map pair_id_mapping.csv.gz \
        --strong-mimicry iedb_mhci_strong_mimicry.csv.gz \
        --protein2ipr protein2ipr.dat.gz \
        --proteome uniprot_human_all.fasta \
        --outdir results/domain_enrichment \
        --fdr 0.05
"""

import argparse
import gzip
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
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
    """
    Parse UniProt FASTA → dict of AC → sequence length.
    """
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
    Stream protein2ipr.dat.gz and extract domain annotations
    for human proteins. All InterPro member databases included.

    Returns:
        ac2domains: dict AC → list of (domain_name, source_db, start, end)
        domain_names: dict (domain_name, source_db) → set of ACs
    """
    print(f"  Loading InterPro annotations from {filepath} ...")
    ac2domains = defaultdict(list)
    domain_info = defaultdict(set)  # (name, source_db) → set of ACs
    n_kept = 0

    with gzip.open(filepath, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            ac = parts[0]
            if ac not in keep_acs:
                continue
            ipr_name = parts[2]
            src_db = parts[3]
            try:
                start = int(parts[5])
                end = int(parts[6])
            except ValueError:
                continue
            ac2domains[ac].append((ipr_name, src_db, start, end))
            domain_info[(ipr_name, src_db)].add(ac)
            n_kept += 1

    print(f"    {n_kept:,} annotation lines kept")
    print(f"    {len(ac2domains):,} proteins with domains")
    print(f"    {len(domain_info):,} unique (domain, source_db) entries")
    return ac2domains, domain_info


def kmer_overlaps_domain(kmer_start, kmer_end, dom_start, dom_end):
    """
    Check if a k-mer has majority overlap (>50% of its residues)
    with a domain.

    All coordinates are 1-based inclusive.
    """
    overlap_start = max(kmer_start, dom_start)
    overlap_end = min(kmer_end, dom_end)
    if overlap_start > overlap_end:
        return False
    overlap_len = overlap_end - overlap_start + 1
    kmer_len = kmer_end - kmer_start + 1
    return overlap_len > kmer_len / 2


# ================================================================
# Background: count k-mer windows overlapping each domain
# ================================================================

def count_background_windows(ac2len, ac2domains, k_values):
    """
    For each domain and each k, count:
        - total k-mer windows across the proteome
        - k-mer windows that overlap the domain by majority

    Returns:
        bg_total: dict k → total windows in proteome
        bg_per_domain: dict (domain_name, source_db) → dict k → window count
    """
    print(f"\n  Computing background k-mer windows ...")

    bg_total = {k: 0 for k in k_values}
    # (domain_name, source_db) → {k: count}
    bg_per_domain = defaultdict(lambda: {k: 0 for k in k_values})

    n_proteins = len(ac2len)
    for idx, (ac, seq_len) in enumerate(ac2len.items()):
        if (idx + 1) % 5000 == 0:
            print(f"    Processed {idx + 1:,} / {n_proteins:,} proteins ...")

        domains = ac2domains.get(ac, [])

        for k in k_values:
            n_windows = seq_len - k + 1
            if n_windows <= 0:
                continue
            bg_total[k] += n_windows

            # For each window, check overlap with each domain
            # Optimization: only iterate windows if protein has domains
            if not domains:
                continue

            for w_start in range(1, n_windows + 1):
                w_end = w_start + k - 1
                for (dom_name, src_db, dom_start, dom_end) in domains:
                    if kmer_overlaps_domain(w_start, w_end, dom_start, dom_end):
                        bg_per_domain[(dom_name, src_db)][k] += 1
                        # Don't break — a window can overlap multiple domains
                        # and we want to count it for each

    for k in k_values:
        print(f"    k={k}: {bg_total[k]:,} total windows, "
              f"{sum(1 for d in bg_per_domain if bg_per_domain[d][k] > 0)} "
              f"domains with ≥1 window")

    return bg_total, bg_per_domain


def count_background_windows_fast(ac2len, ac2domains, k_values):
    """
    Optimized background counting. Instead of iterating every k-mer
    window, for each (protein, domain) pair compute the number of
    k-mer windows with majority overlap analytically.

    For a domain at [dom_start, dom_end] in a protein of length L,
    a k-mer starting at position s (1-based) spans [s, s+k-1].
    Majority overlap means overlap > k/2, i.e., overlap >= floor(k/2)+1.

    The overlap between [s, s+k-1] and [dom_start, dom_end] is:
        max(0, min(s+k-1, dom_end) - max(s, dom_start) + 1)

    We need this >= threshold where threshold = k//2 + 1.

    This can be solved by finding the range of valid s values.
    """
    print(f"\n  Computing background k-mer windows (fast) ...")

    bg_total = {}
    for k in k_values:
        bg_total[k] = sum(max(0, slen - k + 1) for slen in ac2len.values())

    bg_per_domain = defaultdict(lambda: {k: 0 for k in k_values})

    n_proteins_with_domains = len(ac2domains)
    for idx, (ac, domains) in enumerate(ac2domains.items()):
        if (idx + 1) % 5000 == 0:
            print(f"    Processed {idx + 1:,} / {n_proteins_with_domains:,} "
                  f"annotated proteins ...")

        seq_len = ac2len.get(ac, 0)
        if seq_len == 0:
            continue

        for (dom_name, src_db, dom_start, dom_end) in domains:
            dom_len = dom_end - dom_start + 1

            for k in k_values:
                n_windows = seq_len - k + 1
                if n_windows <= 0:
                    continue

                threshold = k // 2 + 1  # minimum overlap for majority

                # Find range of start positions s where overlap >= threshold
                # overlap = min(s+k-1, dom_end) - max(s, dom_start) + 1
                #
                # Case analysis on s:
                # When s <= dom_start and s+k-1 >= dom_end (window contains domain):
                #   overlap = dom_len → valid if dom_len >= threshold
                # When s >= dom_start (window starts inside or after domain):
                #   overlap = min(s+k-1, dom_end) - s + 1
                # When s < dom_start (window starts before domain):
                #   overlap = min(s+k-1, dom_end) - dom_start + 1
                #
                # Easier: enumerate valid s range directly.
                # s_min: earliest start where overlap >= threshold
                # s_max: latest start where overlap >= threshold
                #
                # The overlap as a function of s is a trapezoid.
                # It increases as s approaches the domain, peaks when
                # the window is inside/covering the domain, then decreases.
                #
                # overlap >= threshold requires:
                #   s <= dom_end - threshold + 1   (window must reach into domain)
                #   s + k - 1 >= dom_start + threshold - 1  (window must extend enough)
                #
                # So: s >= dom_start + threshold - k
                #     s <= dom_end - threshold + 1

                s_min = max(1, dom_start + threshold - k)
                s_max = min(n_windows, dom_end - threshold + 1)

                if s_min <= s_max:
                    bg_per_domain[(dom_name, src_db)][k] += (s_max - s_min + 1)

    for k in k_values:
        n_domains_with_windows = sum(
            1 for d in bg_per_domain if bg_per_domain[d][k] > 0
        )
        print(f"    k={k}: {bg_total[k]:,} total windows, "
              f"{n_domains_with_windows:,} domains with ≥1 overlapping window")

    return bg_total, bg_per_domain


# ================================================================
# Foreground: count mimetopes overlapping each domain
# ================================================================

def count_foreground_overlaps(mimetopes_df, ac2domains, k_values):
    """
    For each mimetope, check which domains it overlaps by majority.

    mimetopes_df must have: hu_uniprot, sstart, send, pep_len

    Returns:
        fg_total: dict k → total mimetopes of that length
        fg_per_domain: dict (domain_name, source_db) → dict k → count
        mimetope_domains: list of dicts for per-mimetope annotation
    """
    print(f"\n  Counting foreground mimetope-domain overlaps ...")

    fg_total = {k: 0 for k in k_values}
    fg_per_domain = defaultdict(lambda: {k: 0 for k in k_values})
    mimetope_domain_records = []

    for _, row in mimetopes_df.iterrows():
        ac = row["hu_uniprot"]
        s_start = int(row["sstart"])
        s_end = int(row["send"])
        k = int(row["pep_len"])

        if k not in fg_total:
            continue

        fg_total[k] += 1

        domains = ac2domains.get(ac, [])
        overlapping = []
        for (dom_name, src_db, dom_start, dom_end) in domains:
            if kmer_overlaps_domain(s_start, s_end, dom_start, dom_end):
                fg_per_domain[(dom_name, src_db)][k] += 1
                overlapping.append(f"{dom_name} ({src_db})")

        mimetope_domain_records.append({
            "hu_uniprot": ac,
            "sstart": s_start,
            "send": s_end,
            "pep_len": k,
            "n_domains_overlapped": len(overlapping),
            "domains": "; ".join(overlapping) if overlapping else "none",
        })

    for k in k_values:
        n_domains_hit = sum(
            1 for d in fg_per_domain if fg_per_domain[d][k] > 0
        )
        print(f"    k={k}: {fg_total[k]:,} mimetopes, "
              f"{n_domains_hit:,} domains hit")

    return fg_total, fg_per_domain, mimetope_domain_records


# ================================================================
# Statistical testing
# ================================================================

def run_domain_ora(fg_total, fg_per_domain, bg_total, bg_per_domain,
                   k_values, fdr_threshold, min_fg=2):
    """
    Fisher's exact test per domain per k.

    2×2 table:
                    overlaps D      does not overlap D
    mimetopes       a               b
    background      c               d

    where:
        a = fg_per_domain[D][k]
        b = fg_total[k] - a
        c = bg_per_domain[D][k]
        d = bg_total[k] - c

    Returns DataFrame with test results.
    """
    print(f"\n  Running Fisher's exact tests ...")

    results = []
    for k in k_values:
        fg_t = fg_total[k]
        bg_t = bg_total[k]
        if fg_t == 0:
            print(f"    k={k}: no mimetopes, skipping")
            continue

        # Collect all domains that have either fg or bg windows for this k
        all_domains = set()
        for d in fg_per_domain:
            if fg_per_domain[d][k] > 0:
                all_domains.add(d)

        pvals = []
        test_results = []

        for (dom_name, src_db) in all_domains:
            a = fg_per_domain[(dom_name, src_db)][k]
            if a < min_fg:
                continue
            b = fg_t - a
            c = bg_per_domain[(dom_name, src_db)][k]
            d = bg_t - c

            # Guard against negative values
            c = max(c, 0)
            d = max(d, 0)

            _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

            fg_frac = a / fg_t if fg_t > 0 else 0
            bg_frac = c / bg_t if bg_t > 0 else 0
            fe = fg_frac / bg_frac if bg_frac > 0 else np.inf

            test_results.append({
                "domain": dom_name,
                "source_db": src_db,
                "k": k,
                "fg_mimetopes": a,
                "fg_total": fg_t,
                "bg_windows": c,
                "bg_total": bg_t,
                "fg_fraction": round(fg_frac, 6),
                "bg_fraction": round(bg_frac, 6),
                "fold_enrichment": round(fe, 4),
                "pval": pval,
            })
            pvals.append(pval)

        if not test_results:
            print(f"    k={k}: no domains with >= {min_fg} mimetopes")
            continue

        # BH correction
        _, qvals, _, _ = multipletests(pvals, alpha=fdr_threshold, method="fdr_bh")
        for i, res in enumerate(test_results):
            res["qval"] = qvals[i]
            res["significant"] = qvals[i] < fdr_threshold

        results.extend(test_results)

        n_sig = sum(1 for r in test_results if r["significant"])
        print(f"    k={k}: {len(test_results)} domains tested, "
              f"{n_sig} significant (FDR < {fdr_threshold})")

    return pd.DataFrame(results)


# ================================================================
# Plotting
# ================================================================

def plot_domain_ora(results_df, outdir, fdr_threshold):
    """
    Generate per-k bar chart + volcano plots.
    """
    fig_dir = os.path.join(outdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    k_values = sorted(results_df["k"].unique())

    for k in k_values:
        sub = results_df[results_df["k"] == k].copy()
        if sub.empty:
            continue

        sig = sub[sub["significant"] & (sub["fold_enrichment"] > 1)].copy()
        insig = sub[~(sub["significant"] & (sub["fold_enrichment"] > 1))]

        fig, axes = plt.subplots(1, 2, figsize=(17, 7))

        # Bar chart of top significant domains
        ax = axes[0]
        top = sig.nlargest(20, "fold_enrichment")
        if top.empty:
            ax.text(0.5, 0.5, f"No significant domains\nat FDR < {fdr_threshold}",
                    ha="center", va="center", transform=ax.transAxes, fontsize=12)
            ax.set_title(f"k={k}: Top enriched domains")
        else:
            bar_c = ["#D32F2F" if q < 0.01 else "#FF9800"
                      for q in top["qval"]]
            ax.barh(range(len(top)), top["fold_enrichment"].values[::-1],
                    color=bar_c[::-1], edgecolor="white")
            ax.axvline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
            ax.set_yticks(range(len(top)))
            labels = [f"{r['domain']} ({r['source_db']})"
                      for _, r in top.iloc[::-1].iterrows()]
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_xlabel("Fold enrichment (mimetopes vs proteome k-mer windows)")
            ax.set_title(f"k={k}: Top 20 enriched domains (FDR < {fdr_threshold})")

            for i, (_, row) in enumerate(top.iloc[::-1].iterrows()):
                ax.text(row["fold_enrichment"] + 0.02, i,
                        f"n={int(row['fg_mimetopes'])}, FDR={row['qval']:.2e}",
                        va="center", fontsize=6.5)

            legend_elements = [
                Patch(facecolor="#D32F2F", label="FDR < 0.01"),
                Patch(facecolor="#FF9800", label="FDR 0.01–0.05"),
            ]
            ax.legend(handles=legend_elements, fontsize=8)

        # Volcano
        ax = axes[1]
        ax.scatter(insig["fold_enrichment"],
                   -np.log10(insig["qval"].clip(lower=1e-300)),
                   s=15, alpha=0.3, color="#BDBDBD")
        if not sig.empty:
            ax.scatter(sig["fold_enrichment"],
                       -np.log10(sig["qval"].clip(lower=1e-300)),
                       s=30, alpha=0.8, color="#D32F2F",
                       label=f"FDR < {fdr_threshold}")

            # Label top 10
            sig_c = sig.copy()
            sig_c["score"] = sig_c["fold_enrichment"] * (
                -np.log10(sig_c["qval"].clip(lower=1e-300)))
            for _, row in sig_c.nlargest(10, "score").iterrows():
                ax.annotate(
                    f"{row['domain']}",
                    (row["fold_enrichment"],
                     -np.log10(max(row["qval"], 1e-300))),
                    fontsize=7, xytext=(4, 2), textcoords="offset points")

        ax.axvline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.axhline(-np.log10(fdr_threshold), color="gray", linestyle=":",
                   linewidth=0.8, label=f"FDR = {fdr_threshold}")
        ax.set_xlabel("Fold enrichment")
        ax.set_ylabel("-log10(q-value)")
        ax.set_title(f"k={k}: Domain ORA volcano\n"
                     f"(positional overlap, Fisher's exact, BH FDR)")
        ax.legend(fontsize=8)

        plt.tight_layout()
        fig_path = os.path.join(fig_dir, f"domain_ora_positional_k{k}.png")
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
                        help="protein2ipr.dat.gz from InterPro")
    parser.add_argument("--proteome", required=True,
                        help="Human proteome FASTA (UniProt)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="FDR threshold (default: 0.05)")
    parser.add_argument("--min-fg", type=int, default=2,
                        help="Minimum foreground mimetopes to test a domain "
                             "(default: 2)")
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

    # Compute peptide length from human_mimic_sequence
    mimetopes["pep_len"] = mimetopes["human_mimic_sequence"].str.len()

    # Deduplicate: unique (hu_uniprot, sstart, send) per mimetope
    mimetopes = mimetopes.drop_duplicates(
        subset=["hu_uniprot", "sstart", "send"]
    )
    print(f"    {len(mimetopes):,} unique mimetope positions after dedup")

    # Detect k values
    k_values = sorted(mimetopes["pep_len"].unique())
    print(f"    Peptide lengths: {k_values}")

    # --- Background windows ---
    print("\n" + "=" * 60)
    print("Background computation")
    print("=" * 60)
    bg_total, bg_per_domain = count_background_windows_fast(
        ac2len, ac2domains, k_values
    )

    # --- Foreground overlaps ---
    print("\n" + "=" * 60)
    print("Foreground computation")
    print("=" * 60)
    fg_total, fg_per_domain, mimetope_records = count_foreground_overlaps(
        mimetopes, ac2domains, k_values
    )

    # Save per-mimetope domain annotations
    mimetope_annot_path = os.path.join(args.outdir, "mimetope_domain_annotations.tsv")
    pd.DataFrame(mimetope_records).to_csv(mimetope_annot_path, sep="\t", index=False)
    print(f"  Mimetope annotations: {mimetope_annot_path}")

    # --- Fisher's exact tests ---
    print("\n" + "=" * 60)
    print("Statistical testing")
    print("=" * 60)
    results_df = run_domain_ora(
        fg_total, fg_per_domain, bg_total, bg_per_domain,
        k_values, args.fdr, min_fg=args.min_fg
    )

    if results_df.empty:
        print("\n  No domains tested. Check input data.")
        return

    # Save results
    results_path = os.path.join(args.outdir, "domain_ora_positional_all.tsv")
    results_df.sort_values(["k", "qval"]).to_csv(
        results_path, sep="\t", index=False
    )
    print(f"\n  All results: {results_path}")

    for k in sorted(results_df["k"].unique()):
        k_path = os.path.join(args.outdir, f"domain_ora_positional_k{k}.tsv")
        results_df[results_df["k"] == k].sort_values("qval").to_csv(
            k_path, sep="\t", index=False
        )

    # Summary
    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")
    for k in sorted(results_df["k"].unique()):
        sub = results_df[results_df["k"] == k]
        n_sig = sub["significant"].sum()
        print(f"  k={k}: {len(sub)} domains tested, {n_sig} significant "
              f"(FDR < {args.fdr})")
        if n_sig > 0:
            top3 = sub[sub["significant"]].nlargest(3, "fold_enrichment")
            for _, row in top3.iterrows():
                print(f"    {row['domain']} ({row['source_db']}): "
                      f"FE={row['fold_enrichment']:.1f}, "
                      f"n={int(row['fg_mimetopes'])}, "
                      f"q={row['qval']:.2e}")

    # --- Plots ---
    print(f"\n{'=' * 60}")
    print("Generating plots")
    print(f"{'=' * 60}")
    plot_domain_ora(results_df, args.outdir, args.fdr)

    print(f"\n✅ Done.")


if __name__ == "__main__":
    main()