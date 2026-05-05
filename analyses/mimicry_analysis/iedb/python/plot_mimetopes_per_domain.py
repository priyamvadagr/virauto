#!/usr/bin/env python3
"""
Quick diagnostic: distribution of mimetope counts per domain
from the domain ORA mimetope annotations.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ANNOT_FILE = "/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/iedb/mhc_i/swissprot/ORA_domain_enrichment_full_overlap/mimetope_domain_annotations.tsv"

# Load
df = pd.read_csv(ANNOT_FILE, sep="\t")
print(f"Total mimetopes: {len(df)}")
print(f"Mimetopes with ≥1 domain: {(df['n_domains'] > 0).sum()}")
print(f"Mimetopes with no domain: {(df['n_domains'] == 0).sum()}")

# Explode domains per mimetope
has_domain = df[df["domains"] != "none"].copy()
all_domains = []
for _, row in has_domain.iterrows():
    for d in row["domains"].split("; "):
        all_domains.append(d)

domain_counts = pd.Series(all_domains).value_counts()
print(f"\nUnique domains hit: {len(domain_counts)}")
print(f"\nMimetope count distribution per domain:")
print(f"  Min:    {domain_counts.min()}")
print(f"  Median: {domain_counts.median():.0f}")
print(f"  Mean:   {domain_counts.mean():.1f}")
print(f"  Max:    {domain_counts.max()}")

# Percentiles
for p in [25, 50, 75, 90, 95, 99]:
    val = np.percentile(domain_counts.values, p)
    n_above = (domain_counts >= val).sum()
    print(f"  {p}th percentile: {val:.0f} ({n_above} domains above)")

# Count how many domains at each cutoff
print(f"\nDomains remaining at different min_fg thresholds:")
for t in [2, 3, 5, 8, 10, 15, 20]:
    n = (domain_counts >= t).sum()
    print(f"  min_fg={t}: {n} domains ({n/len(domain_counts)*100:.1f}%)")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.hist(domain_counts.values, bins=50, color="#5C6BC0", edgecolor="white", alpha=0.8)
ax.set_xlabel("Mimetopes per domain")
ax.set_ylabel("Number of domains")
ax.set_title("Distribution of mimetope counts per domain")

ax = axes[1]
bins_log = np.logspace(0, np.log10(domain_counts.max() + 1), 40)
ax.hist(domain_counts.values, bins=bins_log, color="#5C6BC0", edgecolor="white", alpha=0.8)
ax.set_xscale("log")
ax.set_xlabel("Mimetopes per domain (log scale)")
ax.set_ylabel("Number of domains")
ax.set_title("Log scale")

# Mark thresholds
for t, color in [(5, "orange"), (10, "red")]:
    ax.axvline(t, color=color, linestyle="--", linewidth=1.5, label=f"min_fg={t}")
ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig("mimetope_per_domain_distribution.png")
plt.close()
print(f"\nSaved: mimetope_per_domain_distribution.png")