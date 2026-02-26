#%%
import Bio
import pandas as pd
import os
pair_summary = pd.read_csv("/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/virscan/9_mers/Class_I/HLA-A/summaries/pair_binding_summary_HLA_A.csv")
import matplotlib.pyplot as plt
import seaborn as sns
#%%
sns.histplot(pair_summary["Human_dominant_fraction"], 
 bins=100, color="#27ae60", alpha=0.7, label="Human_dominant_fraction")
#%%
def classify_pair(row, thr=0.05):
    if row["Human_dominant_fraction"] >= thr:
        return "Human_dominant"
    elif row["Viral_dominant_fraction"] >= thr:
        return "Viral_dominant"
    else:
        return "Mixed"

pair_summary["Pair_Category"] = pair_summary.apply(classify_pair, axis=1)
pair_summary_filtered = pair_summary[
    (pair_summary["Pair_Category"] == "Human_dominant") 
]
plot_filename = "volcano_peptide_pair_HLA_A.png"
sns.scatterplot(
    data=pair_summary,
    x="Human_dominant_fraction",
    y="Viral_dominant_fraction",
    hue="Pair_Category",
    palette={"Human_dominant": "#27ae60", "Viral_dominant": "#c0392b", "Mixed": "#95a5a6"},
    alpha=0.8
)
plt.title("Mimicry Volcano Plot at Peptide-Pair Level")
plot_path = os.path.join('/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/virscan/9_mers/Class_I/HLA-A/plots', f"{plot_filename}")

plt.savefig(plot_path, dpi=300)
plt.show()
#%%
pair_summary_filtered


# %%
pair_summary["bias_score"] = (
    pair_summary["Human_dominant_fraction"] -
    pair_summary["Viral_dominant_fraction"]
)
sns.histplot(pair_summary["bias_score"], bins=50, kde=True, color="#2c3e50")
plt.xlabel("Binding Bias (Human − Viral Fraction)")
plt.ylabel("Number of Peptide Pairs")
plt.title("Distribution of Allele-Level Binding Bias")
plt.show()
# %%

group_cols = ["uniprot_viral", "uniprot_human"]

# First: aggregate at the pair level (optional if already done)
pair_summary["is_human_dominant"] = pair_summary["Human_dominant_fraction"] >= 0.01

# Now group by viral/human protein pairs and count peptides meeting the threshold
summary = (
    pair_summary.groupby(["uniprot_viral", "uniprot_human"])
    .agg(
        n_peptides=("pair_id", "count"),
        n_human_dominant_peptides=("is_human_dominant", "sum"),
    )
    .reset_index()
)

# Compute fraction of peptides that are human-dominant
summary["fraction_human_dominant_peptides"] = (
    summary["n_human_dominant_peptides"] / summary["n_peptides"]
)


sns.histplot(summary["fraction_human_dominant_peptides"], bins=50, color="#8e44ad")
plt.xlabel("Fraction of Human-Dominant Peptides")
plt.ylabel("Number of Protein Pairs")
plt.title("Distribution of Human-Dominant Peptide Fractions at Protein-Pair Level")
plt.show()
# %%
summary[summary["fraction_human_dominant_peptides"] > 0]['uniprot_human'].value_counts()

# %%
summary[summary["fraction_human_dominant_peptides"] > 0]['uniprot_human'].value_counts().to_csv(
    '/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/virscan/9_mers/Class_I/HLA-A/summaries/human_protein_counts.csv',
    header=['count']
)
# %%
human_dom_peptides = pair_summary[
    pair_summary["Human_dominant_fraction"] >= 0.01
].copy()

print(f"[INFO] Selected {len(human_dom_peptides)} human-dominant peptide pairs")

unique_peptides = human_dom_peptides["peptide_human"].unique()
print(f"[INFO] {len(unique_peptides)} unique human peptides selected")
from itertools import combinations
import numpy as np
import pandas as pd

def jaccard_similarity(seq1, seq2, k=3):
    k1 = {seq1[i:i+k] for i in range(len(seq1)-k+1)}
    k2 = {seq2[i:i+k] for i in range(len(seq2)-k+1)}
    return len(k1 & k2) / len(k1 | k2) if len(k1 | k2) > 0 else 0

# Compute average pairwise similarity
pairs = list(combinations(unique_peptides, 2))
similarities = [jaccard_similarity(a, b, k=3) for a, b in pairs]

avg_similarity = np.mean(similarities)
print(f"[INFO] Average 5-mer Jaccard similarity among human-dominant peptides: {avg_similarity:.3f}")
import seaborn as sns
import matplotlib.pyplot as plt

sns.histplot(similarities, bins=40, color="#138a36")
plt.xlabel("Pairwise 3-mer Jaccard Similarity")
plt.ylabel("Number of Peptide Pairs")
plt.title("Sequence Similarity Among Human-Dominant Peptides")
plt.tight_layout()
plt.show()
# %%
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch

# Convert similarity list to condensed distance matrix
distances = 1 - np.array(similarities)
Z = sch.linkage(distances, method='average')

plt.figure(figsize=(8, 6))
sch.dendrogram(Z, no_labels=True, color_threshold=0.3)
plt.title("Clustering of Human-Dominant Peptides by Sequence Similarity")
plt.show()
# %%
from Bio import pairwise2
from Bio.Align import substitution_matrices

matrix = substitution_matrices.load('BLOSUM62') 

def blosum_score(seq1, seq2):
    alignments = pairwise2.align.globaldx(seq1, seq2, matrix, score_only=True)
    return alignments / max(len(seq1), len(seq2))

# Compute pairwise mean similarity
unique_peptides = human_dom_peptides["peptide_human"].unique()
pairs = list(combinations(unique_peptides, 2))  # sample if large
scores = [blosum_score(a, b) for a, b in pairs]
avg_score = np.mean(scores)
print(f"[INFO] Mean normalized BLOSUM62 similarity: {avg_score:.3f}")

# %%
# Convert similarity list to condensed distance matrix
distances = 1 - np.array(similarities)
Z = sch.linkage(distances, method='average')

plt.figure(figsize=(8, 6))
sch.dendrogram(Z, no_labels=True, color_threshold=0.3)
plt.title("Clustering of Human-Dominant Peptides by Sequence Similarity")
plt.show()
# %%
import logomaker
import pandas as pd

from collections import Counter

# Convert peptide list into position frequency matrix
positions = len(unique_peptides[0])
pfm = pd.DataFrame([{aa: sum(p[i] == aa for p in unique_peptides)
                     for aa in "ACDEFGHIKLMNPQRSTVWY"}
                    for i in range(positions)])
pfm = pfm.div(pfm.sum(axis=1), axis=0)  # normalize

logomaker.Logo(pfm, color_scheme='chemistry')
plt.title("Consensus Motif - Human-Dominant Peptides")
plt.show()
# %%
threshold = 1.5
df_high = similarities.copy()
df_high[df_high < threshold] = 0   # or np.nan if you prefer white cells
