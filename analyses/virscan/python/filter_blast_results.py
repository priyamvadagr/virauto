import pandas as pd

# Load BLAST output
cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
              "qlen","slen","qstart","qend","sstart","send","evalue","bitscore"]

df = pd.read_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/all_peptides_blast_out.tsv", sep="\t", names=cols)

# Compute coverage-normalized similarity
df["similarity"] = df["pident"] * (df["length"] / df["qlen"])

hits = df[(df["similarity"] >= 75) & (df["similarity"] <= 95)]
print(hits.head())

# Save to file
hits.to_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/similarity_filtered_blast_out.csv", index=False)