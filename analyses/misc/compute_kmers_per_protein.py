from Bio import SeqIO
import pandas as pd

# -----------------------------
# Input / output paths
# -----------------------------
fasta_file = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_all.fasta"
remove_file = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/proteins_to_remove_from_UniProtKB.txt"
output_file = "/ix/djishnu/Priyamvada/virauto/data/refs/uniprot/uniprot_human_kmer_counts_filtered.tsv"

k_values = [8, 9, 10, 11]

# -----------------------------
# Read proteins to remove
# -----------------------------
with open(remove_file, "r") as f:
    proteins_to_remove = {line.strip() for line in f if line.strip()}

print(f"Loaded {len(proteins_to_remove)} proteins to remove")

results = []
n_total = 0
n_removed = 0
n_kept = 0

for record in SeqIO.parse(fasta_file, "fasta"):
    n_total += 1

    # UniProt FASTA headers look like:
    # tr|A0A024R1X5|A0A024R1X5_HUMAN
    parts = record.id.split("|")
    if len(parts) >= 3:
        protein_id = parts[1]      # accession
        entry_name = parts[2]      # entry name
    else:
        protein_id = record.id
        entry_name = record.id

    if protein_id in proteins_to_remove:
        n_removed += 1
        continue

    seq = str(record.seq).replace("*", "")
    seq_len = len(seq)

    row = {
        "protein_id": protein_id,
        "entry_name": entry_name,
        "description": record.description,
        "length": seq_len
    }

    for k in k_values:
        row[f"{k}mer_count"] = max(seq_len - k + 1, 0)

    results.append(row)
    n_kept += 1

df = pd.DataFrame(results)
df.to_csv(output_file, sep="\t", index=False)

print("\nFirst few rows:")
print(df.head().to_string(index=False))

print(f"\nSaved filtered k-mer counts to: {output_file}")
print(f"Total proteins in FASTA: {n_total}")
print(f"Removed proteins: {n_removed}")
print(f"Remaining proteins: {n_kept}")

print("\nTotal k-mers across remaining proteins:")
for k in k_values:
    print(f"{k}-mers: {df[f'{k}mer_count'].sum():,}")