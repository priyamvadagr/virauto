import pandas as pd

pairs = pd.read_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla.csv.gz")

# How many unique viral sequences per structure_id?
per_sid = pairs.groupby("structure_id").agg(
    n_viral_seqs=("viral_sequence", "nunique"),
    n_human_seqs=("human_mimic_sequence", "nunique"),
    n_human_prots=("hu_prot_id", "nunique"),
    n_alleles=("mhc_allele", "nunique"),
    viral_seqs=("viral_sequence", lambda x: list(x.unique())),
).reset_index()

# Cases where one structure_id maps to multiple viral sequences
multi_viral = per_sid[per_sid["n_viral_seqs"] > 1]
print(f"Total structure_ids: {len(per_sid)}")
print(f"With multiple viral sequences: {len(multi_viral)}")
if not multi_viral.empty:
    print(f"\nExamples:")
    print(multi_viral[["structure_id", "n_viral_seqs", "n_human_seqs", "viral_seqs"]].head(10))