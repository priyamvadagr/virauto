import pandas as pd

pairs = pd.read_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast/iedb_mhci_pairs_4digit_hla.csv.gz")
short = pairs[pairs["viral_sequence"].str.len() < 8]
print(f"Total pairs: {len(pairs)}")
print(f"Pairs with viral peptide < 8 aa: {len(short)}")
if not short.empty:
    print(f"\nExamples:")
    print(short[["viral_sequence", "human_mimic_sequence", "qlen", "mhc_allele"]].head(10))
    print(f"\nViral sequence lengths:")
    print(short["viral_sequence"].str.len().value_counts().sort_index())