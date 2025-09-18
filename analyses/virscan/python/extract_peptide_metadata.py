# Process virscan dataset to extract peptide metadata
import pandas as pd

# Load the full VirScan dataset (large file)
df = pd.read_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/Virscan_dataset.csv")

# Extract unique peptide identifiers with mapping info
peptide_metadata = df[['pep_id', 'UniProt_acc', 'start', 'end']].drop_duplicates()

# Save as smaller CSV for downstream use
peptide_metadata.to_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/peptide_metadata.csv", index=False)

print(f"Saved {len(peptide_metadata)} unique peptides to peptide_metadata.csv")
