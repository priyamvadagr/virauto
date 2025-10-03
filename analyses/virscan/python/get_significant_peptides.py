#Get the peptide ids are significantly higher in MS cases vs controls 
import pandas as pd

# Read Excel files
post_comp_df = pd.read_excel("/ix/djishnu/Priyamvada/virauto/data/epitopes/post_onset_case_control_compare.xlsx")
pre_comp_df = pd.read_excel("/ix/djishnu/Priyamvada/virauto/data/epitopes/pre_onset_case_control_compare.xlsx")

# Function to annotate each dataset
def annotate_candidates(df, label):
    df = df.copy()
    df["candidate_" + label] = (
        (df["p.value"] < 0.1) & (df["High_in_cases"] == 1)
    )
    return df

post_annot = annotate_candidates(post_comp_df, "post")
pre_annot = annotate_candidates(pre_comp_df, "pre")

# Merge on peptide ID (outer join to keep everything)
merged = pd.merge(
    pre_annot,
    post_annot,
    on=["pep_id","Organism","Protein_name","UniProt_acc","pro_len","start","end"],
    how="outer",
    suffixes=("_pre","_post")
)

# Final candidate flag = true if candidate in pre OR post
merged["candidate"] = merged[["candidate_pre","candidate_post"]].any(axis=1)

# Keep only needed columns
final = merged[["pep_id","Organism","Protein_name","UniProt_acc","pro_len","start","end","candidate"]]

# Save results
final.to_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/annotatced_peptides_post_pre_combined.csv", index=False)

print(final.head(20))

# Plot number of peptides per organism 
import matplotlib.pyplot as plt

# Get count of unique pep_ids per species
subset_df = final[final['candidate'] == True]
pep_counts = subset_df.groupby('Organism')['pep_id'].nunique()

# Plot histogram
plt.figure(figsize=(10,6))
pep_counts.plot(kind='bar')
plt.xlabel('Organism')
plt.ylabel('Number of Unique pep_ids')
plt.title('Unique pep_ids per Organism')
plt.tight_layout()
plt.savefig('/ix/djishnu/Priyamvada/virauto/results/virscan/misc/candidate_pep_id_histogram.png')
plt.close()

# View first few rows
print(post_comp_df.head())

