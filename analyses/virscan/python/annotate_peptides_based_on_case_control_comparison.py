import pandas as pd
from scipy.stats import fisher_exact

# Read annotation and filtered data
annote_df = pd.read_csv('/ix/djishnu/Priyamvada/virauto/data/epitopes/annotatced_peptides_post_pre_combined.csv')
filtered_pep_df = pd.read_csv('/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/similarity_filtered_blast_out.csv')

# Extract numeric pep_id
filtered_pep_df['pep_id'] = pd.to_numeric(filtered_pep_df['qseqid'].str.split('|').str[0], errors='coerce')

# Merge annotations
merged = filtered_pep_df.merge(
    annote_df[['pep_id', 'candidate']], 
    on='pep_id', 
    how='left'
)

# Build contingency table:
#   Rows: In filtered vs not filtered
#   Cols: Candidate vs non-candidate
all_peptides = annote_df.copy()
all_peptides['in_filtered'] = all_peptides['pep_id'].isin(merged['pep_id'])

# Counts
a = ((all_peptides['candidate'] == True) & (all_peptides['in_filtered'])).sum()   # candidate in filtered
b = ((all_peptides['candidate'] == False) & (all_peptides['in_filtered'])).sum()  # non-candidate in filtered
c = ((all_peptides['candidate'] == True) & (~all_peptides['in_filtered'])).sum()  # candidate not filtered
d = ((all_peptides['candidate'] == False) & (~all_peptides['in_filtered'])).sum() # non-candidate not filtered

contingency = [[a, b],
               [c, d]]

# Run Fisher’s exact test
oddsratio, pvalue = fisher_exact(contingency, alternative='greater')  # "greater" = enrichment of candidates

print("Contingency table:")
print(pd.DataFrame(contingency, 
                   index=["Candidate","Non-candidate"], 
                   columns=["In filtered","Not filtered"]))

print(f"\nOdds ratio: {oddsratio:.3f}")
print(f"P-value: {pvalue:.3e}")
