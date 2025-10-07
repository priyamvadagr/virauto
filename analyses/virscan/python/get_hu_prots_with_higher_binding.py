#!/usr/bin/env python3
"""
===============================================================================
Get human peptides with higher binding affinity in candidate allele set 
===============================================================================
"""
#========================================================================
# Import packages and input/output paths
#========================================================================
import os
import pandas as pd
import numpy as np
concat_file_dir = '/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_concat'
out_dir = '/ix/djishnu/Priyamvada/virauto/results/virscan/mimicry_analysis/9_mers/'
candidate_alleles_file = '/ix/djishnu/Priyamvada/virauto/results/virscan/mimicry_analysis/9_mers/data/human_dominant_candidate_alleles.csv'

#========================================================================
# Read input files
#========================================================================
concat_dfs = []
for f in os.listdir(concat_file_dir):
    if f.endswith('_all_predictions_processed.txt'):
        df = pd.read_csv(os.path.join(concat_file_dir, f), sep='\t')
        concat_dfs.append(df)

concat_df = pd.concat(concat_dfs, ignore_index=True)
candidate_alleles = pd.read_csv(candidate_alleles_file)

#========================================================================
# Filter to candidate alleles and human-dominant mimicry class
#========================================================================
candidate_allele_list = candidate_alleles['Allele'].tolist()
filtered_df = concat_df[concat_df['Allele'].isin(candidate_allele_list)]

# --- Define thresholds ---
BA_rank_thresh = 2.0
BA_score_delta = 0.5

# --- Compute ΔBA_score and binding classification ---
filtered_df["BA_score_diff"] = filtered_df["BA_score_viral"] - filtered_df["BA_score_human"]

both_bind = (filtered_df["BA_Rank_viral"] <= BA_rank_thresh) & (
    filtered_df["BA_Rank_human"] <= BA_rank_thresh
)

filtered_df["Mimicry_class"] = "Non_binder"
filtered_df.loc[both_bind & (np.abs(filtered_df["BA_score_diff"]) <= BA_score_delta), "Mimicry_class"] = "Equivalent_binding"
filtered_df.loc[both_bind & (filtered_df["BA_score_diff"] < -BA_score_delta), "Mimicry_class"] = "Viral_dominant"
filtered_df.loc[both_bind & (filtered_df["BA_score_diff"] > BA_score_delta), "Mimicry_class"] = "Human_dominant"

filtered_df['Human_protein'] = filtered_df['Peptide_ID_full_human'].str.split('_').str[3]
filtered_df['Human_protein'].nunique()

#----Filter to human-dominant mimicry class ----
human_dominant_df = filtered_df[filtered_df['Mimicry_class'] == 'Human_dominant']

#========================================================================
# Get uniprot ids of human proteins and plot the number of alleles each protein is enriched
#========================================================================
human_dominant_df['Human_protein'] = human_dominant_df['Peptide_ID_full_human'].str.split('_').str[3]
human_dominant_df_prot_alleles = human_dominant_df[['Human_protein', 'Allele']].drop_duplicates()
# Count the number of alleles per protein
protein_allele_counts = human_dominant_df.groupby('Human_protein').size()
human_dominant_df['Human_protein'].nunique()
import matplotlib.pyplot as plt

# Count alleles per protein and HLA type
protein_hla_counts = human_dominant_df.groupby(['Human_protein', 'Allele']).size().reset_index(name='count')

# Extract HLA type (A, B, or C) from allele name
protein_hla_counts['HLA_type'] = protein_hla_counts['Allele'].str.extract(r'HLA-([ABC])')[0]

# Pivot to get counts by protein and HLA type
hla_breakdown = protein_hla_counts.groupby(['Human_protein', 'HLA_type']).size().unstack(fill_value=0)

# Sort by total count (descending)
hla_breakdown['Total'] = hla_breakdown.sum(axis=1)
hla_breakdown = hla_breakdown.sort_values('Total', ascending=False).drop('Total', axis=1)

# Create stacked bar plot
fig, ax = plt.subplots(figsize=(14, 6))

colors = {'A': '#FF6B6B', 'B': '#4ECDC4', 'C': '#45B7D1'}
hla_breakdown.plot(kind='bar', stacked=True, ax=ax, 
                   color=[colors.get(col, 'gray') for col in hla_breakdown.columns],
                   edgecolor='black', linewidth=0.5)

ax.set_xlabel('Human Protein', fontsize=12)
ax.set_ylabel('Number of Alleles', fontsize=12)
ax.set_title('HLA Allele Distribution by Human Protein', fontsize=14, fontweight='bold')
ax.legend(title='HLA Type', labels=[f'HLA-{col}' for col in hla_breakdown.columns])
plt.xticks(rotation=45, ha='right')
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'plots/protein_allele_distribution_by_HLA.png'), dpi=300, bbox_inches='tight')
plt.show()

# Print summary
print(f"Total proteins: {len(hla_breakdown)}")
print(f"\nTop 10 proteins by allele count:")
print(hla_breakdown.sum(axis=1).head(10))

#========================================================================
# Plot number of unque peptides per protein 
#========================================================================

## Get total unique peptides per protein from filtered_df
total_peptides_per_protein = filtered_df.groupby('Human_protein')['Peptide_human'].nunique()

# Get unique peptides per protein in dominant_df (human_dominant_df)
dominant_peptides_per_protein = human_dominant_df.groupby('Human_protein')['Peptide_human'].nunique()

# Calculate percentage
peptide_percentage = (dominant_peptides_per_protein / total_peptides_per_protein * 100).fillna(0)

# Create summary dataframe
peptide_summary = pd.DataFrame({
    'Total_Peptides': total_peptides_per_protein,
    'Dominant_Peptides': dominant_peptides_per_protein,
    'Percentage': peptide_percentage
}).sort_values('Percentage', ascending=False)

print(peptide_summary)

# Plot the percentages
plt.figure(figsize=(14, 6))
peptide_summary_sorted = peptide_summary.sort_values('Percentage', ascending=False)

plt.bar(range(len(peptide_summary_sorted)), peptide_summary_sorted['Percentage'].values,
        edgecolor='black', alpha=0.7)
plt.xlabel('Protein', fontsize=12)
plt.ylabel('Percentage of Peptides in Dominant Set (%)', fontsize=12)
plt.title('Percentage of Peptides per Protein in Human-Dominant Set', fontsize=14, fontweight='bold')
plt.xticks(range(len(peptide_summary_sorted)), peptide_summary_sorted.index,
           rotation=45, ha='right', fontsize=9)
plt.axhline(y=50, color='r', linestyle='--', alpha=0.5, label='50% threshold')
plt.legend()
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'plots/peptide_percentage_per_protein.png'), dpi=300, bbox_inches='tight')
plt.show()

# Save summary to CSV
peptide_summary.to_csv(os.path.join(out_dir, 'peptide_percentage_summary.csv'))

print(f"\nSummary Statistics:")
print(f"Mean percentage: {peptide_summary['Percentage'].mean():.2f}%")
print(f"Median percentage: {peptide_summary['Percentage'].median():.2f}%")
print(f"Proteins with >50% dominant: {(peptide_summary['Percentage'] > 50).sum()}")


# Analyze HLA allele diversity per peptide
peptide_allele_counts = human_dominant_df.groupby('Peptide_ID_full_human')['Allele'].nunique()

# Plot 3: Distribution of HLA alleles per peptide
plt.figure(figsize=(10, 6))
plt.hist(peptide_allele_counts, bins=range(1, peptide_allele_counts.max() + 2), 
         edgecolor='black', alpha=0.7)
plt.xlabel('Number of HLA Alleles per Peptide', fontsize=12)
plt.ylabel('Number of Peptides', fontsize=12)
plt.title('Distribution of HLA Alleles Associated with Each Peptide', fontsize=14, fontweight='bold')
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'plots/hla_alleles_per_peptide.png'), dpi=300)
plt.show()

print(f"\nPeptide-HLA associations:")
print(f"Total unique peptides: {len(peptide_allele_counts)}")
print(f"Mean alleles per peptide: {peptide_allele_counts.mean():.2f}")
print(f"Max alleles for one peptide: {peptide_allele_counts.max()}")

# Find peptides with most alleles
print(f"\nTop 10 peptides by allele count:")
top_peptides = peptide_allele_counts.sort_values(ascending=False).head(10)
for peptide, count in top_peptides.items():
    print(f"{peptide}: {count} alleles")

import sys

# Open file and redirect stdout
# Open file and write directly
with open(os.path.join(out_dir, 'peptide_analysis_output.txt'), 'w') as f:
    # Detailed analysis: For each protein, show peptide-allele relationship
    f.write("\n" + "="*80 + "\n")
    f.write("Peptide diversity analysis by protein:\n")
    f.write("="*80 + "\n")
    
    for protein in dominant_peptides_per_protein.sort_values(ascending=False).index:
        # Data from human_dominant_df
        protein_data = human_dominant_df[human_dominant_df['Human_protein'] == protein]
        n_peptides_dominant = protein_data['Peptide_human'].nunique()
        n_alleles = protein_data['Allele'].nunique()
        
        # Data from filtered_df (total)
        protein_data_total = filtered_df[filtered_df['Human_protein'] == protein]
        n_peptides_total = protein_data_total['Peptide_human'].nunique()
        
        # Calculate percentage
        percentage = (n_peptides_dominant / n_peptides_total * 100) if n_peptides_total > 0 else 0
        
        f.write(f"\nProtein: {protein}\n")
        f.write(f"  Total unique peptides (filtered_df): {n_peptides_total}\n")
        f.write(f"  Dominant unique peptides (human_dominant_df): {n_peptides_dominant}\n")
        f.write(f"  Percentage dominant: {percentage:.2f}%\n")
        f.write(f"  Unique alleles: {n_alleles}\n")
        
        # Show if different peptides have different alleles
        peptide_allele_breakdown = protein_data.groupby('Peptide_ID_full_human')['Allele'].apply(
            lambda x: f"{x.nunique()} alleles"
        )
        f.write(f"  Peptide breakdown (top 3 dominant peptides):\n")
        for peptide, info in peptide_allele_breakdown.head(3).items():
            f.write(f"    {peptide}: {info}\n")

print("Analysis saved to file!")

# Create a heatmap for a specific protein
protein_of_interest = dominant_peptides_per_protein.idxmax()  # Protein with most peptides

protein_subset = human_dominant_df[human_dominant_df['Human_protein'] == protein_of_interest]

# Create peptide-allele matrix
peptide_allele_matrix = protein_subset.groupby(['Peptide_human', 'Allele']).size().unstack(fill_value=0)
peptide_allele_matrix_v2 = protein_subset.groupby(['Peptide_ID_full_human', 'Allele']).size().unstack(fill_value=0)
peptide_allele_matrix.shape

# Plot heatmap
plt.figure(figsize=(14, 8))
plt.imshow(peptide_allele_matrix.values > 0, cmap='YlOrRd', aspect='auto')
plt.colorbar(label='Association', ticks=[0, 1])
plt.xlabel('HLA Allele', fontsize=12)
plt.ylabel('Peptide', fontsize=12)
plt.title(f'Peptide-Allele Associations for {protein_of_interest}', fontsize=14, fontweight='bold')
plt.xticks(range(len(peptide_allele_matrix.columns)), peptide_allele_matrix.columns, rotation=90, fontsize=8)
plt.yticks(range(len(peptide_allele_matrix.index)), peptide_allele_matrix.index, fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, f'plots/peptide_allele_heatmap_{protein_of_interest}.png'), dpi=300, bbox_inches='tight')
plt.show()