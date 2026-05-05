# %%
import pandas as pd

df = pd.read_csv("/ix/djishnu/Priyamvada/virauto/data/tcr/decoder_tcr/decoder_tcr_input.csv")

print(f"Total rows: {len(df):,}")
print(f"\nNull counts per column:")
print(df.isnull().sum())

# Show the first few rows with any null
null_rows = df[df.isnull().any(axis=1)]
print(f"\nRows with any null: {len(null_rows):,}")
if not null_rows.empty:
    print(null_rows.head())

# Also check for empty strings (which can also cause issues)
for col in df.columns:
    n_empty = (df[col].astype(str).str.strip() == "").sum()
    if n_empty:
        print(f"  Empty strings in {col}: {n_empty:,}")
# %%
