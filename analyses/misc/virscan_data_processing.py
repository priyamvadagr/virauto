import pandas as pd

# Path to your large CSV file
csv_file = '/ix/djishnu/Priyamvada/virauto/data/epitopes/Virscan_dataset.csv'


# Read a small chunk (e.g., first 100 rows)
chunk_size = 10000
df_chunk = pd.read_csv(csv_file, nrows=chunk_size)

print(df_chunk.columns.tolist())
print(df_chunk['serotype'].tolist())
df_chunk.to_csv('/ix/djishnu/Priyamvada/virauto/data/epitopes/Virscan_dataset_chunk.csv', index=False)