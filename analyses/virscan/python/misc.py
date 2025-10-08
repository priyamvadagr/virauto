# Get list of HLAs not processed and chunk list not processed
import os
import shutil
import pandas as pd
import glob
from collections import Counter

type1c_chunks = glob.glob('/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_C_chunks/*.xls')
# get only type1c chunk filenames that start with 'HLA-A'
type1c_chunks_hla_a = [os.path.basename(p) for p in type1c_chunks if os.path.basename(p).startswith('HLA-A')]
type1c_chunks_hla_a_chunks_done = [os.path.basename(x).split('_matched_pairs_chunk_')[1] for x in type1c_chunks_hla_a]
type1c_chunks_hla_a_processed = [x.split('.xls')[0] for x in type1c_chunks_hla_a_chunks_done]
# get counts
counts_chunks_1c = Counter(type1c_chunks_hla_a_processed)
len(type1c_chunks_hla_a)
# print Counter and a sorted DataFrame view
print(counts_chunks_1c)
counts_df = pd.DataFrame.from_dict(counts_chunks_1c, orient='index', columns=['count']).reset_index().rename(columns={'index': 'element'}).sort_values('count', ascending=False)
print(counts_df)
# Filter chunks with count = 200 
hla_a_chunks_processed = counts_df[counts_df['count'] == 200]['element'].tolist()

# Move files to type_1a directory 
src_dir = '/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_C_chunks'
dst_dir = '/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1/type1_A_chunks'
hla_a_chunks_processed = [f'_matched_pairs_chunk_{x}' for x in hla_a_chunks_processed]
moved = 0
os.makedirs(dst_dir, exist_ok=True)
for chunk in hla_a_chunks_processed:
    pattern = os.path.join(src_dir, f"HLA-A_*_*{chunk}")
    matches = glob.glob(pattern)
    print(len(matches))
    for src_path in matches:
        dst_path = os.path.join(dst_dir, os.path.basename(src_path))
        try:
            shutil.move(src_path, dst_path)
            moved += 1
        except Exception as e:
            print(f"Failed to move {src_path}: {e}")

print(f"Moved {moved} files from {src_dir} to {dst_dir}")