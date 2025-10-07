import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

df = pd.read_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/paired_viral_human_9mers.csv")

# Create list to hold alternating viral-human pairs
all_records = []

for _, row in df.iterrows():
    pep_id = str(row["pep_id"])
    offset = str(row["offset"])
    vir_uniprot_id = str(row["viral_protein"]) 
    hu_uniprot_id = str(row["human_protein"])
    v_seq = str(row["viral_kmer"]).upper()
    h_seq = str(row["human_kmer"]).upper()
    
    # Create viral and human records
    viral_record = SeqRecord(Seq(v_seq), id=f"VIRAL_{pep_id}_{offset}_{vir_uniprot_id}", description="")
    human_record = SeqRecord(Seq(h_seq), id=f"HUMAN_{pep_id}_{offset}_{hu_uniprot_id}", description="")
    
    # Add them consecutively (viral followed by human)
    all_records.append(viral_record)
    all_records.append(human_record)

# Write records in chunks of 100
chunk_size = 100
output_dir = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/chunks/"

for i in range(0, len(all_records), chunk_size):
    chunk = all_records[i:i + chunk_size]
    chunk_num = i // chunk_size + 1
    output_file = f"{output_dir}matched_pairs_chunk_{chunk_num}.fasta"
    SeqIO.write(chunk, output_file, "fasta")
    print(f"Wrote chunk {chunk_num} with {len(chunk)} records to {output_file}")

print(f"\nTotal: {len(all_records)} records ({len(all_records)//2} viral-human pairs) written across {(len(all_records)-1)//chunk_size + 1} chunk files.")