import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

df = pd.read_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/paired_viral_human_9mers.csv")

viral_records = []
human_records = []

for _, row in df.iterrows():
    pep_id = str(row["pep_id"])
    offset = str(row["offset"])
    vir_uniprot_id = str(row["viral_protein"]) 
    hu_uniprot_id = str(row["human_protein"])
    v_seq = str(row["viral_kmer"]).upper()
    h_seq = str(row["human_kmer"]).upper()
    viral_records.append(SeqRecord(Seq(v_seq), id=f"VIRAL_{pep_id}_{offset}_{vir_uniprot_id}", description=""))
    human_records.append(SeqRecord(Seq(h_seq), id=f"HUMAN_{pep_id}_{offset}_{hu_uniprot_id}", description=""))

SeqIO.write(viral_records, "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/viral_9mers.fasta", "fasta")
SeqIO.write(human_records, "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers/human_9mers.fasta", "fasta")
print(f"Wrote {len(viral_records)} viral and {len(human_records)} human 9-mers.")
