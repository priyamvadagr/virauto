#!/usr/bin/env python3
"""
Map UniProt accessions via API, extract peptide sequences from full proteins,
and write to FASTA for BLAST or alignment.
"""

import pandas as pd
import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from time import sleep

# ========================
# INPUTS
# ========================
virscan_csv = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/peptides_metadata.csv"  # Your CSV with columns: pep_id, UniProt_acc, start, end
out_fasta = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/VirScan_peptides_api.fasta"

# ========================
# UniProt API Functions
# ========================

def fetch_sequence(accession):
    """Fetch full sequence for a UniProt accession via REST API"""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"⚠️ Could not fetch sequence for {accession}")
        return None

    lines = response.text.splitlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    return sequence

def resolve_accessions(accessions):
    """
    Use UniProt ID mapping API to resolve old/truncated/secondary accessions.
    Returns a dictionary: input_id → current_primary_id
    """
    ids = list(set(accessions))
    print(f"[info] Submitting ID mapping job for {len(ids)} accessions")

    url = "https://rest.uniprot.org/idmapping/run"
    params = {
        "from": "UniProtKB_AC-ID",
        "to": "UniProtKB",
        "ids": ",".join(ids)
    }
    r = requests.post(url, data=params)
    r.raise_for_status()
    job_id = r.json()["jobId"]

    # Poll job status
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status = requests.get(status_url).json()
        if status.get("jobStatus") == "RUNNING":
            print("[wait] Mapping job still running...")
            sleep(2)
        else:
            break

    # Fetch results
    result_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
    result = requests.get(result_url).json()

    mapping = {}
    mapping_rows = []

    for entry in result["results"]:
        old_id = entry["from"]
        new_id = entry["to"]["primaryAccession"]
        mapping[old_id] = new_id
        mapping_rows.append((old_id, new_id))

    # Print and save mapping
    print("\n[✅ UniProt ID Mapping Results]")
    for old_id, new_id in sorted(mapping.items()):
        print(f"{old_id} → {new_id}")

    # Save to CSV
    pd.DataFrame(mapping_rows, columns=["Original_ID", "Mapped_ID"]) \
        .to_csv("uniprot_id_mapping.csv", index=False)
    print("📄 Saved mapping table to uniprot_id_mapping.csv")

    return mapping
# ========================
# Main peptide extraction
# ========================

# Load metadata
df = pd.read_csv(virscan_csv)

# Resolve accessions via API
resolved_map = resolve_accessions(df["UniProt_acc"].dropna().astype(str).tolist())

peptide_records = []
failures = []

for _, row in df.iterrows():
    orig_acc = str(row["UniProt_acc"])
    acc = resolved_map.get(orig_acc, orig_acc)
    start = int(row["start"])
    end = int(row["end"])
    pep_id = str(row["pep_id"])

    full_seq = fetch_sequence(acc)
    if full_seq is None:
        print(f"⚠️ {pep_id}: Sequence not found for UniProt ID {acc}")
        continue

    if end > len(full_seq):
        print(f"⚠️ {pep_id}: Range {start}-{end} exceeds sequence length {len(full_seq)} for UniProt {acc}")
        continue

    if full_seq is None or len(full_seq) == 0:
        failures.append((pep_id, acc, start, end, 0, "Sequence not found"))
        continue
    if end > len(full_seq):
        failures.append((pep_id, acc, start, end, len(full_seq), "Out of bounds"))
        continue

    subseq = full_seq[start-1:end]
    header = f"{pep_id}|{acc}|{start}-{end}"
    record = SeqRecord(Seq(subseq), id=header, description="")
    peptide_records.append(record)

# Write to FASTA
SeqIO.write(peptide_records, out_fasta, "fasta")
print(f"✅ Wrote {len(peptide_records)} peptide sequences to {out_fasta}")
# Log failures
pd.DataFrame(failures, columns=["pep_id", "UniProt_acc", "start", "end", "seq_len", "reason"]).to_csv("/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/peptide_sequence_failures.csv", index=False)
