#!/usr/bin/env python3
"""
Build a k-mer index from a protein FASTA.

Usage:
  python build_kmer_index.py --db-fasta uniprot_human_all.fasta --k 8 --out-index human_k8.pkl
"""

import argparse
import gzip
import pickle
from collections import defaultdict
from Bio import SeqIO

AA20 = set("ACDEFGHIKLMNPQRSTVWY")

def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def valid_kmer(kmer):
    return set(kmer).issubset(AA20)

def iter_kmers(seq, k):
    L = len(seq)
    for i in range(0, L - k + 1):
        kmer = seq[i:i+k]
        if valid_kmer(kmer):
            yield i, kmer

def build_index(db_fasta, k, out_index):
    index = defaultdict(list)
    n_prots, n_kmers = 0, 0

    with open_maybe_gzip(db_fasta, "rt") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            pid = rec.id.split("|")[1] if "|" in rec.id else rec.id
            seq = str(rec.seq).upper()
            for pos0, kmer in iter_kmers(seq, k):
                index[kmer].append((pid, pos0))
                n_kmers += 1
            n_prots += 1

    with open(out_index, "wb") as f:
        pickle.dump({"k": k, "index": index}, f, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"[index] Proteins: {n_prots}")
    print(f"[index] Stored k-mers: {n_kmers}")
    print(f"[index] Wrote index: {out_index}")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--db-fasta", required=True)
    p.add_argument("--k", type=int, default=8)
    p.add_argument("--out-index", required=True)
    args = p.parse_args()

    build_index(args.db_fasta, args.k, args.out_index)
