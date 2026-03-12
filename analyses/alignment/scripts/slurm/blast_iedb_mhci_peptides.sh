#!/bin/bash
#SBATCH -J blast_iedb_mhci
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 1-00:00
#SBATCH --mem=16G
#SBATCH --output=/ix/djishnu/Priyamvada/virauto/analyses/alignment/logs/blast_iedb_mhci_epitopes.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

# ============================================================
# BLAST IEDB MHC Class I epitopes against the human proteome
#
# These are short peptides (8-14 aa), so we use settings 
# optimized for short query sequences:
#   -word_size 2        : smallest seed for maximum sensitivity
#   -matrix BLOSUM62    : standard for short peptide comparisons
#   -evalue 1000        : lenient; short queries need high e-value
#                         to capture weak but real similarities
#   -seg no             : disable low-complexity masking (would 
#                         mask parts of short peptides)
#   -comp_based_stats 0 : disable composition-based statistics 
#                         which are unreliable for very short queries
#   -ungapped           : for 8-11mer epitopes, gapped alignments 
#                         are rarely meaningful — we want direct
#                         substitution matches, not indel alignments
#   -max_target_seqs 500: keep more hits since we filter downstream
#
# Output format 6 columns:
#   qseqid   : query epitope ID (structure_id|seq|org|mhc)
#   sseqid   : subject (human protein) ID
#   pident   : percent identity of alignment
#   length   : alignment length
#   mismatch : number of mismatches
#   gapopen  : number of gap openings
#   qlen     : query length (epitope length)
#   slen     : subject length (human protein length)
#   qstart   : alignment start position in query
#   qend     : alignment end position in query
#   sstart   : alignment start position in subject
#   send     : alignment end position in subject
#   evalue   : expect value
#   bitscore : bit score
# ============================================================

module load blast-plus/2.14.1

QUERY="/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/fasta/iedb_epitopes_mhc_i.fasta"
DB="/ix/djishnu/Priyamvada/virauto/data/refs/blastdb/uniprot_human_all/uniprot_human_all_db"
OUTDIR="/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/blast"
mkdir -p ${OUTDIR}

OUT="${OUTDIR}/iedb_mhci_vs_human_proteome.tsv"

echo "============================================================"
echo "BLAST: IEDB MHC-I epitopes vs human proteome"
echo "Query: ${QUERY}"
echo "Database: ${DB}"
echo "Output: ${OUT}"
echo "============================================================"

# Count input sequences
N_SEQS=$(grep -c "^>" ${QUERY})
echo "Input sequences: ${N_SEQS}"
echo "Starting BLAST..."

blastp \
    -query ${QUERY} \
    -db ${DB} \
    -out ${OUT} \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen slen qstart qend sstart send evalue bitscore" \
    -evalue 1000 \
    -word_size 2 \
    -matrix BLOSUM62 \
    -seg no \
    -comp_based_stats 0 \
    -ungapped \
    -max_target_seqs 500 \
    -num_threads 8

echo ""
echo "BLAST complete."
echo "Output lines: $(wc -l < ${OUT})"
echo "Output file: ${OUT}"
echo "============================================================"
echo ""
echo "Next: Filter hits with filter_blast_hits.py"
echo "  - Remove self-hits (same sequence)"  
echo "  - Filter by alignment coverage (alignment_length / qlen)"
echo "  - Run NetMHCpan on human hits using epitope-specific HLAs"
echo "============================================================"