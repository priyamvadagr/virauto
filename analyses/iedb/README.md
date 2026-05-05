
## Overview

This pipeline identifies microbial epitopes that mimic human self-peptides and could drive autoimmune T cell cross-reactivity. Starting from experimentally validated T cell epitopes in IEDB, it finds sequence-similar human peptides, predicts MHC binding for both, and categorizes mimicry candidates by binding equivalence. The candidates are then linked to TCR sequence data for downstream analysis with DecoderTCR.

### Pipeline Architecture

```
IEDB T cell epitopes
        │
        ├── Query API (MHC Class I + II separately)
        │
        ├── Link to TCR sequences (bulk export)
        │
        ├── BLAST against human proteome
        │
        ├── Filter hits (coverage, mismatches, immune proteins)
        │
        ├── NetMHCpan binding prediction (per-allele, paired viral/human)
        │
        ├── Parse results, compute ΔBA, categorize
        │
        ├── Summary plots and enrichment analysis
        │
        └── [Next] DecoderTCR embedding analysis
```

---

## Step 1: Retrieve Microbial T Cell Epitopes from IEDB

**Script:** `analyses/iedb/scripts/python/query_iedb_epitopes.py`

**What it does:**  
Queries the IEDB IQ-API (`https://query-api.iedb.org/tcell_search`) to retrieve experimentally validated T cell epitopes. Produces separate MHC Class I and Class II datasets.

**API filters (server-side):**
- Host organism: *Homo sapiens* (`NCBITaxon:9606`)
- Assay result: Positive (includes Positive-Low, Positive-High)
- Structure type: Linear peptide
- MHC class: `I` or `II`
- No disease context filter applied — captures infectious, autoimmune, transplant, allergy contexts

**Post-query filters (client-side):**
- Peptide length: 8–14 aa (Class I) or 12–25 aa (Class II)
- Source organism: non-human (microbial)
- MHC allele: must be annotated

**Output:**
```
data/epitopes/iedb/
├── mhc_i/
│   ├── iedb_tcell_mhc_i_epitopes.csv          # ~29,700 records, 7,744 unique epitopes
│   └── iedb_unique_epitopes_mhc_i.csv
└── mhc_ii/
    ├── iedb_tcell_mhc_ii_epitopes.csv          # ~44,000 records, 24,780 unique epitopes
    └── iedb_unique_epitopes_mhc_ii.csv
```

---

## Step 2: Link Epitopes to TCR Sequences

**Script:** `analyses/iedb/scripts/python/merge_iedb_epitopes_with_tcr.py`

**What it does:**  
Merges epitope data with TCR sequences from the IEDB bulk export (`tcr_full_v3.csv`). The bulk export has a two-row header (category + field name) that is parsed into unique column names like `Receptor - IEDB Receptor ID`, `Chain 1 - CDR3 Curated`, etc.

**TCR data source:**  
Downloaded from IEDB Database Export page → CSV Metric Exports → `tcr_full_v3.zip`

**Merge strategy:**
1. Extract `structure_id` from epitope IRI in TCR export (e.g., `https://www.iedb.org/epitope/69921` → `69921`)
2. Join to epitope data on `structure_id`
3. Secondary join on `receptor_id` via `receptor_ids` list in epitope data

**CDR sequence handling:**  
Consolidates curated and calculated versions — prefers curated, falls back to calculated. Output has single columns: `alpha_cdr3`, `beta_cdr3`, `alpha_cdr1`, etc.

**DecoderTCR readiness flags:**
- `decoderTCR_paired_cdr3`: both α and β CDR3 + peptide + MHC allele
- `decoderTCR_full_ready`: above + V/J genes for both chains (Stitchr-ready)

**Create fasta for epitopes to use for BLAST:**
- All epitopes irrespective of TCR availability are written into the fasta file to blast against the human proteome
- FASTA header format: `>structure_id|sequence|organism|mhc_allele`


**Output:**
```
data/epitopes/iedb/
├── mhc_i/
│   ├── iedb_epitopes_with_tcr_mhc_i.csv.gz
│   └── iedb_epitopes_no_tcr_mhc_i.csv.gz
├── mhc_ii/
│   ├── iedb_epitopes_with_tcr_mhc_ii.csv.gz
│   └── iedb_epitopes_no_tcr_mhc_ii.csv.gz
├── all_receptor_data/
│   ├── tcr_full_v3.csv
│   └── tcr_full_v3_parsed.csv.gz
└── fasta/
    ├── iedb_epitopes_mhc_i.fasta 
    ├── iedb_epitopes_mhc_ii.fasta (used iedb_epitopes_mhc_ii_microbial.fasta as there were a lot of other non-human epitopes  ~ 14%)
    └── iedb_epitopes_all.fasta
```

---

## Step 3: BLAST Epitopes Against Human Proteome

**Script:** `analyses/alignment/scripts/slurm/blast_iedb_mhci.sh` (SLURM job)

**What it does:**  
BLASTs MHC Class I epitope sequences against the full UniProt human proteome using settings optimized for short peptides.

**BLAST parameters:**
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `-word_size` | 2 | Maximum sensitivity for short queries |
| `-evalue` | 1000 | Short queries produce high e-values; filter downstream |
| `-matrix` | BLOSUM62 | Standard for peptide comparisons |
| `-seg` | no | Prevents masking short peptide regions |
| `-comp_based_stats` | 0 | Unreliable for very short queries |
| `-ungapped` | — | Gaps change MHC binding register; not biologically meaningful |
| `-max_target_seqs` | 500 | Comprehensive coverage; filter downstream |

**Input:**
- Query: `data/epitopes/iedb/fasta/iedb_epitopes_mhc_i.fasta`
- Database: `data/refs/blastdb/uniprot_human_all/uniprot_human_all_db`

**Output:** `data/epitopes/iedb/blast/iedb_mhci_vs_human_proteome.tsv`

---

## Step 4: Filter BLAST Hits

**Script:** `analyses/iedb/scripts/python/filter_iedb_blast_hits.py`

**What it does:**  
Progressively filters BLAST results to identify biologically plausible mimicry candidates.

**Filters applied (in order):**

| Filter |                        Criterion        | Rationale                                |
|--------|-----------------------------------------|------------------------------------------|
| Coverage | alignment_length / query_length ≥ 80% | Alignment must cover most of the epitope |
| Mismatches | 1 ≤ mismatches ≤ qlen // 3 | Excludes identical (not mimicry) and too distant; ~33% max mismatches scaled by length |
| Deprecated UniProt | Remove IDs from `proteins_to_remove_from_UniProtKB.txt` | Data quality |
| Immune proteins | Exclude HLA, immunoglobulin, TCR, MHC proteins | Spurious matches from peptide-binding domains |

**HLA resolution classification:**
- 4-digit (e.g., `HLA-A*02:01`): NetMHCpan-ready
- Low-resolution (e.g., `HLA-A2`, `human`): excluded from current analysis 

**Output:**
```
data/epitopes/iedb/blast/
├── iedb_mhci_filtered_blast_hits.csv.gz        # All filtered hits
├── iedb_mhci_human_mimic_seqs.fasta            # Human subsequences
├── iedb_mhci_pairs_for_netmhcpan.csv.gz        # All pairs with MHC info
├── iedb_mhci_pairs_4digit_hla.csv.gz           # 33,250 pairs with 4-digit HLA
└── iedb_mhci_pairs_high_confidence.csv.gz
```

---

## Step 5: Prepare NetMHCpan Input

**Script:** `analyses/netmhcpan/iedb/scripts/python/prepare_netmhcpan_iedb.py`

**What it does:**  
Assigns unique 4-character alphanumeric pair IDs, groups pairs by HLA allele, writes per-allele FASTA files, and generates SLURM submission scripts.

**Pair ID system:**
- Each unique viral-human pair (structure_id + hu_prot_id + both sequences) → 4-char code (e.g., `A01B`)
- FASTA headers: `V_A01B` (viral, 6 chars) and `H_A01B` (human, 6 chars)
- Mapping table saved for downstream reconstruction

**Sequence filters:**
- Both sequences ≥ 8 aa
- Standard amino acids only (ACDEFGHIKLMNPQRSTVWY)
- Viral and human can be different lengths (MHC accommodates different lengths at same anchor positions)

**SLURM batching:**
- One array task per HLA allele
- Maximum 100 tasks per batch (cluster limit)
- Master script chains batches with `--dependency=afterany`
- Skip-if-done logic for re-runs

**Output:**
```
data/epitopes/iedb/netmhcpan/
├── fasta_by_allele/
│   ├── HLA_A_01_01.fasta
│   ├── HLA_A_02_01.fasta
│   └── ...
├── allele_manifest.tsv
└── pair_id_map.csv.gz

analyses/netmhcpan/iedb/scripts/slurm/
├── submit_netmhcpan_iedb_batch1.sh
├── submit_netmhcpan_iedb_batch2.sh (if >100 alleles)
└── submit_netmhcpan_iedb_all.sh    (master script)
```

---

## Step 6: Run NetMHCpan

**Script:** `analyses/netmhcpan/iedb/scripts/slurm/submit_netmhcpan_iedb*.sh`

**What it does:**  
Runs NetMHCpan 4.2 on each allele's peptide FASTA. Each array task reads the manifest to get its allele and FASTA file.

**NetMHCpan parameters:**
- `-a`: HLA allele from manifest (converted to `HLA-A02:01` format)
- `-f`: Per-allele FASTA file
- `-l 8,9,10,11,12,13,14`: Score all Class I peptide lengths (prevents default 9-mer windowing)
- `-BA`: Include binding affinity prediction
- `-xls`: Tabular output

**Output format (tab-separated):**
```
Pos  Peptide  ID  core  icore  Score  Rank  BA_score  BA_Rank  Ave  NB
```

Where `Score`/`Rank` = eluted ligand prediction, `BA_score`/`BA_Rank` = binding affinity prediction, `NB` = strong binder flag (EL-based).

**Output:** `results/netmhcpan/iedb/mhc_i/*.xls` (one per allele)

---

## Step 7: Parse Results and Categorize

**Script:** `analyses/netmhcpan/iedb/scripts/python/parse_netmhcpan_iedb_results.py`

**What it does:**  
Parses NetMHCpan output, matches viral-human pairs, computes binding affinity differences, and categorizes mimicry candidates.

**Processing steps:**
1. Parse all `.xls` files; extract allele name from first line
2. Identify VIRAL (`V_`) vs HUMAN (`H_`) from ID column; extract pair_id
3. For each peptide ID, keep only best-binding prediction (lowest BA_Rank) across all length windows
4. Match viral-human pairs on pair_id + allele_file
5. Compute ΔBA_score = BA_score(viral) − BA_score(human)
6. Categorize pairs
7. Merge metadata via pair_id_map

**Categorization criteria:**

| Category | Criteria | Biological meaning |
|----------|----------|--------------------|
| **Non-binder** | Viral BA_Rank > 2 | Viral epitope not presented |
| **Viral-dominant** | Viral BA_Rank ≤ 2 and ΔBA_score ≥ +0.5 | Viral binds stronger — possible immune evasion |
| **Human-dominant** | Viral BA_Rank ≤ 2 and ΔBA_score ≤ −0.5 | Human binds stronger — mimicry/tolerance risk |
| **Equivalent-binding** | Viral BA_Rank ≤ 2 and \|ΔBA_score\| < 0.5 | Similar binding — strongest mimicry candidates |

**Strong mimicry candidates** = Equivalent-binding AND human BA_Rank ≤ 2 (both peptides are binders).

**Output:**
```
results/netmhcpan/iedb/mhc_i/parsed/
├── iedb_mhci_netmhcpan_results.csv.gz         # All scored pairs
├── iedb_mhci_binders.csv.gz                    # Viral binds (BA_Rank ≤ 2)
├── iedb_mhci_mimicry_candidates.csv.gz         # Both bind
├── iedb_mhci_strong_mimicry.csv.gz             # Equivalent + both bind
└── figures/                                     # Summary plots
```

---

## Step 8: Summary Plots

**Script:** `analyses/netmhcpan/iedb/scripts/python/plot_mimicry_summary.py`

**Figures generated:**

| Figure | Content |
|--------|---------|
| 01 | HLA allele distribution (top 20 bar chart + locus pie chart) |
| 02 | Source organism enrichment (total pairs + unique epitopes per organism) |
| 03 | Binding comparison (viral vs human BA_score scatter, ΔBA distribution, BA_Rank scatter) |
| 04 | Peptide properties (length distribution, mismatch counts, percent identity) |
| 05 | HLA × organism heatmap (top alleles vs top organisms) |
| 06 | Human protein enrichment (top 20 human proteins by unique viral epitopes) |
| 07 | Fold enrichment over input (allele + organism frequency in mimicry vs input, controls for representation bias) |

---

## Key Data Files

| File | Description | Records |
|------|-------------|---------|
| `iedb_tcell_mhc_i_epitopes.csv` | All Class I T cell assay records | ~29,700 |
| `iedb_epitopes_with_tcr_mhc_i.csv.gz` | Epitopes linked to TCR sequences | varies |
| `pair_id_map.csv.gz` | 4-char pair ID → structure_id, sequences, protein IDs | ~33,000 |
| `iedb_mhci_pairs_4digit_hla.csv.gz` | Viral-human pairs with 4-digit HLA | 33,250 |
| `iedb_mhci_strong_mimicry.csv.gz` | Strong mimicry candidates | varies |
| `tcr_full_v3_parsed.csv.gz` | Parsed IEDB TCR export | ~190,000 |

---

## Column Reference: Strong Mimicry Output

| Column | Source | Description |
|--------|--------|-------------|
| `pair_id` | prepare script | Unique 4-char alphanumeric pair identifier |
| `allele_file` | NetMHCpan | HLA allele (safe filename format) |
| `viral_peptide` | NetMHCpan | Viral peptide sequence (best-binding window) |
| `viral_rank_BA` | NetMHCpan | Viral binding affinity percentile rank |
| `viral_score_BA` | NetMHCpan | Viral binding affinity score (higher = stronger) |
| `human_peptide` | NetMHCpan | Human mimic sequence (best-binding window) |
| `human_rank_BA` | NetMHCpan | Human binding affinity percentile rank |
| `human_score_BA` | NetMHCpan | Human binding affinity score |
| `delta_BA_score` | parser | BA_score(viral) − BA_score(human) |
| `category` | parser | Viral-dominant / Human-dominant / Equivalent-binding / Non-binder |
| `human_is_binder` | parser | Human BA_Rank ≤ 2 |
| `structure_id` | pair_id_map | IEDB epitope identifier |
| `hu_prot_id` | BLAST | UniProt accession of human mimic protein |
| `hu_prot_name` | BLAST | UniProt entry name |
| `source_organism` | IEDB | Microbial source of viral epitope |
| `mhc_allele` | IEDB | HLA restriction (4-digit) |
| `n_mismatches` | BLAST | Amino acid differences between viral and human |
| `pident` | BLAST | Percent identity of alignment |

---

## Next Steps

1. **Cross-reference mimicry candidates with TCR data** — identify which strong mimicry pairs have associated paired αβ TCR sequences (from IEDB merge or VDJdb)
2. **DecoderTCR embedding analysis** — for pairs with TCR data, extract Stage 2 embeddings conditioned on autoimmune-associated TCRs and compare viral vs human epitope positions in embedding space
3. **Positive controls** — validate embedding approach using experimentally confirmed cross-reactive pairs (MAGE-A3/Titin, EBV-BALF5/MBP, Klebsiella/PPI)
4. **Human protein expression** — check whether human mimic proteins are expressed in autoimmune-relevant tissues (Human Protein Atlas, GTEx)
5. **MHC Class II analysis** — extend pipeline to class II epitopes using NetMHCIIpan

---

## Dependencies

- Python 3.10 (conda env: `virauto`)
- pandas, biopython, requests, tqdm, matplotlib, seaborn
- BLAST+ 2.14.1 (module: `blast-plus/2.14.1`)
- NetMHCpan 4.2 (`/ix/djishnu/Priyamvada/auto_immune/NetMHCpan/`)
- DecoderTCR (conda env: `decoder_tcr`, 650M checkpoint)

## External Data Sources

| Source | URL | Data used |
|--------|-----|-----------|
| IEDB | iedb.org | T cell epitopes, TCR sequences |
| IEDB IQ-API | query-api.iedb.org | Programmatic epitope queries |
| UniProt | uniprot.org | Human proteome FASTA |
| IPD-IMGT/HLA | github.com/ANHIG/IMGTHLA | HLA nomenclature (WMDA `rel_dna_ser.txt`) |
| VDJdb | vdjdb.cdr3.net | TCR-epitope specificity data |
| McPAS-TCR | friedmanlab.weizmann.ac.il/McPAS-TCR/ | Pathology-associated TCR sequences |

---

## References

- DecoderTCR: CZ Biohub, Feb 2026 preprint (bioRxiv 10.64898/2026.02.04.703820)
- IEDB: Vita et al., Nucleic Acids Res 2025, 53:D436-D443
- NetMHCpan 4.2: Reynisson et al., Nucleic Acids Res 2020
- WMDA HLA Dictionary: Holdsworth et al., Tissue Antigens 2009, 73:95-170
- VDJdb: Bagaev et al., Nucleic Acids Res 2020, 48:D1057-D1062