#!/usr/bin/env python3
"""
======================================================================
Script: query_iedb_epitopes_and_receptors.py
Description:
    Programmatically query the IEDB IQ-API to retrieve:
      1. Microbial T cell epitopes with positive assay results
         and MHC restriction data, from human hosts.
      2. Associated TCR (receptor) sequences for those epitopes.
    
    Produces separate MHC class I and class II datasets for 
    downstream molecular mimicry analysis.

API Reference:
    Base URL: https://query-api.iedb.org
    Docs: https://query-api.iedb.org (Swagger/OpenAPI)
    Syntax: PostgREST (https://postgrest.org/en/stable/references/api.html)
    
    Key endpoints:
      /epitope_search   - query epitopes
      /tcell_search     - query T cell assays (one row per assay)
      /receptor_search  - query receptor (TCR/BCR) data

Notes:
    - No API key required
    - Default page size limit is 10,000 records; use offset for pagination
    - Results returned as JSON (default) or TSV (via Accept header)

Dependencies:
    requests, pandas, tqdm

Usage:
    python query_iedb_epitopes_receptors.py
======================================================================
"""

import requests
import pandas as pd
import time
import os
from tqdm import tqdm

# ====================================================
# Config
# ====================================================
BASE_URL = "https://query-api.iedb.org"
OUT_DIR = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb"
os.makedirs(OUT_DIR, exist_ok=True)

# Separate output directories for each MHC class
for subdir in ["mhc_i", "mhc_ii"]:
    os.makedirs(os.path.join(OUT_DIR, subdir), exist_ok=True)

PAGE_SIZE = 1  # IEDB max per request
RATE_LIMIT_DELAY = 1  # seconds between requests

# MHC class-specific peptide length ranges
MHC_CLASS_CONFIG = {
    "mhc_i": {
        "label": "MHC Class I",
        "min_length": 8,
        "max_length": 14,   # most are 8-11, allow up to 14
        "mhc_class_filter": "Class I",
    },
    "mhc_ii": {
        "label": "MHC Class II",
        "min_length": 12,
        "max_length": 25,   # typically 13-25
        "mhc_class_filter": "Class II",
    }
}


# ====================================================
# Helper: paginated API query
# ====================================================
def query_iedb(endpoint, params=None, max_records=None):
    """
    Query an IEDB IQ-API endpoint with pagination.
    
    Parameters
    ----------
    endpoint : str
        API endpoint (e.g., 'tcell_search', 'receptor_search')
    params : dict
        PostgREST query parameters
    max_records : int or None
        Stop after this many records (None = get all)
    
    Returns
    -------
    pd.DataFrame
    """
    if params is None:
        params = {}
    
    all_records = []
    offset = 0
    
    while True:
        paginated_params = {**params, "limit": PAGE_SIZE, "offset": offset}
        url = f"{BASE_URL}/{endpoint}"
        
        print(f"    Querying {endpoint} (offset={offset})...")
        
        try:
            resp = requests.get(url, params=paginated_params, timeout=120)
            resp.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"    ⚠️ Request failed: {e}")
            break
        
        data = resp.json()
        
        if not data:
            break
        
        all_records.extend(data)
        print(f"    Retrieved {len(data)} records (total: {len(all_records)})")
        
        if len(data) < PAGE_SIZE:
            break
        
        if max_records and len(all_records) >= max_records:
            all_records = all_records[:max_records]
            break
        
        offset += PAGE_SIZE
        time.sleep(RATE_LIMIT_DELAY)
    
    if not all_records:
        return pd.DataFrame()
    
    return pd.DataFrame(all_records)


# ====================================================
# Step 1: Query T cell assays for microbial epitopes
# ====================================================
def fetch_tcell_epitopes(mhc_class_key):
    """
    Fetch T cell assay data for linear peptide epitopes from 
    microbial organisms, tested in human hosts, with positive results,
    restricted to a specific MHC class.
    
    Parameters
    ----------
    mhc_class_key : str
        Either 'mhc_i' or 'mhc_ii'
    
    Returns
    -------
    pd.DataFrame
    """
    config = MHC_CLASS_CONFIG[mhc_class_key]
    
    print(f"\n{'='*60}")
    print(f"Step 1: Fetching T cell assays — {config['label']}")
    print(f"{'='*60}")
    
    select_fields = ",".join([
        "structure_id",
        "linear_sequence",
        "structure_type",
        "parent_source_antigen_names",
        "parent_source_antigen_source_org_names",
        "parent_source_antigen_source_org_iris",
        "host_organism_names",
        "host_organism_iris",
        "mhc_allele_names",
        "mhc_class",
        "qualitative_measure",
        "assay_type",
        "reference_iris",
        "receptor_ids"
    ])
    
    params = {
        "select": select_fields,
        "host_organism_iris": "cs.{http://purl.obolibrary.org/obo/NCBITaxon_9606}",
        "qualitative_measure": "eq.Positive",
        "structure_type": "eq.Linear peptide",
        "mhc_class": f"eq.{config['mhc_class_filter']}",
    }
    
    df = query_iedb("tcell_search", params)
    
    if df.empty:
        print(f"  ⚠️ No records returned for {config['label']}.")
        return df
    
    print(f"\n  Raw records: {len(df):,}")
    
    # --------------------------------------------------
    # Post-query filters
    # --------------------------------------------------
    
    # Filter by peptide length appropriate for MHC class
    df["peptide_length"] = df["linear_sequence"].str.len()
    df = df[
        (df["peptide_length"] >= config["min_length"]) & 
        (df["peptide_length"] <= config["max_length"])
    ]
    print(f"  After length filter ({config['min_length']}-{config['max_length']} aa): {len(df):,}")
    
    # Exclude epitopes from Homo sapiens (we want microbial sources)
    def is_non_human_source(org_iris):
        if not isinstance(org_iris, list):
            return False
        human_iri = "http://purl.obolibrary.org/obo/NCBITaxon_9606"
        return any(iri != human_iri for iri in org_iris)
    
    df = df[df["parent_source_antigen_source_org_iris"].apply(is_non_human_source)]
    print(f"  After excluding human-source epitopes: {len(df):,}")
    
    # Require MHC restriction data
    df = df[df["mhc_allele_names"].apply(lambda x: isinstance(x, list) and len(x) > 0)]
    print(f"  After requiring MHC restriction: {len(df):,}")
    
    # Create string versions for dedup and downstream use
    df["mhc_allele_str"] = df["mhc_allele_names"].apply(
        lambda x: "|".join(sorted(x)) if isinstance(x, list) else ""
    )
    df["source_org_str"] = df["parent_source_antigen_source_org_names"].apply(
        lambda x: "|".join(sorted(set(x))) if isinstance(x, list) else ""
    )
    df["antigen_str"] = df["parent_source_antigen_names"].apply(
        lambda x: "|".join(sorted(set(x))) if isinstance(x, list) else ""
    )
    
    # Flag records with associated TCR data
    df["has_receptor"] = df["receptor_ids"].apply(
        lambda x: isinstance(x, list) and len(x) > 0
    )
    
    # Add MHC class label
    df["mhc_class_label"] = config["label"]
    
    print(f"  Records with TCR data: {df['has_receptor'].sum():,}")
    print(f"  Unique epitope sequences: {df['linear_sequence'].nunique():,}")
    
    # Save
    out_path = os.path.join(OUT_DIR, mhc_class_key, f"iedb_tcell_{mhc_class_key}_epitopes.csv")
    df.to_csv(out_path, index=False)
    print(f"\n  ✅ T cell data saved to {out_path}")
    
    return df


# ====================================================
# Step 2: Fetch receptor (TCR) sequences
# ====================================================
def fetch_receptor_data(tcell_df, mhc_class_key):
    """
    For epitopes that have associated receptor IDs, query the
    receptor_search endpoint to get full TCR sequence data.
    """
    config = MHC_CLASS_CONFIG[mhc_class_key]
    
    print(f"\n{'='*60}")
    print(f"Step 2: Fetching TCR sequences — {config['label']}")
    print(f"{'='*60}")
    
    # Collect unique receptor IDs
    all_receptor_ids = set()
    for ids in tcell_df["receptor_ids"].dropna():
        if isinstance(ids, list):
            all_receptor_ids.update(ids)
    
    all_receptor_ids = sorted(all_receptor_ids)
    print(f"  Unique receptor IDs to query: {len(all_receptor_ids):,}")
    
    if not all_receptor_ids:
        print("  ⚠️ No receptor IDs found.")
        return pd.DataFrame()
    
    select_fields = ",".join([
        "receptor_id",
        "receptor_type",
        "receptor_species_names",
        "chain_1_cdr3_sequence",
        "chain_2_cdr3_sequence",
        "chain_1_v_gene",
        "chain_1_j_gene",
        "chain_2_v_gene",
        "chain_2_j_gene",
        "chain_1_full_sequence",
        "chain_2_full_sequence",
        "chain_1_cdr1_sequence",
        "chain_1_cdr2_sequence",
        "chain_2_cdr1_sequence",
        "chain_2_cdr2_sequence",
        "epitope_iris",
        "epitope_descriptions"
    ])
    
    all_receptors = []
    batch_size = 100
    
    for i in tqdm(range(0, len(all_receptor_ids), batch_size),
                  desc=f"  Fetching receptors ({config['label']})"):
        batch_ids = all_receptor_ids[i:i + batch_size]
        id_list = ",".join(str(rid) for rid in batch_ids)
        
        params = {
            "select": select_fields,
            "receptor_id": f"in.({id_list})"
        }
        
        try:
            resp = requests.get(
                f"{BASE_URL}/receptor_search",
                params=params,
                timeout=120
            )
            resp.raise_for_status()
            data = resp.json()
            if data:
                all_receptors.extend(data)
        except requests.exceptions.RequestException as e:
            print(f"    ⚠️ Batch failed: {e}")
        
        time.sleep(RATE_LIMIT_DELAY)
    
    receptor_df = pd.DataFrame(all_receptors)
    
    if receptor_df.empty:
        print("  ⚠️ No receptor data returned.")
        return receptor_df
    
    print(f"\n  Total receptor records: {len(receptor_df):,}")
    print(f"  Unique receptors: {receptor_df['receptor_id'].nunique():,}")
    
    has_chain1 = receptor_df["chain_1_cdr3_sequence"].notna()
    has_chain2 = receptor_df["chain_2_cdr3_sequence"].notna()
    print(f"  With chain 1 (alpha/gamma) CDR3: {has_chain1.sum():,}")
    print(f"  With chain 2 (beta/delta) CDR3: {has_chain2.sum():,}")
    print(f"  With paired chains: {(has_chain1 & has_chain2).sum():,}")
    
    out_path = os.path.join(OUT_DIR, mhc_class_key, f"iedb_receptors_{mhc_class_key}.csv")
    receptor_df.to_csv(out_path, index=False)
    print(f"\n  ✅ Receptor data saved to {out_path}")
    
    return receptor_df


# ====================================================
# Step 3: Merge epitope and receptor data
# ====================================================
def merge_epitope_receptor(tcell_df, receptor_df, mhc_class_key):
    """
    Link T cell epitope records to their TCR sequences.
    """
    config = MHC_CLASS_CONFIG[mhc_class_key]
    
    print(f"\n{'='*60}")
    print(f"Step 3: Merging epitope + receptor — {config['label']}")
    print(f"{'='*60}")
    
    if receptor_df.empty:
        print("  ⚠️ No receptor data to merge.")
        return pd.DataFrame()
    
    tcell_with_receptors = tcell_df[tcell_df["has_receptor"]].copy()
    tcell_exploded = tcell_with_receptors.explode("receptor_ids")
    tcell_exploded = tcell_exploded.rename(columns={"receptor_ids": "receptor_id"})
    
    merged = tcell_exploded.merge(
        receptor_df,
        on="receptor_id",
        how="inner",
        suffixes=("_epitope", "_receptor")
    )
    
    print(f"  Merged records: {len(merged):,}")
    print(f"  Unique epitopes with TCR: {merged['linear_sequence'].nunique():,}")
    print(f"  Unique TCRs: {merged['receptor_id'].nunique():,}")
    
    # Deduplicate
    dedup_cols = ["linear_sequence", "mhc_allele_str", "receptor_id"]
    merged_dedup = merged.drop_duplicates(subset=dedup_cols)
    print(f"  After dedup: {len(merged_dedup):,}")
    
    # Save merged
    out_path = os.path.join(OUT_DIR, mhc_class_key, f"iedb_epitopes_with_tcr_{mhc_class_key}.csv")
    merged_dedup.to_csv(out_path, index=False)
    print(f"\n  ✅ Merged data saved to {out_path}")
    
    return merged_dedup


# ====================================================
# Step 4: Save epitopes without TCR (for VDJdb cross-ref)
# ====================================================
def save_epitopes_without_tcr(tcell_df, mhc_class_key):
    """Save epitopes that lack TCR data in IEDB for VDJdb matching."""
    config = MHC_CLASS_CONFIG[mhc_class_key]
    
    no_receptor = tcell_df[~tcell_df["has_receptor"]].copy()
    no_receptor_dedup = no_receptor.drop_duplicates(
        subset=["linear_sequence", "mhc_allele_str"]
    )
    
    out_path = os.path.join(OUT_DIR, mhc_class_key, f"iedb_epitopes_no_tcr_{mhc_class_key}.csv")
    no_receptor_dedup.to_csv(out_path, index=False)
    
    print(f"  📌 {config['label']} epitopes without TCR (for VDJdb cross-ref): "
          f"{no_receptor_dedup['linear_sequence'].nunique():,} unique sequences")
    print(f"     Saved to {out_path}")
    
    return no_receptor_dedup


# ====================================================
# Summary
# ====================================================
def print_summary(merged_df, no_tcr_df, mhc_class_key):
    """Print summary statistics for a given MHC class dataset."""
    config = MHC_CLASS_CONFIG[mhc_class_key]
    
    print(f"\n{'='*60}")
    print(f"Summary: {config['label']}")
    print(f"{'='*60}")
    
    # Epitopes with TCR
    if not merged_df.empty:
        has_paired = (
            merged_df["chain_1_cdr3_sequence"].notna() &
            merged_df["chain_2_cdr3_sequence"].notna()
        )
        
        print(f"  Epitope-TCR pairs: {len(merged_df):,}")
        print(f"    With paired chains: {has_paired.sum():,}")
        print(f"    Single chain only: {(~has_paired).sum():,}")
        print(f"  Unique epitopes: {merged_df['linear_sequence'].nunique():,}")
        print(f"  Unique TCRs: {merged_df['receptor_id'].nunique():,}")
        
        print(f"\n  Top 10 source organisms:")
        orgs = merged_df["source_org_str"].value_counts().head(10)
        for org, count in orgs.items():
            print(f"    {org}: {count}")
        
        print(f"\n  Top 10 MHC alleles:")
        alleles = merged_df["mhc_allele_str"].value_counts().head(10)
        for allele, count in alleles.items():
            print(f"    {allele}: {count}")
        
        print(f"\n  Peptide length distribution:")
        print(merged_df["peptide_length"].value_counts().sort_index().to_string())
    else:
        print("  No epitope-TCR pairs found.")
    
    # Epitopes without TCR
    if not no_tcr_df.empty:
        print(f"\n  Epitopes without TCR (available for VDJdb matching):")
        print(f"    Unique sequences: {no_tcr_df['linear_sequence'].nunique():,}")


# ====================================================
# Main
# ====================================================
if __name__ == "__main__":
    print("="*60)
    print("IEDB IQ-API Query Pipeline")
    print("Microbial T cell epitopes with TCR data")
    print("Separate MHC Class I and Class II datasets")
    print("="*60)
    
    results = {}
    
    for mhc_class_key in ["mhc_i", "mhc_ii"]:
        config = MHC_CLASS_CONFIG[mhc_class_key]
        print(f"\n\n{'#'*60}")
        print(f"# Processing: {config['label']}")
        print(f"{'#'*60}")
        
        # Step 1: Fetch T cell epitope data
        tcell_df = fetch_tcell_epitopes(mhc_class_key)
        
        if tcell_df.empty:
            print(f"\n  ❌ No data for {config['label']}. Skipping.")
            results[mhc_class_key] = {
                "tcell": pd.DataFrame(),
                "receptor": pd.DataFrame(),
                "merged": pd.DataFrame(),
                "no_tcr": pd.DataFrame()
            }
            continue
        
        # Step 2: Fetch receptor sequences
        receptor_df = fetch_receptor_data(tcell_df, mhc_class_key)
        
        # Step 3: Merge
        merged_df = merge_epitope_receptor(tcell_df, receptor_df, mhc_class_key)
        
        # Step 4: Save epitopes without TCR
        no_tcr_df = save_epitopes_without_tcr(tcell_df, mhc_class_key)
        
        results[mhc_class_key] = {
            "tcell": tcell_df,
            "receptor": receptor_df,
            "merged": merged_df,
            "no_tcr": no_tcr_df
        }
    
    # --------------------------------------------------
    # Print summaries for both classes
    # --------------------------------------------------
    print(f"\n\n{'#'*60}")
    print(f"# FINAL SUMMARIES")
    print(f"{'#'*60}")
    
    for mhc_class_key in ["mhc_i", "mhc_ii"]:
        r = results[mhc_class_key]
        print_summary(r["merged"], r["no_tcr"], mhc_class_key)
    
    # --------------------------------------------------
    # Combined stats
    # --------------------------------------------------
    print(f"\n{'='*60}")
    print("Combined Statistics")
    print(f"{'='*60}")
    
    total_epitopes_with_tcr = sum(
        r["merged"]["linear_sequence"].nunique()
        for r in results.values() if not r["merged"].empty
    )
    total_epitopes_no_tcr = sum(
        r["no_tcr"]["linear_sequence"].nunique()
        for r in results.values() if not r["no_tcr"].empty
    )
    total_tcrs = sum(
        r["merged"]["receptor_id"].nunique()
        for r in results.values() if not r["merged"].empty
    )
    
    print(f"  Total unique epitopes with TCR data: {total_epitopes_with_tcr:,}")
    print(f"  Total unique epitopes without TCR (for VDJdb): {total_epitopes_no_tcr:,}")
    print(f"  Total unique TCRs: {total_tcrs:,}")
    
    print(f"\n  Output directory: {OUT_DIR}")
    print(f"    mhc_i/  — MHC Class I epitopes, receptors, merged")
    print(f"    mhc_ii/ — MHC Class II epitopes, receptors, merged")
    
    print(f"\n  Next steps:")
    print(f"    1. BLAST epitope sequences against human proteome")
    print(f"    2. Cross-reference epitopes_no_tcr files with VDJdb")
    print(f"    3. Run NetMHCpan (Class I) / NetMHCIIpan (Class II) on mimicry pairs")
    print(f"    4. Feed into DecoderTCR embedding analysis")
    
    print(f"\n✅ Pipeline complete.")