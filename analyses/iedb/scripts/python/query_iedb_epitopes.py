#!/usr/bin/env python3
"""
======================================================================
Script: query_iedb_epitopes.py
Description:
    Programmatically query the IEDB IQ-API to retrieve microbial T cell 
    epitopes with positive assay results and MHC restriction data, from 
    human hosts. Produces separate MHC class I and class II datasets.

    TCR linkage is handled separately via the IEDB bulk receptor_full_v3 
    CSV export, which maps every receptor to every epitope it binds.

API Reference:
    Base URL: https://query-api.iedb.org
    Endpoint: /tcell_search
    Syntax: PostgREST

Dependencies:
    requests, pandas

Usage:
    python query_iedb_epitopes.py
======================================================================
"""

import requests
import pandas as pd
import time
import os

# ====================================================
# Config
# ====================================================
BASE_URL = "https://query-api.iedb.org"
OUT_DIR = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb"
os.makedirs(OUT_DIR, exist_ok=True)

for subdir in ["mhc_i", "mhc_ii"]:
    os.makedirs(os.path.join(OUT_DIR, subdir), exist_ok=True)

PAGE_SIZE = 10000
RATE_LIMIT_DELAY = 1

HUMAN_IRI = "NCBITaxon:9606"

MHC_CLASS_CONFIG = {
    "mhc_i": {
        "label": "MHC Class I",
        "min_length": 8,
        "max_length": 14,
        "mhc_class_value": "I",
    },
    "mhc_ii": {
        "label": "MHC Class II",
        "min_length": 12,
        "max_length": 25,
        "mhc_class_value": "II",
    }
}


# ====================================================
# Helper: paginated API query
# ====================================================
def query_iedb(endpoint, params=None, order_by=None):
    """
    Query an IEDB IQ-API endpoint with pagination.
    IEDB requires 'order' parameter whenever 'offset' is used.
    """
    if params is None:
        params = {}
    
    if order_by is None:
        order_defaults = {
            "tcell_search": "tcell_id",
            "epitope_search": "structure_id",
            "reference_search": "reference_id",
        }
        order_by = order_defaults.get(endpoint, "structure_id")
    
    all_records = []
    offset = 0
    
    while True:
        paginated_params = {
            **params,
            "limit": PAGE_SIZE,
            "offset": offset,
            "order": order_by,
        }
        url = f"{BASE_URL}/{endpoint}"
        
        print(f"    Querying {endpoint} (offset={offset})...")
        
        try:
            resp = requests.get(url, params=paginated_params, timeout=300)
            resp.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"    ⚠️ Request failed: {e}")
            if hasattr(e, 'response') and e.response is not None:
                try:
                    print(f"    Response: {e.response.json()}")
                except:
                    print(f"    Response: {e.response.text[:300]}")
            break
        
        data = resp.json()
        
        if not data:
            break
        
        all_records.extend(data)
        print(f"    Retrieved {len(data)} records (total: {len(all_records)})")
        
        if len(data) < PAGE_SIZE:
            break
        
        offset += PAGE_SIZE
        time.sleep(RATE_LIMIT_DELAY)
    
    if not all_records:
        return pd.DataFrame()
    
    return pd.DataFrame(all_records)


# ====================================================
# Fetch T cell assays for microbial epitopes
# ====================================================
def fetch_tcell_epitopes(mhc_class_key):
    """
    Fetch T cell assay data for linear peptide epitopes from 
    microbial organisms, tested in human hosts, with positive results.
    """
    config = MHC_CLASS_CONFIG[mhc_class_key]
    
    print(f"\n{'='*60}")
    print(f"Fetching T cell assays — {config['label']}")
    print(f"{'='*60}")
    
    select_fields = ",".join([
        "structure_id",
        "structure_iri",
        "linear_sequence",
        "linear_sequence_length",
        "structure_type",
        "parent_source_antigen_name",
        "parent_source_antigen_iri",
        "parent_source_antigen_source_org_name",
        "parent_source_antigen_source_org_iri",
        "source_organism_name",
        "source_organism_iri",
        "host_organism_name",
        "host_organism_iri",
        "mhc_allele_name",
        "mhc_restriction",
        "mhc_class",
        "qualitative_measure",
        "assay_names",
        "reference_iri",
        "reference_id",
        "pubmed_id",
        "receptor_ids",
        "tcell_id",
    ])
    
    params = {
        "select": select_fields,
        "host_organism_iri": f"eq.{HUMAN_IRI}",
        "qualitative_measure": "like.Positive*",
        "structure_type": "eq.Linear peptide",
        "mhc_class": f"eq.{config['mhc_class_value']}",
    }
    
    df = query_iedb("tcell_search", params)
    
    if df.empty:
        print(f"  ⚠️ No records returned for {config['label']}.")
        return df
    
    print(f"\n  Raw records: {len(df):,}")
    
    # --------------------------------------------------
    # Post-query filters
    # --------------------------------------------------
    
    if "linear_sequence_length" in df.columns:
        df["peptide_length"] = df["linear_sequence_length"]
    else:
        df["peptide_length"] = df["linear_sequence"].str.len()
    
    df = df[
        (df["peptide_length"] >= config["min_length"]) & 
        (df["peptide_length"] <= config["max_length"])
    ]
    print(f"  After length filter ({config['min_length']}-{config['max_length']} aa): {len(df):,}")
    
    def is_non_human_source(org_iri):
        if pd.isna(org_iri):
            return False
        return str(org_iri) != HUMAN_IRI
    
    df = df[df["parent_source_antigen_source_org_iri"].apply(is_non_human_source)]
    print(f"  After excluding human-source epitopes: {len(df):,}")
    
    df = df[df["mhc_allele_name"].notna() & (df["mhc_allele_name"] != "")]
    print(f"  After requiring MHC allele: {len(df):,}")
    
    df["has_receptor"] = df["receptor_ids"].notna()
    df["mhc_class_label"] = config["label"]
    
    print(f"  Records with TCR data: {df['has_receptor'].sum():,}")
    print(f"  Unique epitope sequences: {df['linear_sequence'].nunique():,}")
    print(f"  Unique source organisms: {df['parent_source_antigen_source_org_name'].nunique():,}")
    print(f"  Unique MHC alleles: {df['mhc_allele_name'].nunique():,}")
    
    # Save full assay data
    out_path = os.path.join(OUT_DIR, mhc_class_key, f"iedb_tcell_{mhc_class_key}_epitopes.csv")
    df.to_csv(out_path, index=False)
    print(f"\n  ✅ T cell data saved to {out_path}")
    
    # Save deduplicated epitope-level file for downstream BLAST
    dedup = df.drop_duplicates(subset=["linear_sequence", "mhc_allele_name"])
    dedup_path = os.path.join(OUT_DIR, mhc_class_key, f"iedb_unique_epitopes_{mhc_class_key}.csv")
    dedup.to_csv(dedup_path, index=False)
    print(f"  ✅ Unique epitope-MHC pairs: {len(dedup):,} → {dedup_path}")
    
    return df


# ====================================================
# Summary
# ====================================================
def print_summary(tcell_df, mhc_class_key):
    config = MHC_CLASS_CONFIG[mhc_class_key]
    
    print(f"\n{'='*60}")
    print(f"Summary: {config['label']}")
    print(f"{'='*60}")
    
    print(f"  Total T cell assay records: {len(tcell_df):,}")
    print(f"  Unique epitope sequences: {tcell_df['linear_sequence'].nunique():,}")
    print(f"  Unique source organisms: {tcell_df['parent_source_antigen_source_org_name'].nunique():,}")
    print(f"  Unique MHC alleles: {tcell_df['mhc_allele_name'].nunique():,}")
    print(f"  Records with receptor_ids: {tcell_df['has_receptor'].sum():,}")
    
    print(f"\n  Top 10 source organisms:")
    orgs = tcell_df["parent_source_antigen_source_org_name"].value_counts().head(10)
    for org, count in orgs.items():
        print(f"    {org}: {count}")
    
    print(f"\n  Top 10 MHC alleles:")
    alleles = tcell_df["mhc_allele_name"].value_counts().head(10)
    for allele, count in alleles.items():
        print(f"    {allele}: {count}")
    
    print(f"\n  Peptide length distribution:")
    print(tcell_df["peptide_length"].value_counts().sort_index().to_string())


# ====================================================
# Main
# ====================================================
if __name__ == "__main__":
    print("=" * 60)
    print("IEDB IQ-API Query Pipeline")
    print("Microbial T cell epitopes (MHC Class I and II)")
    print("=" * 60)
    
    results = {}
    
    for mhc_class_key in ["mhc_i", "mhc_ii"]:
        config = MHC_CLASS_CONFIG[mhc_class_key]
        print(f"\n\n{'#'*60}")
        print(f"# Processing: {config['label']}")
        print(f"{'#'*60}")
        
        tcell_df = fetch_tcell_epitopes(mhc_class_key)
        
        if tcell_df.empty:
            print(f"\n  ❌ No data for {config['label']}. Skipping.")
            results[mhc_class_key] = pd.DataFrame()
            continue
        
        results[mhc_class_key] = tcell_df
    
    # --------------------------------------------------
    # Print summaries
    # --------------------------------------------------
    print(f"\n\n{'#'*60}")
    print(f"# FINAL SUMMARIES")
    print(f"{'#'*60}")
    
    for mhc_class_key in ["mhc_i", "mhc_ii"]:
        if not results[mhc_class_key].empty:
            print_summary(results[mhc_class_key], mhc_class_key)
    
    # Combined
    print(f"\n{'='*60}")
    print("Combined Statistics")
    print(f"{'='*60}")
    
    total_records = sum(len(df) for df in results.values() if not df.empty)
    total_epitopes = sum(df["linear_sequence"].nunique() for df in results.values() if not df.empty)
    total_with_receptor = sum(df["has_receptor"].sum() for df in results.values() if not df.empty)
    
    print(f"  Total T cell assay records: {total_records:,}")
    print(f"  Total unique epitope sequences: {total_epitopes:,}")
    print(f"  Total records with receptor_ids: {total_with_receptor:,}")
    
    print(f"\n  Output directory: {OUT_DIR}")
    print(f"    mhc_i/  — MHC Class I epitopes")
    print(f"    mhc_ii/ — MHC Class II epitopes")
    
    print(f"\n  Next steps:")
    print(f"    1. Download IEDB receptor_full_v3.csv bulk export for TCR linkage")
    print(f"    2. BLAST epitope sequences against human proteome")
    print(f"    3. Cross-reference with VDJdb for additional TCR sequences")
    print(f"    4. Run NetMHCpan / NetMHCIIpan on mimicry pairs")
    print(f"    5. Feed into DecoderTCR embedding analysis")
    
    print(f"\n✅ Pipeline complete.")