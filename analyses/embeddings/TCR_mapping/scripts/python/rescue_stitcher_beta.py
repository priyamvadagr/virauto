#!/usr/bin/env python3
"""
======================================================================
Script: rescue_stitchr_beta.py

Description:
    Rescue failed beta chain stitchr reconstructions caused by
    ambiguous TRBV family names (e.g. TRBV6 instead of TRBV6-5).

    Strategy:
        1. Identify failed rows from stitchr_output_beta.tsv
        2. Apply most-common-subgroup fallback mapping to TRBV column
        3. Write a patched input TSV for the failed rows only
        4. Run thimble on the patched input
        5. Merge rescued rows back with original successes

    NOTE: Subgroup fallback introduces a small inaccuracy — the
    reconstructed V gene germline sequence may not exactly match the
    original clonotype. The CDR3 sequence is preserved exactly.
    Affected TCR_names are logged to rescue_fallback_log.tsv for
    downstream sensitivity analysis.

Output:
    stitchr_input_beta_rescue.tsv      → patched thimble input
    stitchr_output_beta_rescue.tsv     → thimble output for rescued rows
    stitchr_output_beta_merged.tsv     → final merged beta output
    rescue_fallback_log.tsv            → log of all fallback substitutions

Usage:
    python rescue_stitchr_beta.py
    # Then run thimble manually (printed at end of script), then
    python rescue_stitchr_beta.py --merge
======================================================================
"""

import os
import argparse
import pandas as pd

# ====================================================================
# Config
# ====================================================================
STITCHR_DIR  = "/ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/mimicry_candidates/mhc_i/stitchr"

BETA_INPUT   = os.path.join(STITCHR_DIR, "stitchr_input_beta.tsv")
BETA_OUTPUT  = os.path.join(STITCHR_DIR, "stitchr_output_beta.tsv")
RESCUE_INPUT = os.path.join(STITCHR_DIR, "stitchr_input_beta_rescue.tsv")
RESCUE_OUT   = os.path.join(STITCHR_DIR, "stitchr_output_beta_rescue.tsv")
MERGED_OUT   = os.path.join(STITCHR_DIR, "stitchr_output_beta_merged.tsv")
FALLBACK_LOG = os.path.join(STITCHR_DIR, "rescue_fallback_log.tsv")

# Most common subgroup per ambiguous TRBV family
# Based on healthy human peripheral blood repertoire frequency data.
# Applied ONLY when the family name has no subgroup suffix (e.g. TRBV6)
# and the original stitchr run produced an Error.
# ⚠️  NOTE: These are population-level approximations. The true subgroup
#     for any individual clonotype is unknown without full sequence data.
TRBV_FALLBACK = {
    "TRBV6"  : "TRBV6-5",
    "TRBV7"  : "TRBV7-2",
    "TRBV5"  : "TRBV5-1",
    "TRBV20" : "TRBV20-1",
    "TRBV4"  : "TRBV4-1",
    "TRBV10" : "TRBV10-3",
    "TRBV3"  : "TRBV3-1",
    "TRBV29" : "TRBV29-1",
    "TRBV11" : "TRBV11-2",
    "TRBV12" : "TRBV12-3",
    "TRBV24" : "TRBV24-1",
    "TRBV25" : "TRBV25-1",
}

# Allele-level failures — strip allele suffix and retry with *01 default
# e.g. TRBV19*01 → TRBV19  (stitchr will default to *01 anyway but
# the explicit *01 on a non-*01 allele caused the lookup to fail)
STRIP_ALLELE_PATTERN = r"\*\d+$"


def prepare_rescue():
    """Identify failed rows, apply fallbacks, write rescue input TSV."""
    print("Loading stitchr beta output...")
    beta_out = pd.read_csv(BETA_OUTPUT, sep="\t", dtype=str).fillna("")
    print(f"  Total rows: {len(beta_out):,}")

    # Failed rows = those with "Error:" in Warnings/Errors column
    err_col = "Warnings/Errors"
    failed  = beta_out[beta_out[err_col].str.contains("Error:", na=False)].copy()
    success = beta_out[~beta_out[err_col].str.contains("Error:", na=False)].copy()
    print(f"  Succeeded : {len(success):,}")
    print(f"  Failed    : {len(failed):,}")

    # Load original input to get the TRBV values for failed rows
    beta_in = pd.read_csv(BETA_INPUT, sep="\t", dtype=str).fillna("")
    failed_input = beta_in[beta_in["TCR_name"].isin(failed["TCR_name"])].copy()
    print(f"\n  Failed input rows to rescue: {len(failed_input):,}")

    # Apply fallback mapping
    fallback_records = []
    import re

    def apply_fallback(trbv):
        original = trbv
        # Strip allele suffix first (e.g. TRBV19*01 → TRBV19)
        trbv_base = re.sub(STRIP_ALLELE_PATTERN, "", trbv)
        if trbv_base in TRBV_FALLBACK:
            substituted = TRBV_FALLBACK[trbv_base]
            fallback_type = "subgroup_fallback"
        elif trbv != trbv_base:
            # Had allele suffix but not in fallback map — use base gene
            # and let stitchr default to *01
            substituted = trbv_base
            fallback_type = "allele_strip"
        else:
            substituted = trbv
            fallback_type = "unchanged"
        return substituted, original, fallback_type

    patched_trbv      = []
    original_trbv     = []
    fallback_types    = []

    for trbv in failed_input["TRBV"]:
        sub, orig, ftype = apply_fallback(trbv)
        patched_trbv.append(sub)
        original_trbv.append(orig)
        fallback_types.append(ftype)

    failed_input = failed_input.copy()
    failed_input["TRBV"] = patched_trbv

    # Log fallback substitutions
    log = failed_input[["TCR_name"]].copy()
    log["original_TRBV"]  = original_trbv
    log["fallback_TRBV"]  = patched_trbv
    log["fallback_type"]  = fallback_types
    log.to_csv(FALLBACK_LOG, sep="\t", index=False)
    print(f"\n  Fallback substitutions:")
    print(log["fallback_type"].value_counts().to_string())
    print(f"  Fallback log → {FALLBACK_LOG}")

    # Write rescue input TSV
    failed_input.to_csv(RESCUE_INPUT, sep="\t", index=False)
    print(f"\n  Rescue input → {RESCUE_INPUT}  ({len(failed_input):,} rows)")
    print(f"""
  Now run thimble on the rescue input:

    thimble \\
        -in {RESCUE_INPUT} \\
        -o  {RESCUE_OUT} \\
        -r b -s HUMAN

  Then rerun this script with --merge:
    python rescue_stitchr_beta.py --merge
""")


def merge_rescue():
    """Merge original successes with rescue output into final beta TSV."""
    print("Merging rescue output with original successes...")

    beta_out   = pd.read_csv(BETA_OUTPUT,  sep="\t", dtype=str).fillna("")
    rescue_out = pd.read_csv(RESCUE_OUT,   sep="\t", dtype=str).fillna("")

    err_col = "Warnings/Errors"
    success = beta_out[~beta_out[err_col].str.contains("Error:", na=False)].copy()

    # From rescue: keep only newly succeeded rows
    rescue_success = rescue_out[~rescue_out[err_col].str.contains("Error:", na=False)].copy()
    rescue_failed  = rescue_out[rescue_out[err_col].str.contains("Error:", na=False)].copy()

    print(f"  Original successes  : {len(success):,}")
    print(f"  Rescue successes    : {len(rescue_success):,}")
    print(f"  Rescue still failed : {len(rescue_failed):,}")
    if not rescue_failed.empty:
        print(f"  Still-failed TRBV   : "
              f"{rescue_failed['TRBV'].value_counts().to_dict()}")

    merged = pd.concat([success, rescue_success], ignore_index=True)
    merged.to_csv(MERGED_OUT, sep="\t", index=False)
    print(f"\n  Final merged beta output → {MERGED_OUT}  ({len(merged):,} rows)")
    print(f"\n  Summary:")
    print(f"    Total TCRs input        : {len(beta_out):,}")
    print(f"    Successfully stitched   : {len(merged):,}  "
          f"({100*len(merged)/len(beta_out):.1f}%)")
    print(f"    Unrecoverable failures  : "
          f"{len(beta_out) - len(merged):,}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--merge", action="store_true",
                        help="Merge rescue output with original successes")
    args = parser.parse_args()

    if args.merge:
        merge_rescue()
    else:
        prepare_rescue()


