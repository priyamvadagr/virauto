import pandas as pd
import argparse

HLA_RISK_FILE = "/ix/djishnu/Priyamvada/pepSource/data/Ogishi_biorxiv_2019/HLA_autoimmunity_classification.txt"


def risk_to_standard(hla_str):
    """Convert e.g. HLA_DRB1_0301 -> HLA-DRB1*03:01, HLA_B_2705 -> HLA-B*27:05"""
    parts = hla_str.split("_")
    if len(parts) != 3:
        return None
    gene = parts[1]
    digits = parts[2]
    if len(digits) == 4:
        return f"HLA-{gene}*{digits[:2]}:{digits[2:]}"
    elif len(digits) == 2:
        return f"HLA-{gene}*{digits}"
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Check overlap between HLA risk annotations and prediction alleles"
    )
    parser.add_argument(
        "--blast-input",
        required=True,
        help="Path to BLAST input CSV (e.g. iedb_mhc_i_pairs_4digit_hla.csv.gz)",
    )
    parser.add_argument(
        "--mhc-class",
        choices=["I", "II"],
        required=True,
        help="MHC class to check (I or II)",
    )
    parser.add_argument(
        "--hla-col",
        required=True,
        help="HLA column name in the BLAST input file",
    )
    args = parser.parse_args()

    class_key = "ClassI" if args.mhc_class == "I" else "ClassII"

    # ── parse risk file ──
    risk = pd.read_csv(HLA_RISK_FILE, sep=r'\s+', quotechar='"')
    risk["hla_standard"] = risk["HLA"].apply(risk_to_standard)

    print(f"=== MHC Class {args.mhc_class} | HLA column: {args.hla_col} ===\n")

    # filter to relevant class
    risk_cls = risk[risk["HLA_Class"] == class_key].copy()
    risk_std = set(risk_cls["hla_standard"].dropna().unique())
    print(f"Risk alleles ({class_key}): {len(risk_cls)} rows, "
          f"{len(risk_std)} unique standardized")
    print(f"  Examples: {sorted(risk_std)[:10]}")

    # ── load prediction alleles ──
    blast = pd.read_csv(args.blast_input)
    if args.hla_col not in blast.columns:
        print(f"\nERROR: column '{args.hla_col}' not found.")
        print(f"Available columns: {list(blast.columns)}")
        return
    pred_alleles = blast[args.hla_col].dropna().unique()
    print(f"\nPrediction alleles ({args.hla_col}): {len(pred_alleles)}")
    print(f"  Examples: {sorted(pred_alleles)[:10]}")

    # ── exact match ──
    exact = risk_std & set(pred_alleles)
    print(f"\nExact matches: {len(exact)}")
    for a in sorted(exact):
        assoc = risk_cls.loc[risk_cls["hla_standard"] == a, "Association"].values
        print(f"  {a:30s}  ({', '.join(assoc[:3])})")

    # ── substring match ──
    substring_matches = []
    for risk_a in sorted(risk_std):
        short = risk_a.replace("HLA-", "")
        for pred_a in pred_alleles:
            if short in pred_a and pred_a not in exact:
                substring_matches.append((risk_a, pred_a))

    print(f"\nSubstring matches (chain in paired allele): {len(substring_matches)}")
    for risk_a, pred_a in substring_matches[:25]:
        print(f"  risk: {risk_a:30s}  matches pred: {pred_a}")
    if len(substring_matches) > 25:
        print(f"  ... and {len(substring_matches) - 25} more")

    # ── unmatched ──
    matched = exact | {r for r, _ in substring_matches}
    unmatched = risk_std - matched
    print(f"\nUnmatched risk alleles: {len(unmatched)}")
    for a in sorted(unmatched)[:15]:
        print(f"  {a}")
    if len(unmatched) > 15:
        print(f"  ... and {len(unmatched) - 15} more")


if __name__ == "__main__":
    main()