import pandas as pd

HLA_RISK_FILE = "/ix/djishnu/Priyamvada/pepSource/data/Ogishi_biorxiv_2019/HLA_autoimmunity_classification.txt"
OUTPUT_FILE = "/ix/djishnu/Priyamvada/pepSource/data/Ogishi_biorxiv_2019/HLA_autoimmunity_classification_standard.tsv"

# parse
risk = pd.read_csv(HLA_RISK_FILE, sep=r"\s+", engine="python")
risk.columns = risk.columns.str.strip('"')
for col in risk.select_dtypes(include="object").columns:
    risk[col] = risk[col].str.strip('"')

print(f"Loaded: {risk.shape[0]} rows")
print(f"Columns: {list(risk.columns)}")
print(f"\nOriginal HLA examples:")
print(risk["HLA"].head(10).tolist())


def convert_hla_name(name):
    """HLA_DRB1_0301 -> HLA-DRB1*03:01, HLA_B_2705 -> HLA-B*27:05"""
    parts = name.split("_")
    if len(parts) != 3:
        return name
    gene = parts[1]
    digits = parts[2]
    if len(digits) == 4:
        return f"HLA-{gene}*{digits[:2]}:{digits[2:]}"
    elif len(digits) == 2:
        return f"HLA-{gene}*{digits}"
    return name


risk["HLA"] = risk["HLA"].apply(convert_hla_name)

print(f"\nConverted HLA examples:")
print(risk["HLA"].head(10).tolist())
print(f"\nUnique alleles: {risk['HLA'].nunique()}")

risk.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"\nSaved to: {OUTPUT_FILE}")