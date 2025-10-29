import pandas as pd
import os
pair_summary = pd.read_csv("/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/virscan/Type1_NR_all/delta_binding/9_mers/summaries/pair_binding_summary_HLA_A.csv")
import matplotlib.pyplot as plt
import seaborn as sns
def classify_pair(row, thr=0.1):
    if row["Human_dominant_fraction"] >= thr:
        return "Human_dominant"
    elif row["Viral_dominant_fraction"] >= thr:
        return "Viral_dominant"
    else:
        return "Mixed"

pair_summary["Pair_Category"] = pair_summary.apply(classify_pair, axis=1)

pair_summary[pair_summary["Pair_Category"] == "Human_dominant"]
