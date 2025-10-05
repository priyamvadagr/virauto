import os

viral_dir = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1_viral_9_mers"
human_dir = "/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/9_mers/type1_human_9_mers"

# List all result files
viral_files = [f for f in os.listdir(viral_dir) if f.endswith(".txt")]
human_files = [f for f in os.listdir(human_dir) if f.endswith(".txt")]

# Extract allele names (remove chunk info + extension)
def extract_allele(filename):
    base = os.path.splitext(filename)[0]
    # Remove trailing _chunk### if present
    return base.split("_chunk")[0]

viral_alleles = {extract_allele(f) for f in viral_files}
human_alleles = {extract_allele(f) for f in human_files}

# Find intersection (alleles done for both viral and human)
common_alleles = sorted(viral_alleles & human_alleles)
only_viral = sorted(viral_alleles - human_alleles)
only_human = sorted(human_alleles - viral_alleles)

print(f"✅ Alleles processed in BOTH viral & human: {len(common_alleles)}")
print(f"⚠️ Alleles only in viral: {len(only_viral)}")
print(f"⚠️ Alleles only in human: {len(only_human)}")

# Write to files for convenience
outdir = "/ix/djishnu/Priyamvada/virauto/data/epitopes/virscan/paired_k_mers/9_mers"
with open(os.path.join(outdir, "alleles_in_both.txt"), "w") as f:
    f.write("\n".join(common_alleles))
with open(os.path.join(outdir, "alleles_only_viral.txt"), "w") as f:
    f.write("\n".join(only_viral))
with open(os.path.join(outdir, "alleles_only_human.txt"), "w") as f:
    f.write("\n".join(only_human))

print("\nSaved results to:")
print(f"  - alleles_in_both.txt")
print(f"  - alleles_only_viral.txt")
print(f"  - alleles_only_human.txt")
