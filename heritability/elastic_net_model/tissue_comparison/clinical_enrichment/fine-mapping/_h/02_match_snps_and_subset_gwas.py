#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) != 5:
    print("Usage: python 02_match_snps_and_subset_gwas.py <rsid_list.txt> <gwas.tsv> <output_snp_list.txt> <output_gwas.tsv>")
    sys.exit(1)

rsid_file = sys.argv[1]
gwas_file = sys.argv[2]
output_snp_file = sys.argv[3]
output_gwas_file = sys.argv[4]

# Load SNP list (preserve order)
with open(rsid_file) as f:
    snps = [line.strip() for line in f]

# Load GWAS
gwas = pd.read_csv(gwas_file, sep="\t", dtype=str)

if "rsid" not in gwas.columns:
    raise ValueError("GWAS file must contain a column named 'rsid'")

# Keep only SNPs present in GWAS
gwas_filtered = gwas[gwas["rsid"].isin(snps)].copy()

# Reorder GWAS to match SNP list order
gwas_filtered["order"] = gwas_filtered["rsid"].map({s:i for i,s in enumerate(snps)})
gwas_filtered = gwas_filtered.sort_values("order")
gwas_filtered = gwas_filtered.drop(columns=["order"])

# Final SNP list (perfect match and order)
matched_snps = gwas_filtered["rsid"].tolist()

# Write ordered SNP list
with open(output_snp_file, "w") as f:
    for snp in matched_snps:
        f.write(snp + "\n")

# Write subset GWAS
gwas_filtered.to_csv(output_gwas_file, sep="\t", index=False)

print(f"Matched {len(matched_snps)} SNPs")
print(f"Filtered SNP list: {output_snp_file}")
print(f"Subset GWAS: {output_gwas_file}")