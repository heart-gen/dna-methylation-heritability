#!/usr/bin/env python3

import sys

if len(sys.argv) != 5:
    print("Usage: python 03_regenerate_plink_subset.py <input_prefix> <matched_rsid_list.txt> <output_prefix> <plink_path>")
    sys.exit(1)

input_prefix = sys.argv[1]
matched_rsid_file = sys.argv[2]
output_prefix = sys.argv[3]
plink_path = sys.argv[4]

bim_file = input_prefix + ".bim"

# Load matched rsIDs
with open(matched_rsid_file) as f:
    matched_rsids = set(line.strip() for line in f)

# Map rsID -> original SNP ID
extract_ids = []

with open(bim_file) as f:
    for line in f:
        cols = line.strip().split()
        original_id = cols[1]
        rsid = original_id.split("_")[-1]

        if rsid in matched_rsids:
            extract_ids.append(original_id)

# Write extract file
extract_file = output_prefix + "_extract.txt"

with open(extract_file, "w") as f:
    for snp in extract_ids:
        f.write(snp + "\n")

print(f"Writing PLINK extract file: {extract_file}")

# Run PLINK
import os
cmd = f"{plink_path} --bfile {input_prefix} --extract {extract_file} --make-bed --out {output_prefix}"
print("Running:", cmd)
os.system(cmd)

print("Done.")