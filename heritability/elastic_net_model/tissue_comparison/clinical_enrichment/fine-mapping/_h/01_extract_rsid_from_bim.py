#!/usr/bin/env python3

import sys

if len(sys.argv) != 3:
    print("Usage: python 01_extract_rsid_from_bim.py <input.bim> <output_rsid_list.txt>")
    sys.exit(1)

bim_file = sys.argv[1]
output_file = sys.argv[2]

with open(bim_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        cols = line.strip().split()
        snp_id = cols[1]

        # Extract rsID (assumes rsID is last underscore-separated field)
        rsid = snp_id.split("_")[-1]

        outfile.write(rsid + "\n")

print(f"Extracted rsIDs written to {output_file}")