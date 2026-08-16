#!/usr/bin/env python
"""
Create HapMap3 SNP list with alleles (SNP, A1, A2) from 1000G plink bim files.

This creates a proper --merge-alleles file for LDSC munge_sumstats.py by:
1. Reading all bim files from 1000G_EUR_Phase3_plink
2. Filtering to only HapMap3 SNPs (from hm3_no_MHC.list.txt)
3. Outputting SNP, A1, A2 columns

Usage:
    python create_hm3_snplist.py
"""
import pandas as pd
import os
import sys

# Paths
RESOURCE_DIR = "/projects/b1213/resources/ldsc"
BIM_DIR = os.path.join(RESOURCE_DIR, "1000G_EUR_Phase3_plink")
HM3_SNPS_FILE = os.path.join(RESOURCE_DIR, "hm3_no_MHC.list.txt")
OUTPUT_FILE = os.path.join(RESOURCE_DIR, "w_hm3.snplist")


def main():
    print(f"Loading HapMap3 SNP list from {HM3_SNPS_FILE}...")
    hm3_snps = set(pd.read_csv(HM3_SNPS_FILE, header=None, names=['SNP'])['SNP'].values)
    print(f"  Loaded {len(hm3_snps):,} HapMap3 SNPs")

    # Read all bim files and combine
    print(f"\nReading bim files from {BIM_DIR}...")
    bim_dfs = []
    for chrom in range(1, 23):
        bim_file = os.path.join(BIM_DIR, f"1000G.EUR.QC.{chrom}.bim")
        if os.path.exists(bim_file):
            df = pd.read_csv(bim_file, sep='\t', header=None,
                           names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
            bim_dfs.append(df[['SNP', 'A1', 'A2']])
            print(f"  Chr {chrom}: {len(df):,} SNPs")
        else:
            print(f"  Warning: {bim_file} not found")

    all_bim = pd.concat(bim_dfs, ignore_index=True)
    print(f"\nTotal SNPs from bim files: {len(all_bim):,}")

    # Filter to HapMap3 SNPs
    hm3_bim = all_bim[all_bim['SNP'].isin(hm3_snps)].copy()
    print(f"HapMap3 SNPs with allele info: {len(hm3_bim):,}")

    # Remove duplicates (keep first)
    hm3_bim = hm3_bim.drop_duplicates(subset='SNP', keep='first')
    print(f"After removing duplicates: {len(hm3_bim):,}")

    # Report SNPs in HM3 list but not in bim files
    missing = hm3_snps - set(hm3_bim['SNP'].values)
    if missing:
        print(f"\nNote: {len(missing):,} HapMap3 SNPs not found in bim files")

    # Write output
    print(f"\nWriting {OUTPUT_FILE}...")
    hm3_bim.to_csv(OUTPUT_FILE, sep='\t', index=False)
    print(f"Done! Created {OUTPUT_FILE}")

    # Also create gzipped version
    print(f"Creating gzipped version...")
    hm3_bim.to_csv(OUTPUT_FILE + '.gz', sep='\t', index=False, compression='gzip')
    print(f"Done! Created {OUTPUT_FILE}.gz")

    # Show sample
    print(f"\nSample output:")
    print(hm3_bim.head(10).to_string(index=False))


if __name__ == '__main__':
    main()
