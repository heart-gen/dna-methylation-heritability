#!/usr/bin/env python
# Author: Francois Aguet
## Edited by Kynon J Benjamin (05/05/2021)

import argparse
import os
import subprocess
import numpy as np
import pandas as pd


def read_vmr(vmr_file, sample_ids=None):
    """
    Load VMR methylation matrix (feature_id x samples).
    Assumes first column = feature_id, remaining columns = samples.
    """
    df = pd.read_csv(vmr_file, sep='\t', index_col=0)
    if sample_ids is not None:
        sample_ids = [s for s in sample_ids if s in df.columns]
        df = df.loc[:, sample_ids]
    df.index.name = 'feature_id'
    return df


def sort_samples(vmr_df, vcf_lookup_s, sort=True):
    """
    Ensure VMR matrix columns match VCF sample IDs.
    vcf_lookup_s: pandas.Series mapping VMR sample IDs -> VCF IDs
    """
    if sort:
        ix = np.intersect1d(vmr_df.columns, vcf_lookup_s.index)
        vmr_df = vmr_df[ix]
    return vmr_df

def get_bed(bed_file):
    """
    Load VMR BED annotation file.
    Expected columns: seqnames, start, end, feature_id
    """
    bed_template_df = (
        pd.read_csv(bed_file, sep='\t', index_col=0)
          .rename(columns={'seqnames': 'chr'})
          .loc[:, ["chr", "start", "end", "feature_id"]]
    )
    return bed_template_df

def prepare_bed(df, bed_template_df, chr_subset=None):
    bed_df = pd.merge(bed_template_df, df, left_on="feature_id",
                      right_index=True)
    # sort by start position within each chromosome
    bed_df = bed_df.groupby('chr', sort=False, group_keys=False)\
                   .apply(lambda x: x.sort_values('start'),
                          include_groups=True)
    if chr_subset is not None:
        bed_df = bed_df[bed_df.chr.isin(chr_subset)]
    return bed_df

def write_bed(bed_df, output_name):
    """
    Write DataFrame to BED format (bgzipped + tabix indexed).
    """
    assert bed_df.columns[0] == 'chr'
    assert bed_df.columns[1] == 'start'
    assert bed_df.columns[2] == 'end'
    # header must be commented in BED format
    header = bed_df.columns.values.copy()
    header[0] = '#chr'
    bed_df.to_csv(output_name, sep='\t', index=False, header=header)
    subprocess.check_call(f'bgzip -f {output_name}', shell=True,
                          executable='/bin/bash')
    subprocess.check_call(f'tabix -f {output_name}.gz', shell=True,
                          executable='/bin/bash')

def main():
    parser = argparse.ArgumentParser(
        description='Generate normalized methylation BED files for mQTL analyses'
    )
    parser.add_argument('vmr_matrix',
                        help='Tab-delimited file with VMR methylation values (feature_id x samples)')
    parser.add_argument('sample_participant_lookup',
                        help='Lookup table linking VMR sample IDs to participant IDs')
    parser.add_argument('vcf_chr_list',
                        help='File with list of chromosomes in VCF')
    parser.add_argument('prefix',
                        help='Prefix for output file names')
    parser.add_argument('--output_dir', '-o', default='.',
                        help='Output directory')
    parser.add_argument('--sample_id_list', default=None,
                        help='Optional file listing sample IDs to include')
    parser.add_argument('--bed_file',
                        help='BED file annotation with feature_id coordinates',
                        required=True)
    args = parser.parse_args()

    print('Loading VMR methylation data', flush=True)
    sample_ids = None
    if args.sample_id_list is not None:
        with open(args.sample_id_list) as f:
            sample_ids = f.read().strip().split('\n')
            print(f'  * Restricting to {len(sample_ids)} samples', flush=True)

    vmr_df = read_vmr(args.vmr_matrix, sample_ids)
    print(f'  * Loaded {vmr_df.shape[0]} features across {vmr_df.shape[1]} samples', flush=True)

    print('Mapping sample IDs to participants', flush=True)
    sample_participant_lookup_s = pd.read_csv(
        args.sample_participant_lookup,
        sep='\t', index_col=0, dtype=str
    ).squeeze("columns")

    vmr_df = sort_samples(vmr_df, sample_participant_lookup_s, sort=False)
    vmr_df.rename(columns=sample_participant_lookup_s.to_dict(), inplace=True)

    bed_template_df = get_bed(args.bed_file)
    print('  * Annotation loaded:', bed_template_df.shape, flush=True)

    with open(args.vcf_chr_list) as f:
        chr_list = f.read().strip().split('\n')

    vmr_bed_df = prepare_bed(vmr_df, bed_template_df, chr_subset=chr_list)
    print(f'  * {vmr_bed_df.shape[0]} features remain after restricting to VCF contigs', flush=True)

    out_file = os.path.join(args.output_dir, args.prefix + '.methylation.bed')
    print('Writing BED file:', out_file, flush=True)
    write_bed(vmr_bed_df, out_file)


if __name__ == '__main__':
    main()
