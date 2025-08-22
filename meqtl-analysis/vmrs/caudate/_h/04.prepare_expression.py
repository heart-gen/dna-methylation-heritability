#!/usr/bin/env python
# Author: Francois Aguet
## Edited by Kynon J Benjamin (05/05/2021)

import argparse
import gzip, os
import subprocess
import numpy as np
import pandas as pd
import scipy.stats as stats


def read_gct(gct_file, sample_ids=None, dtype=None):
    """
    Load GCT as DataFrame. The first two columns must be 'Name' and 'Description'.
    """
    if sample_ids is not None:
        sample_ids = ['Name']+list(sample_ids)
    if gct_file.endswith('.gct.gz') or gct_file.endswith('.gct'):
        if dtype is not None:
            with gzip.open(gct_file, 'rt') as gct:
                gct.readline()
                gct.readline()
                sample_ids = gct.readline().strip().split()
            dtypes = {i:dtype for i in sample_ids[2:]}
            dtypes['Name'] = str
            dtypes['Description'] = str
            df = pd.read_csv(gct_file, sep='\t', skiprows=2, usecols=sample_ids,
                             index_col=0, dtype=dtypes)
        else:
            df = pd.read_csv(gct_file, sep='\t', skiprows=2, usecols=sample_ids,
                             index_col=0)
    elif gct_file.endswith('.parquet'):
        df = pd.read_parquet(gct_file, columns=sample_ids)
    elif gct_file.endswith('.ft'):  # feather format
        df = pd.read_feather(gct_file, columns=sample_ids)
        df = df.set_index('Name')
    else:
        raise ValueError('Unsupported input format.')
    df.index.name = 'feature_id'
    if 'Description' in df.columns:
        df = df.drop('Description', axis=1)
    return df


def sort_samples(norm_df, vcf_lookup_s, sort=True):
    """
    vcf_lookup: lookup table mapping sample IDs to VCF IDs
    """
    if sort:
        ix = np.intersect1d(norm_df.columns, vcf_lookup_s.index)
        norm_df = norm_df[ix]
    return norm_df


def get_bed(bed_file, FLIP):
    bed_template_df = pd.read_csv(bed_file, sep='\t', index_col=0)\
                        .rename(columns={'seqnames':'chr'})\
                        .loc[:, ["chr", "start", "end", "feature_id", "strand"]]
    if FLIP:
        # fix strand TSS assignment
        gene_order = bed_template_df.index
        bed_df_minus = bed_template_df[(bed_template_df["strand"] == "-")].copy()
        bed_df_minus.loc[:, "start"] = bed_df_minus.end
        bed_df_minus.loc[:, "end"] = bed_df_minus.end + 1
        bed_df_plus = bed_template_df[(bed_template_df["strand"] == "+")].copy()
        bed_df_plus.loc[:, "end"] = bed_df_plus.start + 1
        bed_template_df = pd.concat([bed_df_plus, bed_df_minus], axis=0)\
                            .drop(["strand"], axis=1)\
                            .loc[gene_order, :]
    else:
        bed_template_df.loc[:, "end"] = bed_template_df.start + 1
    return bed_template_df


def prepare_bed(df, bed_template_df, chr_subset=None):
    bed_df = pd.merge(bed_template_df, df, left_on="feature_id",
                      right_index=True)
    # sort by start position
    bed_df = bed_df.groupby('chr', sort=False, group_keys=False)\
                   .apply(lambda x: x.sort_values('start'),
                          include_groups=True)
    if chr_subset is not None:
        # subset chrs from VCF
        bed_df = bed_df[bed_df.chr.isin(chr_subset)]
    return bed_df


def write_bed(bed_df, output_name):
    """
    Write DataFrame to BED format
    """
    assert bed_df.columns[0]=='chr' and bed_df.columns[1]=='start' and bed_df.columns[2]=='end'
    # header must be commented in BED format
    header = bed_df.columns.values.copy()
    header[0] = '#chr'
    bed_df.to_csv(output_name, sep='\t', index=False, header=header)
    subprocess.check_call(f'bgzip -f {output_name}', shell=True,
                          executable='/bin/bash')
    subprocess.check_call(f'tabix -f {output_name}.gz', shell=True,
                          executable='/bin/bash')


def main():
    parser = argparse.ArgumentParser(description='Generate normalized expression BED files for eQTL analyses')
    parser.add_argument('norm_gct',
                        help='GCT file with normalized expression')
    parser.add_argument('sample_participant_lookup',
                        help='Lookup table linking samples to participants')
    parser.add_argument('vcf_chr_list', help='List of chromosomes in VCF')
    parser.add_argument('prefix', help='Prefix for output file names')
    parser.add_argument('-o', '--output_dir', default='.',
                        help='Output directory')
    parser.add_argument('--sample_id_list', default=None,
                        help='File listing sample IDs to include')
    parser.add_argument('--feature', default='gene', help='gene, transcript or exon')
    parser.add_argument('--bed_file', help='this is the bed file annotation')
    parser.add_argument('--flip', action="store_false",
                        help='Flip TSS using strand information')
    args = parser.parse_args()

    print('Loading expression data', flush=True)
    sample_ids = None

    if args.sample_id_list is not None:
        with open(args.sample_id_list) as f:
            sample_ids = f.read().strip().split('\n')
            print(f'  * Loading {len(sample_ids)} samples', flush=True)

    norm_df = read_gct(args.norm_gct, sample_ids)
    print(f'Map data', flush=True)
    sample_participant_lookup_s = pd.read_csv(args.sample_participant_lookup,
                                              sep='\t', index_col=0, dtype=str)\
                                    .squeeze("columns")
    norm_df = sort_samples(norm_df, sample_participant_lookup_s, sort=False)
    print(f'  * {norm_df.shape[0]} genes.', flush=True)

    # change sample IDs to participant IDs
    norm_df.rename(columns=sample_participant_lookup_s.to_dict(), inplace=True)
    bed_template_df = get_bed(args.bed_file, args.flip)
    print('bed_template_df.shape', bed_template_df.shape, flush=True)
    with open(args.vcf_chr_list) as f:
        chr_list = f.read().strip().split('\n')
    norm_bed_df = prepare_bed(norm_df, bed_template_df, chr_subset=chr_list)
    print(f'  * {norm_bed_df.shape[0]} genes remain after removing contigs absent from VCF.', flush=True)
    print('Writing BED file', flush=True)
    write_bed(norm_bed_df, os.path.join(args.output_dir,
                                        args.prefix+'.expression.bed'))


if __name__=='__main__':
    main()
