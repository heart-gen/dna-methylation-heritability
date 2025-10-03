"""
This script performs the Fisher's pairwise enrichment analysis
between significant meQTL and heritable sites estimated by
elastic-net.
"""
import numpy as np
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection
import session_info

@lru_cache()
def get_enet(tissue):
    # get vmrs
    enet_fn = here(f"heritability/elastic_net_model/{tissue.lower()}/_m/{tissue.lower()}_summary_elastic-net.tsv")
    df = pd.read_csv(enet_fn, sep='\t')
    df['chrom'] = 'chr' + df['chrom'].astype(str)

    # map vmr ids
    bed_fn = here(f"meqtl-analysis/vmrs/{tissue.lower()}/_m/feature.bed")
    bed = pd.read_csv(bed_fn, sep="\t", usecols=[0, 1, 2, 3], header=0)
    bed.columns = ['phenotype_id', 'chrom', 'start', 'end']
    merged_enet = pd.merge(df, bed, on=['chrom', 'start', 'end'], how='left')

    # assign h2 categories
    merged_enet['h2_category'] = merged_enet.apply(
    lambda row: (
        'Heritable' if row['h2_unscaled'] >= 0.1 and row['r_squared_cv'] > 0.75 else
        'Non-heritable' if row['h2_unscaled'] < 0.1 and row['r_squared_cv'] > 0.75 else
        'Low prediction'
    ),
    axis=1
    )
    return merged_enet

@lru_cache()
def get_meqtl(tissue):
    fn = here(f"meqtl-analysis/vmrs/{tissue.lower()}/cis_analysis/_m/TOPMed_LIBD.permutation.txt.gz")
    return pd.read_csv(fn, sep='\t')

@lru_cache()
def merge_dataframe(tissue):
    return get_enet(tissue).merge(get_meqtl(tissue), 
                                  on='phenotype_id', how='inner')


def cal_fishers_direction(direction, h2_cat, tissue):
    df = merge_dataframe(tissue)
    if direction == 'Up':
        df = df[(df['slope'] > 0)].copy()
    elif direction == 'Down':
        df = df[(df['slope'] < 0)].copy()
    else:
        df = df
    
    table = [[np.sum((df['qval']<0.05) & (df['h2_category'] == h2_cat)), 
              np.sum((df['qval']<0.05) & (df['h2_category'] != h2_cat))],
             [np.sum((df['qval']>0.05) & (df['h2_category'] == h2_cat)), 
              np.sum((df['qval']>0.05) & (df['h2_category'] != h2_cat))]]
    print(table)
    return fisher_exact(table)


def calculate_enrichment():
    region_lt = []; h2_lt = []; dir_lt = []; fdr_lt = []; pval_lt = []; oddratio_lt = []
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        for h2_cat in ["Heritable", "Non-heritable", "Low prediction"]:
            pvals = []
            for direction in ['Up', 'Down', 'All']:
                    odd_ratio, pval = cal_fishers_direction(direction, h2_cat, tissue)
                    pvals.append(pval); h2_lt.append(h2_cat)
                    oddratio_lt.append(odd_ratio); dir_lt.append(direction)
                    region_lt.append(tissue)
            _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
            pval_lt = np.concatenate((pval_lt, pvals))
            fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Tissue': region_lt, 'h2_Category': h2_lt, 
                         'OR': oddratio_lt, 'PValue': pval_lt, 
                         "FDR": fdr_lt, 'Direction': dir_lt})


def main():
    calculate_enrichment().to_csv('meqtl_vmr_enrichment_analysis.txt', sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()