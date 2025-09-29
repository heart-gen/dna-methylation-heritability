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
def get_enet():
    # get vmrs
    enet_fn = here("elastic_net_model/caudate/_m/caudate_betas_elastic-net.tsv.gz")
    df = pd.read_csv(enet_fn, sep='\t')

    # map vmr ids
    bed_fn = here("meqtl-analysis/vmrs/caudate/_m/feature.bed")
    bed = pd.read_csv(bed_fn, sep="\t", usecols=[0, 1, 2, 3], header=0)
    bed.columns = ['phenotype_id', 'chr', 'start', 'stop']
    merged_enet = pd.merge(df, bed, on=['chr', 'start', 'stop'], how='left')

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
def get_meqtl():
    fn = "../../_m/TOPMed_LIBD.permutation.txt.gz"
    return pd.read_csv(fn, sep='\t')

@lru_cache()
def merge_dataframe():
    return get_enet().merge(get_meqtl(), left_on='phenotype_id', how='inner')


def cal_fishers_direction(direction, h2_cat):
    df = merge_dataframe()
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
    h2_lt = []; dir_lt = []; fdr_lt = []; pval_lt = []; oddratio_lt = []
    for h2_cat in ["Heritable", "Non-heritable", "Low prediction"]:
        pvals = []
        for direction in ['Up', 'Down', 'All']:
                odd_ratio, pval = cal_fishers_direction(direction, h2_cat)
                pvals.append(pval); h2_lt.append(h2_cat)
                oddratio_lt.append(odd_ratio); dir_lt.append(direction)
        _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
        pval_lt = np.concatenate((pval_lt, pvals))
        fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'h2_Category': h2_lt, 'OR': oddratio_lt,
                         'PValue': pval_lt, "FDR": fdr_lt, 'Direction': dir_lt})


def main():
    calculate_enrichment().to_csv('meqtl_vmr_enrichment_analysis.txt', sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()