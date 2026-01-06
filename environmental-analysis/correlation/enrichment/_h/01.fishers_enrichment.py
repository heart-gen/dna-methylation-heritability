"""
This script performs Fisher's pairwise enrichment analysis
between environmental features and heritable sites estimated by
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
    enet_fn = here(f"heritability/elastic_net_model/{tissue.lower()}/_m_to_del/{tissue.lower()}_summary_elastic-net.tsv")
    df = pd.read_csv(enet_fn, sep='\t')
    df.dropna()
    df['chrom'] = 'chr' + df['chrom'].astype(str)

    # assign h2 categories
    df['h2_category'] = df.apply(
    lambda row: (
        'Heritable' if row['h2_unscaled'] >= 0.1 and row['r_squared_cv'] > 0.75 else
        'Non-heritable' if row['h2_unscaled'] < 0.1 and row['r_squared_cv'] > 0.75 else
        'Low prediction'
    ),
    axis=1
    )
    return df

@lru_cache()
def get_vmrs():
    fn = here(f"environmental-analysis/correlation/_m/smoking_logit.csv.gz")
    vmrs = pd.read_csv(fn, sep=',') 
    vmrs['test'] = 'logit'
    vmrs['sig'] = vmrs['p'] < 0.05
    return vmrs

@lru_cache()
def get_dmrs():
    fn = here(f"environmental-analysis/correlation/_m/smoking_dmr.csv.gz")
    dmrs = pd.read_csv(fn, sep=',')
    dmrs['test'] = 'dmr'
    dmrs['sig'] = dmrs['p.value.x2'] < 0.05
    return dmrs

@lru_cache()
def concat_dataframe():
    return pd.concat([get_vmrs(), get_dmrs()])

@lru_cache()
def merge_dataframe(tissue):
    env = concat_dataframe()
    env['chr'] = 'chr' + env['chr'].astype(str)

    df = env.pivot_table(values='sig', columns='test', fill_value=0,
                         index=['chr', 'start', 'end'])
    
    return df.merge(get_enet(tissue), 
                    left_on=['chr', 'start', 'end'],
                    right_on=['chrom', 'start', 'end'],
                    how='left')
@lru_cache()
def cal_fishers_annot(tissue, test, h2_cat):
    df = merge_dataframe(tissue)
    print(df.head())
    
    table = [[np.sum((df[test.lower()] == 0) & (df['h2_category'] == h2_cat)), 
              np.sum((df[test.lower()] == 0) & (df['h2_category'] != h2_cat))],
             [np.sum((df[test.lower()] == 1) & (df['h2_category'] == h2_cat)), 
              np.sum((df[test.lower()] == 1) & (df['h2_category'] != h2_cat))]]
    print(table)
    return fisher_exact(table)

def calculate_enrichment():
    region_lt = []; h2_lt = []; test_lt = []; fdr_lt = []; pval_lt = []; oddratio_lt = []; env_lt = []

    for tissue in ["Caudate"]:
        for h2_cat in ["Heritable", "Non-heritable", "Low prediction"]:
            pvals = []
            for test in ["Logit", "DMR"]:
                for env in ["Smoking"]:
                    odd_ratio, pval = cal_fishers_annot(tissue, test, h2_cat)
                    pvals.append(pval); h2_lt.append(h2_cat)
                    oddratio_lt.append(odd_ratio); test_lt.append(test)
                    region_lt.append(tissue); env_lt.append(env)
            _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
            pval_lt = np.concatenate((pval_lt, pvals))
            fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Tissue': region_lt, 'h2_Category': h2_lt, 
                         'OR': oddratio_lt, 'PValue': pval_lt, 
                         "FDR": fdr_lt, 'Test': test_lt, 'Env': env_lt})

def main():
    calculate_enrichment().to_csv('environmental_vmr_enrichment_analysis.txt', sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()