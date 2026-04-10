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

def filter_sites(enet_df, pop):
    enet_df = enet_df.dropna()
    enet_df['chrom'] = 'chr' + enet_df['chrom'].astype(str)

    # assign h2 categories
    enet_df['h2_category'] = enet_df.apply(
    lambda row: (
        'Heritable' if row['h2_unscaled'] >= 0.1 and row['r_squared_cv'] > 0.75 else
        'Non-heritable' if row['h2_unscaled'] < 0.1 and row['r_squared_cv'] > 0.75 else
        'Low prediction'
    ),
    axis=1
    )

    # assign feature ids
    enet_df['feature_id'] = (
        enet_df['chrom'].astype(str) + "_" +
        enet_df['start'].astype(str) + "_" +
        enet_df['end'].astype(str)
    )

    print(enet_df[['chrom', 'start', 'end', 'h2_category']].head())
    return enet_df

@lru_cache()
def get_matching_enet(tissue):
    # get vmrs
    enet_AA_fn = here(f"heritability/elastic_net_model/all_individuals/{tissue.lower()}/_m/{tissue.lower()}_summary_elastic-net_AA.tsv")
    enet_EA_fn = here(f"heritability/elastic_net_model/all_individuals/{tissue.lower()}/_m/{tissue.lower()}_summary_elastic-net_EA.tsv")

    df_AA = pd.read_csv(enet_AA_fn, sep='\t')
    df_EA = pd.read_csv(enet_EA_fn, sep='\t')

    enet_AA = filter_sites(df_AA, 'AA')
    enet_EA = filter_sites(df_EA, 'EA')
    
    enet_combined = enet_AA.merge(enet_EA, on='feature_id', 
                                  suffixes=('_AA', '_EA'))
    
    print(f"[ENET] Combined shape (after merge): {enet_combined.shape}")
    
    enet_combined = enet_combined[
        enet_combined['h2_category_AA'] == enet_combined['h2_category_EA']
    ]

    enet_combined = enet_combined.rename(columns={
        'h2_category_AA': 'h2_category',
        'chrom_AA': 'chrom',
        'start_AA': 'start',
        'end_AA': 'end'})
    
    enet_combined = enet_combined.drop(columns=[
        'h2_category_EA', 'chrom_EA', 'start_EA', 
        'end_EA', 'race_EA', 'race_AA'])
 
    return enet_combined

@lru_cache()
def get_vmrs(tissue, env, cat=None):
    fn = here(f"environmental-analysis/all_individuals/{tissue.lower()}/high_r2/correlation/_m/{env.lower()}_logit.csv.gz")
    vmrs = pd.read_csv(fn, sep=',') 

    if cat is not None:
        vmrs = vmrs[vmrs['var'] == f"{env}{cat}"]

    vmrs['test'] = 'logit'
    vmrs['sig'] = vmrs['p'] < 0.05

    print(f"[VMR] sig counts:\n{vmrs['sig'].value_counts()}")
    return vmrs

@lru_cache()
def get_dmrs(tissue, env, cat=None):
    fn = here(f"environmental-analysis/all_individuals/{tissue.lower()}/high_r2/correlation/_m/{env.lower()}_dmr.csv.gz")
    dmrs = pd.read_csv(fn, sep=',')

    if cat is not None:
        dmrs = dmrs[dmrs['level'] == cat]
    
    dmrs['test'] = 'dmr'
    dmrs['sig'] = dmrs['p.value'] < 0.05

    print(f"[DMR] sig counts:\n{dmrs['sig'].value_counts()}")
    return dmrs

@lru_cache()
def concat_dataframe(tissue, env, cat=None):
    return pd.concat([get_vmrs(tissue, env, cat), 
                      get_dmrs(tissue, env)])

@lru_cache()
def merge_dataframe(tissue, env, cat=None):
    merged = concat_dataframe(tissue, env, cat)
    merged['chr'] = 'chr' + merged['chr'].astype(str)

    df = merged.pivot_table(values='sig', columns='test', fill_value=0,
                         index=['chr', 'start', 'end'])
    
    df["both"] = ((df['logit'] == 1) & (df['dmr'] == 1)).astype(int)
    
    return df.merge(get_matching_enet(tissue), 
                    left_on=['chr', 'start', 'end'],
                    right_on=['chrom', 'start', 'end'],
                    how='left')
@lru_cache()
def cal_fishers_annot(tissue, test, h2_cat, env, cat=None):
    df = merge_dataframe(tissue, env, cat)
    
    table = [[np.sum((df[test.lower()] == 1) & (df['h2_category'] == h2_cat)), 
              np.sum((df[test.lower()] == 1) & (df['h2_category'] != h2_cat))],
             [np.sum((df[test.lower()] == 0) & (df['h2_category'] == h2_cat)), 
              np.sum((df[test.lower()] == 0) & (df['h2_category'] != h2_cat))]]
    print(table)
    return fisher_exact(table)

def calculate_enrichment():
    region_lt = []; h2_lt = []; test_lt = []; fdr_lt = []; pval_lt = []; oddratio_lt = []; env_lt = []

    env_vars = ["smoking", "codeine", "morphine", "cocaine", "ethanol", 
                "antipsychotics","nicotine","amphetamines", "hx_sexual_abuse","hx_physical_abuse"]

    categorical_vars = {
        "education": ['less_than_hs', 'more_than_hs'],
        "marital_status": ['single', 'previously_married']
    }

    for tissue in ["Caudate", "Hippocampus", "DLPFC"]:
        for h2_cat in ["Heritable", "Non-heritable", "Low prediction"]:
            pvals = []
            for test in ["Logit", "DMR", "Both"]:
                for env in env_vars:
                    odd_ratio, pval = cal_fishers_annot(tissue, test, h2_cat, env, None)
                    pvals.append(pval); h2_lt.append(h2_cat)
                    oddratio_lt.append(odd_ratio); test_lt.append(test)
                    region_lt.append(tissue); env_lt.append(env)
                    
                for env, vars in categorical_vars.items():
                    for cat in vars:
                        odd_ratio, pval = cal_fishers_annot(tissue, test, h2_cat, env, cat)
                        pvals.append(pval); h2_lt.append(h2_cat)
                        oddratio_lt.append(odd_ratio); test_lt.append(test)
                        region_lt.append(tissue); env_lt.append(cat)
                        
            _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
            pval_lt = np.concatenate((pval_lt, pvals))
            fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Tissue': region_lt, 'h2_Category': h2_lt, 
                         'OR': oddratio_lt, 'PValue': pval_lt, 
                         "FDR": fdr_lt, 'Test': test_lt, 'Env': env_lt})

def main():
    calculate_enrichment().to_csv('vmr_enrichment_analysis.txt', 
                                  sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
