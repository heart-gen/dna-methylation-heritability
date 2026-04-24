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

BINARY_ENVS = [
    "smoking", "codeine", "morphine", "cocaine", "ethanol",
    "nicotine", "amphetamines"
]

CATEGORICAL_ENVS = {
    "education": ["less_than_hs", "more_than_hs"],
    "marital_status": ["single", "previously_married"],
}

@lru_cache()
def get_enet(tissue, n_quintiles=5):
    # get vmrs
    enet_fn = here(f"heritability/elastic_net_model/BA_only/{tissue.lower()}/_m/{tissue.lower()}_summary_elastic-net.tsv")
    df = pd.read_csv(enet_fn, sep='\t')
    df = df.dropna()
    df['chrom'] = 'chr' + df['chrom'].astype(str)

    # Remove low prediction sites
    df = df[df['r_squared_cv'] > 0.3]

    # Define h2 quartiles
    breaks = np.unique(np.quantile(df["h2_unscaled"],
                                   q=np.linspace(0, 1, n_quintiles + 1)))
    n_bins = len(breaks) -1 

    if n_bins < 1:
        raise ValueError("Insufficient unique h2 values to compute quantile bins")

    if n_bins < n_quintiles: 
        raise ValueError(f"Reduced h2 bins from {n_quintiles} to {n_bins} because quantile breaks were not unique")
    
    labels = [f"Q{i}" for i in range(1, n_bins + 1)]

    df["h2_quintile"] = pd.cut(
        df["h2_unscaled"],
        bins=breaks,
        labels=labels,
        include_lowest=True
    )
    
    df = df[df["h2_quintile"].notna()]

    return df

@lru_cache()
def get_vmrs(tissue, quintile, env, cat=None):
    fn = here(f"environmental-analysis/BA_only/{tissue.lower()}/correlation/_m/{quintile}/{env.lower()}_linear.csv.gz")
    vmrs = pd.read_csv(fn, sep=',') 

    if cat is not None:
        vmrs = vmrs[vmrs['var'] == f"{env}{cat}"]

    vmrs['test'] = 'lm'
    vmrs['sig'] = vmrs['p'] < 0.05
    return vmrs

@lru_cache()
def get_dmrs(tissue, quintile, env, cat=None):
    fn = here(f"environmental-analysis/BA_only/{tissue.lower()}/correlation/_m/{quintile}/{env.lower()}_dmr.csv.gz")
    dmrs = pd.read_csv(fn, sep=',')

    if cat is not None:
        dmrs = dmrs[dmrs['level'] == cat]
    
    dmrs['test'] = 'dmr'
    dmrs['sig'] = dmrs['p.value'] < 0.05
    return dmrs

@lru_cache()
def concat_dataframe(tissue, quintile, env, cat=None):
    return pd.concat([get_vmrs(tissue, quintile, env, cat), 
                      get_dmrs(tissue, quintile, env, cat)])

@lru_cache()
def merge_dataframe(tissue, quintile, env, cat=None):
    merged = concat_dataframe(tissue, quintile, env, cat)
    merged['chr'] = 'chr' + merged['chr'].astype(str)

    df = merged.pivot_table(values='sig', columns='test', fill_value=0,
                         index=['chr', 'start', 'end'])
    
    df["both"] = ((df['lm'] == 1) & (df['dmr'] == 1)).astype(int)
    
    return df.merge(get_enet(tissue), 
                    left_on=['chr', 'start', 'end'],
                    right_on=['chrom', 'start', 'end'],
                    how='left')
@lru_cache()
def cal_fishers_annot(tissue, test, quintile, env, cat=None):
    df = merge_dataframe(tissue, quintile, env, cat)
    
    table = [[np.sum((df[test.lower()] == 1) 
                     & (df['h2_quintile'] == quintile)), 
              np.sum((df[test.lower()] == 1) 
                     & (df['h2_quintile'] != quintile))],
             [np.sum((df[test.lower()] == 0) 
                     & (df['h2_quintile'] == quintile)), 
              np.sum((df[test.lower()] == 0) 
                     & (df['h2_quintile'] != quintile))]]
    print(table)
    return fisher_exact(table)

def calculate_enrichment():
    region_lt = []; test_lt = []; fdr_lt = []; pval_lt = []; oddratio_lt = []; env_lt = []; quintile_lt = []

    for quintile in ["Q1", "Q2", "Q3", "Q4", "Q5"]:
        for tissue in ["Caudate", "Hippocampus", "DLPFC"]:
            pvals = []
            for test in ["LM", "DMR", "Both"]:
                for env in BINARY_ENVS:
                    odd_ratio, pval = cal_fishers_annot(tissue, test, quintile, env, None)
                    pvals.append(pval); oddratio_lt.append(odd_ratio) 
                    test_lt.append(test); region_lt.append(tissue)
                    env_lt.append(env); quintile_lt.append(quintile)
                        
                for env, vars in CATEGORICAL_ENVS.items():
                    for cat in vars:
                        odd_ratio, pval = cal_fishers_annot(tissue, test, quintile, env, cat)
                        pvals.append(pval); oddratio_lt.append(odd_ratio)
                        test_lt.append(test); region_lt.append(tissue); 
                        env_lt.append(cat ); quintile_lt.append(quintile)
                            
            _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
            pval_lt = np.concatenate((pval_lt, pvals))
            fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Tissue': region_lt, "Quintile": quintile_lt,
                         'OR': oddratio_lt, 'PValue': pval_lt, 
                         "FDR": fdr_lt, 'Test': test_lt, 'Env': env_lt})

def main():
    calculate_enrichment().to_csv('vmr_enrichment_analysis.txt', 
                                  sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
