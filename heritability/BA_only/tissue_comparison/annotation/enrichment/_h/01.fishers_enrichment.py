"""
This script performs Fisher's pairwise enrichment analysis
between genomic annotation type and heritable sites estimated by
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
    enet_fn = here(f"heritability/elastic_net_model/BA_only/{tissue.lower()}/_m/{tissue.lower()}_summary_elastic-net.tsv")
    df = pd.read_csv(enet_fn, sep='\t')
    df['chrom'] = 'chr' + df['chrom'].astype(str)

    # assign h2 categories
    df['h2_category'] = df.apply(
    lambda row: (
        'Heritable' if row['h2_unscaled'] >= 0.1 and row['r_squared_cv'] > 0.3 else
        'Non-heritable' if row['h2_unscaled'] < 0.1 and row['r_squared_cv'] > 0.3 else
        'Low prediction'
    ),
    axis=1
    )
    return df

@lru_cache()
def get_annotation(tissue):
    fn = here(f"heritability/elastic_net_model/BA_only/tissue_comparison/annotation/_m/{tissue.lower()}_vmr_annotations_hg38_wide.tsv")
    return pd.read_csv(fn, sep='\t')

@lru_cache()
def merge_dataframe(tissue):
    return get_enet(tissue).merge(get_annotation(tissue), 
                                  left_on=['chrom', 'start', 'end', 'h2_category'],
                                  right_on=['seqnames', 'start', 'end', 'h2_category'],
                                  how='inner')

def cal_fishers_annot(annot, h2_cat, tissue):
    df = merge_dataframe(tissue)
    
    table = [[np.sum((df[annot] == 1) & (df['h2_category'] == h2_cat)),
              np.sum((df[annot] == 1) & (df['h2_category'] != h2_cat))],
             [np.sum((df[annot] == 0) & (df['h2_category'] == h2_cat)),
              np.sum((df[annot] == 0) & (df['h2_category'] != h2_cat))]]
    print(table)
    return fisher_exact(table)

def calculate_enrichment():
    region_lt = []; h2_lt = []; annot_lt = []; fdr_lt = []; pval_lt = []; oddratio_lt = []

    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        for h2_cat in ["Heritable", "Non-heritable", "Low prediction"]:
            pvals = []
            df = merge_dataframe(tissue)
            annot_cols = [col for col in df.columns if 'hg38' in col]
            for annot in annot_cols:
                    odd_ratio, pval = cal_fishers_annot(annot, h2_cat, tissue)
                    pvals.append(pval); h2_lt.append(h2_cat)
                    oddratio_lt.append(odd_ratio); annot_lt.append(annot)
                    region_lt.append(tissue)
            _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
            pval_lt = np.concatenate((pval_lt, pvals))
            fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Tissue': region_lt, 'h2_Category': h2_lt, 
                         'OR': oddratio_lt, 'PValue': pval_lt, 
                         "FDR": fdr_lt, 'Annotation': annot_lt})

def main():
    calculate_enrichment().to_csv('annotation_vmr_enrichment_analysis.txt', sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()