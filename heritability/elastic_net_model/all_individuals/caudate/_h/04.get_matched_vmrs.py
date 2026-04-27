"""
This script performs Fisher's pairwise enrichment analysis
between environmental features and heritable sites estimated by
elastic-net.
"""
import numpy as np
import pandas as pd
from pyhere import here
from functools import lru_cache
import sys
import session_info

def filter_sites(enet_df, r2_threshold):
    enet_df = enet_df.dropna()
    enet_df['chrom'] = 'chr' + enet_df['chrom'].astype(str)

    # assign h2 categories
    enet_df['h2_category'] = enet_df.apply(
    lambda row: (
        'Heritable' if row['h2_unscaled'] >= 0.1 and row['r_squared_cv'] > r2_threshold 
        else
        'Non-heritable' if row['h2_unscaled'] < 0.1 and row['r_squared_cv'] > r2_threshold 
        else
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
def get_matching_enet(tissue, r2_threshold):
    # get vmrs
    enet_AA_fn = here(f"heritability/elastic_net_model/all_individuals/{tissue.lower()}/_m/{tissue.lower()}_summary_elastic-net_AA.tsv")
    enet_EA_fn = here(f"heritability/elastic_net_model/all_individuals/{tissue.lower()}/_m/{tissue.lower()}_summary_elastic-net_EA.tsv")

    df_AA = pd.read_csv(enet_AA_fn, sep='\t')
    df_EA = pd.read_csv(enet_EA_fn, sep='\t')

    print(f"[ENET] AA shape (raw): {df_AA.shape}")
    print(f"[ENET] EA shape (raw): {df_EA.shape}")

    enet_AA = filter_sites(df_AA, r2_threshold)
    enet_EA = filter_sites(df_EA, r2_threshold)
    
    enet_combined = enet_AA.merge(enet_EA, on='feature_id', 
                                  suffixes=('_AA', '_EA'))
    
    print(f"[ENET] Combined shape (after merge): {enet_combined.shape}")
    
    enet_combined = enet_combined[
        enet_combined['h2_category_AA'] == enet_combined['h2_category_EA']
    ]

    print(f"[ENET] After h2_category match filter: {enet_combined.shape}")

    enet_combined = enet_combined.rename(columns={
        'h2_category_AA': 'h2_category',
        "region_AA": "region",
        'chrom_AA': 'chrom',
        'start_AA': 'start',
        'end_AA': 'end'})
    
    enet_combined = enet_combined.drop(columns=[
        'h2_category_EA', 'chrom_EA', 'start_EA', 
        'end_EA', 'race_EA', 'race_AA', "region_EA"])
    
    print(f"[ENET] Final columns: {enet_combined.columns.tolist()}")
    print(f"[ENET] Final shape: {enet_combined.shape}")
 
    return enet_combined

def main():
    # Get population from command line
    tissue = sys.argv[1]

    for r2_threshold in 0.3, 0.75:
        get_matching_enet(tissue, r2_threshold).to_csv(f"{tissue.lower()}_summary_elastic-net_matched_r2_{r2_threshold}.tsv", sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
