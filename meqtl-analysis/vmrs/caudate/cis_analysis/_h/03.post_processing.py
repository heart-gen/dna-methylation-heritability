"""
This script runs tensorQTL in python.
"""

import pandas as pd
import session_info
import glob
from functools import lru_cache
from tensorqtl.post import get_significant_pairs
from statsmodels.stats.multitest import fdrcorrection

print(f"Pandas {pd.__version__}")

@lru_cache()
def get_nominal(prefix):
   
    # Load parquet files
    parquet_files = glob.glob(f'{prefix}.cis_qtl_pairs.*.parquet')
    nominal_df = pd.concat([pd.read_parquet(f) for f in parquet_files])

    # FDR correction
    _, fdr = fdrcorrection(nominal_df["pval_nominal"])
    nominal_df["fdr"] = fdr

    return nominal_df[nominal_df["fdr"] < 0.05]

def get_permutation_results(prefix):
    return pd.read_csv(f"{prefix}.permutation.txt.gz", sep="\t", index_col=0)


def get_eqtl(prefix, perm_df):
    return get_significant_pairs(perm_df, prefix)


def main():
    prefix = "TOPMed_LIBD"

    # Load and save nominal results
    get_nominal(prefix)\
        .to_csv(f"{prefix}.nominal_signif_variants.txt.gz" , 
                sep='\t', index=False)

    # Load the permutation results
    perm_df = get_permutation_results(prefix)

    # Load and save eQTL results
    get_eqtl(prefix, perm_df)\
        .to_csv(f"{prefix}.permutation_signif_variants.txt.gz" , sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
