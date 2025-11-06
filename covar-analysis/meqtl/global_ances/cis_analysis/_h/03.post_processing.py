"""
This script runs tensorQTL in python.
"""

import pandas as pd
import session_info
from tensorqtl.post import get_significant_pairs

print(f"Pandas {pd.__version__}")

def get_permutation_results(prefix):
    return pd.read_csv(f"{prefix}.permutation.txt.gz", sep="\t", index_col=0)


def get_eqtl(prefix, perm_df):
    return get_significant_pairs(perm_df, prefix)


def main():
    prefix = "TOPMed_LIBD"

    # Load the permutation results
    perm_df = get_permutation_results(prefix)

    # Load and save eQTL results
    get_eqtl(prefix, perm_df)\
        .to_csv(f"{prefix}.signif_variants.txt.gz" , sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
