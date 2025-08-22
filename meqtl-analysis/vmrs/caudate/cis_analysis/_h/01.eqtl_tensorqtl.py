"""This script runs tensorQTL in Python."""
import pandas as pd
import session_info
from functools import lru_cache
from tensorqtl import cis, pgen, read_phenotype_bed

@lru_cache()
def get_genotypes():
    plink_prefix_path = "../../_m/protected_data/TOPMed_LIBD"
    pgr = pgen.PgenReader(plink_prefix_path)
    variant_df = pgr.variant_df
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pgr.load_genotypes(), variant_df


@lru_cache()
def get_covars(feature = "genes"):
    covar_file = f"../../_m/{feature}.combined_covariates.txt"
    return pd.read_csv(covar_file, sep='\t', index_col=0).T


@lru_cache()
def get_phenotype(feature = "genes"):
    expr_bed = f"../../_m/{feature}.expression.bed.gz"
    return read_phenotype_bed(expr_bed)


def main():
    # Load data
    feature = "genes"; prefix = "TOPMed_LIBD"
    phenotype_df, phenotype_pos_df = get_phenotype(feature)
    genotype_df, variant_df = get_genotypes()

    # Nominal
    cis.map_nominal(genotype_df, variant_df, phenotype_df,
                    phenotype_pos_df, prefix, covariates_df=get_covars(feature),
                    maf_threshold=0.01, window=500000, output_dir=".")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
