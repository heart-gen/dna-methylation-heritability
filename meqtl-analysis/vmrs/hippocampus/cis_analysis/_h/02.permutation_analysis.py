"""This script runs tensorQTL in Python."""
import pandas as pd
import session_info
from sys import stdout
from functools import lru_cache
from tensorqtl import cis, pgen, read_phenotype_bed

class SimpleLogger(object):
    def __init__(self, logfile=None, verbose=True):
        self.console = stdout
        self.verbose = verbose
        if logfile is not None:
            self.log = open(logfile, 'w')
        else:
            self.log = None

    def write(self, message):
        if self.verbose:
            self.console.write(message+'\n')
        if self.log is not None:
            self.log.write(message+'\n')
            self.log.flush()


@lru_cache()
def get_genotypes():
    plink_prefix_path = "../../_m/protected_data/TOPMed_LIBD"
    pgr = pgen.PgenReader(plink_prefix_path)
    variant_df = pgr.variant_df
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pgr.load_genotypes(), variant_df


@lru_cache()
def get_covars(feature = "vmrs"):
    covar_file = f"../../_m/{feature}.combined_covariates.txt"
    return pd.read_csv(covar_file, sep='\t', index_col=0).T


@lru_cache()
def get_phenotype(feature = "vmrs"):
    meth_bed = f"../../_m/{feature}.methylation.bed.gz"
    return read_phenotype_bed(meth_bed)


def calculate_qvalues(res_df, fdr=0.05, qvalue_lambda=None, logger=None):
    """Annotate permutation results with q-values, p-value threshold"""
    from py_qvalue import qvalue
    from scipy.stats import pearsonr, beta

    logger = logger or SimpleLogger()
    logger.write('Computing q-values')
    logger.write(f'  * Number of phenotypes tested: {res_df.shape[0]}')

    # Determine which p-value column to use
    if res_df['pval_beta'].notnull().any():
        pval_col = 'pval_beta'
        corr = pearsonr(res_df['pval_perm'], res_df['pval_beta'])[0]
        logger.write(f'  * Correlation between Beta-approximated and empirical p-values: {corr:.4f}')
    else:
        pval_col = 'pval_perm'
        logger.write(f'  * WARNING: no beta-approximated p-values found, using permutation p-values instead.')

    # Calculate q-values
    if qvalue_lambda is not None:
        logger.write(f'  * Calculating q-values with lambda = {qvalue_lambda:.3f}')
        qval_res = qvalue(res_df[pval_col], lambda_=qvalue_lambda)
    else:
        qval_res = qvalue(res_df[pval_col])

    res_df['qval'] = qval_res["qvalues"]
    pi0 = qval_res["pi0"]

    logger.write(f'  * Proportion of significant phenotypes (1-pi0): {1-pi0:.2f}')
    logger.write(f"  * QTL phenotypes @ FDR {fdr:.2f}: {(res_df['qval'] <= fdr).sum()}")

    # Calculate nominal p-value threshold if using beta p-values
    if pval_col == 'pval_beta':
        sig = res_df.loc[res_df['qval'] <= fdr, 'pval_beta']
        nonsig = res_df.loc[res_df['qval'] > fdr, 'pval_beta']
        if not sig.empty:
            lb = sig.max()
            ub = nonsig.min() if not nonsig.empty else lb
            pthreshold = (lb + ub) / 2 if ub != lb else lb
            logger.write(f'  * min p-value threshold @ FDR {fdr}: {pthreshold:.6g}')
            res_df['pval_nominal_threshold'] = beta.ppf(pthreshold, res_df['beta_shape1'], res_df['beta_shape2'])


def main():
    # Load data
    feature = "vmrs"; prefix = "TOPMed_LIBD"
    phenotype_df, phenotype_pos_df = get_phenotype(feature)
    genotype_df, variant_df = get_genotypes()

    # Permutation
    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df,
                         phenotype_pos_df, covariates_df=get_covars(feature),
                         maf_threshold=0.01, window=500000, seed=13131313)

    calculate_qvalues(cis_df, fdr=0.05)
    cis_df.to_csv(f"{prefix}.permutation.txt.gz", sep='\t')

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
