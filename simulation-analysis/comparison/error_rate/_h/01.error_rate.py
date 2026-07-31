import numpy as np
import polars as pl
import session_info
from pyhere import here
from functools import lru_cache
from sklearn.metrics import confusion_matrix
from statsmodels.stats.multitest import multipletests

# Constants
THRESHOLDS = {
    "heritability": 0.1,
    "r_squared_cv": 0.75,
    "fdr": 0.05
}

@lru_cache(maxsize=None)
def load_file_cached(path):
    return pl.read_csv(path, separator="\t")


def get_ground_truth_path(n_samples):
    return here("inputs/simulated-data/_m", f"sim_{n_samples}_indiv",
                "snp_phenotype_mapping.tsv")

def get_predicted_path(method, n_samples):
    if method == "elastic-net":
        return here(f"simulation-analysis/{method}/sim_{n_samples}_indiv/_m",
                    f"simulation_{n_samples}_summary_elastic-net.tsv")
    elif method == "gcta":
        return here(f"simulation-analysis/{method}/_m",
                    "summary/greml_summary.tsv")
    else:
        raise ValueError(f"Unsupported method: {method}")


def load_ground_truth(n_samples, heritable):
    df = load_file_cached(get_ground_truth_path(n_samples))
    df = df.with_columns((pl.col("phenotype_id").alias("pheno_id")))\
           .select(pl.col(["pheno_id", "simulated_heritability"]))
    condition = pl.col("simulated_heritability") >= THRESHOLDS["heritability"] if heritable \
                else pl.col("simulated_heritability") < THRESHOLDS["heritability"]
    return df.with_columns(condition.alias("truth"))


def load_predicted(n_samples, method, heritable):
    fn = get_predicted_path(method, n_samples)
    df = load_file_cached(fn)

    if method == "elastic-net":
        df = df.select(["pheno_id", "h2_unscaled", "r_squared_cv"])
        condition = (
            (pl.col("h2_unscaled") >= THRESHOLDS["heritability"]) &
            (pl.col("r_squared_cv") >= THRESHOLDS["r_squared_cv"])
        ) if heritable else (
            (pl.col("h2_unscaled") < THRESHOLDS["heritability"]) &
            (pl.col("r_squared_cv") < THRESHOLDS["r_squared_cv"])
        )
        return df.with_columns(condition.alias("predicted"))
    elif method == "gcta":
        df = df.select(["N","pheno","Sum of V(G)_Vp_Variance","Pval_Variance"])\
               .rename({"pheno": "pheno_id",
                        "Sum of V(G)_Vp_Variance": "h2_unscaled"})
        pvals = df["Pval_Variance"].to_numpy()
        df = df.with_columns(
            pl.Series("pval_fdr_adjusted", multipletests(pvals, method="fdr_bh")[1])
        )
        condition = (
            (pl.col("h2_unscaled") >= THRESHOLDS["heritability"]) &
            (pl.col("pval_fdr_adjusted") < THRESHOLDS["fdr"]) &
            (pl.col("N") == n_samples)
        ) if heritable else (
            (pl.col("h2_unscaled") < THRESHOLDS["heritability"]) &
            (pl.col("pval_fdr_adjusted") < THRESHOLDS["fdr"]) &
            (pl.col("N") == n_samples)
        )
        return df.with_columns(condition.alias("predicted"))


def merge_dataframes(n_samples, method, heritable):
    true_df = load_ground_truth(n_samples, heritable)
    pred_df = load_predicted(n_samples, method, heritable)
    merged_df = true_df.join(pred_df, on="pheno_id", how="inner")
    return merged_df.with_columns(pl.col("simulated_heritability").alias("h2_simulated"))\
                    .select(["pheno_id","h2_simulated","h2_unscaled","truth","predicted"])


def get_confusion_matrix_elements(n_samples, method, heritable):
    df = merge_dataframes(n_samples, method, heritable)
    y_true = df["truth"].to_numpy()
    y_pred = df["predicted"].to_numpy()
    cm = confusion_matrix(y_true, y_pred, labels=[False, True])
    return cm.ravel() if cm.size == 4 else (0, 0, 0, 0)


def calculate_error_and_power(n_samples, method="elastic-net", heritable=True):
    tn, fp, fn, tp = get_confusion_matrix_elements(n_samples, method, heritable)
    type1_error = fp / (fp + tn) if (fp + tn) > 0 else np.nan
    type2_error = fn / (fn + tp)  if (fn + tp) > 0 else np.nan
    power = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    return power, type1_error, type2_error


def main():
    results = []
    for method in ["elastic-net", "gcta"]:
        for n in [100, 150, 200, 250, 500, 1_000, 5_000]:
            power, type1, type2 = calculate_error_and_power(n, method)
            n_tested = (1_000 if method != "elastic-net" else
                        merge_dataframes(n, method, heritable=True).shape[0])
            results.append({
                'sample_size': n,
                'power': power,
                'type1_error': type1,
                'type2_error': type2,
                'n_tested': n_tested,
                'method': method
            })
    power_df = pl.DataFrame(results)
    power_df.write_csv("power-analysis.tsv", separator="\t")
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
