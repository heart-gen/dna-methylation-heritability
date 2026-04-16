"""
Error rate analysis across LD decay levels at fixed N=200.

Compares power / type-1 / type-2 error for:
  - boosting_hybrid
  - joint_ridge
across LD decay levels 0.5, 0.6, 0.7, 0.8.

The LD=0.8 baseline uses the original sim_200_indiv pipeline output
(no METHOD label in filename), representing the current approach.
"""

import numpy as np
import polars as pl
import session_info
from pyhere import here
from functools import lru_cache
from sklearn.metrics import confusion_matrix

# Constants
LD_LEVELS = ["0.8", "0.7", "0.6", "0.5"]
METHODS   = ["boosting_hybrid", "joint_ridge"]

THRESHOLDS = {
    "heritability": 0.025,
    "r_squared_cv": 0.3,
}


@lru_cache(maxsize=None)
def load_file_cached(path):
    return pl.read_csv(path, separator="\t")


def get_ground_truth_path(ld_decay):
    if ld_decay == "0.8":
        return here(
            "inputs/simulated-data/_m/sim_200_indiv",
            "snp_phenotype_mapping.tsv"
        )
    return here(
        "inputs/simulated-data/_m",
        f"ld_{ld_decay}_sim_200_indiv",
        "snp_phenotype_mapping.tsv"
    )


def get_predicted_path(ld_decay, model_method):
    """LD=0.8 baseline uses sim_200_indiv (original pipeline, no METHOD label)."""
    if ld_decay == "0.8":
        return here(
            "simulation-analysis/elastic-net/sim_200_indiv/_m",
            "simulation_200_summary_elastic-net.tsv"
        )
    return here(
        "simulation-analysis/elastic-net/ld_decay/_m",
        f"simulation_200_{ld_decay}_{model_method}_summary_elastic-net.tsv"
    )


def load_ground_truth(ld_decay, heritable):
    df = load_file_cached(get_ground_truth_path(ld_decay))
    # Use simulated_heritability (realized h2) for classification threshold
    df = df.with_columns(pl.col("phenotype_id").alias("pheno_id")) \
           .select(["pheno_id", "simulated_heritability"])
    condition = (
        pl.col("simulated_heritability") >= THRESHOLDS["heritability"]
        if heritable
        else pl.col("simulated_heritability") < THRESHOLDS["heritability"]
    )
    return df.with_columns(condition.alias("truth"))


def load_predicted(ld_decay, model_method, heritable):
    fn = get_predicted_path(ld_decay, model_method)
    df = load_file_cached(fn)
    df = df.select(["pheno_id", "h2_unscaled", "r_squared_cv"])
    # Null comparisons in Polars return null, not False.
    # fill_null(False) treats any row with a missing value as "not predicted"
    # (equivalent to "low prediction" in the R classification framework).
    condition = (
        (pl.col("h2_unscaled") >= THRESHOLDS["heritability"]) &
        (pl.col("r_squared_cv") >= THRESHOLDS["r_squared_cv"])
        if heritable
        else (
            (pl.col("h2_unscaled") < THRESHOLDS["heritability"]) |
            (pl.col("r_squared_cv") < THRESHOLDS["r_squared_cv"])
        )
    ).fill_null(False)
    return df.with_columns(condition.alias("predicted"))


def merge_dataframes(ld_decay, model_method, heritable):
    true_df = load_ground_truth(ld_decay, heritable)
    pred_df = load_predicted(ld_decay, model_method, heritable)
    merged  = true_df.join(pred_df, on="pheno_id", how="inner")
    return merged.with_columns(
        pl.col("simulated_heritability").alias("h2_simulated")
    ).select(["pheno_id", "h2_simulated", "h2_unscaled", "truth", "predicted"])


def get_confusion_matrix_elements(ld_decay, model_method, heritable):
    df     = merge_dataframes(ld_decay, model_method, heritable)
    y_true = df["truth"].to_numpy()
    y_pred = df["predicted"].to_numpy()
    cm     = confusion_matrix(y_true, y_pred, labels=[False, True])
    return cm.ravel() if cm.size == 4 else (0, 0, 0, 0)


def calculate_error_and_power(ld_decay, model_method, heritable=True):
    tn, fp, fn, tp = get_confusion_matrix_elements(ld_decay, model_method, heritable)
    type1_error = fp / (fp + tn) if (fp + tn) > 0 else np.nan
    type2_error = fn / (fn + tp) if (fn + tp) > 0 else np.nan
    power       = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    return power, type1_error, type2_error


def main():
    results = []
    for ld in LD_LEVELS:
        for method in METHODS:
            power, type1, type2 = calculate_error_and_power(ld, method, heritable=True)
            n_tested = merge_dataframes(ld, method, heritable=True).shape[0]
            results.append({
                "ld_decay":    ld,
                "method":      method,
                "power":       power,
                "type1_error": type1,
                "type2_error": type2,
                "n_tested":    n_tested,
            })

    power_df = pl.DataFrame(results)
    print(power_df)
    power_df.write_csv("ld_decay_power-analysis.tsv", separator="\t")

    session_info.show()


if __name__ == "__main__":
    main()
