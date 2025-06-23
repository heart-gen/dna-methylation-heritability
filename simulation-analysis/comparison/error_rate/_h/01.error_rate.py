import numpy as np
import polars as pl
import session_info
from pyhere import here
from sklearn.metrics import confusion_matrix

def load_ground_truth(n_samples, heritable):
    fn = here("inputs/simulated-data/_m",
              f"sim_{n_samples}_indiv",
              "snp_phenotype_mapping.tsv")
    df = pl.read_csv(fn, separator="\t")\
           .with_columns((pl.col("phenotype_id").alias("pheno_id")))\
           .select(pl.col(["pheno_id", "simulated_heritability"]))
    if heritable:
        return df.with_columns((pl.col("simulated_heritability") >= 0.1).alias("truth"))
    else:
        return df.with_columns((pl.col("simulated_heritability") < 0.1).alias("truth"))


def load_predicted(n_samples, method, heritable):
    if method == "elastic-net":
        fn = f"../../../{method}/"+\
            f"sim_{n_samples}_indiv/_m/"+\
            f"simulation_{n_samples}_summary_elastic-net.tsv"
        df = pl.read_csv(fn, separator="\t")\
               .select(["pheno_id", "h2_unscaled", "r_squared_cv"])
        if heritable:
            return df.with_columns((
                (pl.col("h2_unscaled") >= 0.1) &
                (pl.col("r_squared_cv") >= 0.75)
            ).alias("predicted"))
        else:
            return df.with_columns((
                (pl.col("h2_unscaled") < 0.1) &
                (pl.col("r_squared_cv") < 0.75)
            ).alias("predicted"))
    else:
        fn = None
        return None


def merge_dataframes(n_samples, method, heritable):
    true_df = load_ground_truth(n_samples, heritable)
    pred_df = load_predicted(n_samples, method, heritable)
    return true_df.join(pred_df, on="pheno_id", how="inner")\
                  .with_columns(pl.col("simulated_heritability").alias("h2_simulated"))\
                  .select(["pheno_id", "h2_simulated",
                           "h2_unscaled", "truth",
                           "predicted"])


def get_confusion_matrix_elements(n_samples, method, heritable):
    df = merge_dataframes(n_samples, method, heritable)
    y_true = df["truth"].to_numpy()
    y_pred = df["predicted"].to_numpy()
    return confusion_matrix(y_true, y_pred).ravel()


def calculate_error_rate(n_samples, method="elastic-net", heritable = True):
    tn, fp, fn, tp = get_confusion_matrix_elements(n_samples, method, heritable)
    type1_error_rate = fp / (fp + tn) if (fp + tn) > 0 else np.nan
    type2_error_rate = fn / (fn + tp)  if (fn + tp) > 0 else np.nan
    return type1_error_rate, type2_error_rate


def empirical_power(n_samples, method="elastic-net", heritable = True):
    tn, fp, fn, tp = get_confusion_matrix_elements(n_samples, method, heritable)
    power = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    type1_error = fp / (fp + tn) if (fp + tn) > 0 else np.nan
    return power, type1_error


def main():
    results = []; method = "elastic-net"
    for n in [100, 150, 200, 250, 500, 1000, 5000, 10_000]:
        power, type1_error = empirical_power(n, method)
        results.append({
            'sample_size': n,
            'power': power,
            'type1_error': type1_error,
            'n_tested': 1000,
            'method': method
        })
    power_df = pl.DataFrame(results)


if __name__ == "__main__":
    main()
