"""
This script aggregates data by diagnosis and prints
summary tables to the log file.
"""
import session_info
import polars as pl
from pyhere import here

def get_pheno():
    return pl.read_csv(here("sample_summary","_m/phenotype_data.tsv"),
                       separator="\t", has_header=True, infer_schema_length=10000)

def main():
    ## Main
    with open("sample_breakdown.log", "w") as f:
        dx = get_pheno().select(["race", "region", "primarydx", "sex", "agedeath"])
        print(dx.group_by(["region"]).agg(pl.len()).to_pandas().to_string(),
              file=f)
        print(dx.group_by(["race", "region"]).agg(pl.len()).to_pandas().to_string(),
              file=f)
        print(dx.group_by(["race", "region", "primarydx"]).agg(pl.len()).to_pandas().to_string(),
              file=f)
        print(dx.group_by(["race", "region", "primarydx", "sex"]).agg(pl.len()).to_pandas().to_string(),
              file=f)
        print("Mean:", file=f)
        print(dx.group_by(["race", "region"]).agg([pl.col("agedeath").mean().alias("Age_mean")])\
              .to_pandas().to_string(), file=f)
        print("Standard deviation", file=f)
        print(dx.group_by(["race", "region"]).agg([pl.col("agedeath").std().alias("Age_sd")])\
              .to_pandas().to_string(), file=f)
        ## Session information
        session_info.show()


if __name__ == "__main__":
    main()
