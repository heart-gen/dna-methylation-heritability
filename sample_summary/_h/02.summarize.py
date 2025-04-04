"""
This script aggregates data by diagnosis and prints
summary tables to the log file.
"""
import session_info
import polars as pl
from pyhere import here

def get_pheno():
    return pl.read_csv(here("sample_summary","_m/phenotype_data.tsv"),
                       separator="\t", has_header=True)


def main():
    ## Main
    with open("sample_breakdown.log", "w") as f:
        dx = get_pheno().select(["Dx", "Race", "Sex", "Age", "RIN"])
        print(dx.group_by(["Dx"]).agg(pl.len()).to_pandas().to_string(),
              file=f)
        print(dx.group_by(["Dx", "Sex"]).agg(pl.len()).to_pandas().to_string(),
              file=f)
        print(dx.group_by(["Dx", "Race"]).agg(pl.len()).to_pandas().to_string(),
              file=f)
        print("Mean:", file=f)
        print(dx.group_by(["Dx"]).agg([pl.col("Age").mean().alias("Age_mean"),
                                       pl.col("RIN").mean().alias("RIN_mean")])\
              .to_pandas().to_string(), file=f)
        print("Standard deviation", file=f)
        print(dx.group_by(["Dx"]).agg([pl.col("Age").std().alias("Age_mean"),
                                       pl.col("RIN").std().alias("RIN_mean")])\
              .to_pandas().to_string(), file=f)
        ## Session information
        session_info.show()


if __name__ == "__main__":
    main()
