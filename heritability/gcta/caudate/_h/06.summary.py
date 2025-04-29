import session_info
import pandas as pd
from glob import glob

import os
import argparse
from pyhere import here

def read_hsq_file(fn):
    return pd.read_csv(fn, sep="\t")


def merge_data():
    basepath = "/projects/p32505/users/alexis/projects/" +\
        "dna-methylation-heritability/heritability/gcta/caudate/_m"
    dt = pd.DataFrame()
    for fn in iglob(f"{basepath}/summary/chr_*/*csv"):
        dt = pd.concat([dt, read_hsq_file(fn)], axis=0)
    return dt


def process_hsq_file(filename):
    greml_df = pd.read_csv(filename, sep="\t")\
                 .melt(id_vars="Source", value_vars=["Variance", "SE"],
                       var_name="Stat", value_name="Value")\
                 .pivot_table(index=None, columns=["Source", "Stat"],
                              values="Value").dropna(axis=1)
    greml_df.columns = [f"{src.replace('/', '_')}_{stat}" for src, stat in df.columns]
    return df


def main():
    ## Main
    basepath = Path("../_m/") / "h2"
    dt = pd.DataFrame()
    for fn in basepath.glob("*/*hsq"):
        greml_df = process_hsq_file(fn)
        # extract chromsome from filename
        # Extract start and end from filename
        # Add to data frame
        dt = pd.concat([dt, greml_df], axis=0)

    ## Save
    ## Session information
    session_info.show()


if __name__ == '__main__':
    main()
