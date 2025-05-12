"""
This script merges phenotypes into a single file.
"""
import pandas as pd

def main():
    pd.concat((pd.read_csv('../../_m/%s_phenotypes.csv' % x, index_col=0)
               for x in ['caudate', 'dlpfc', 'hippo', 'dg']))\
      .to_csv("merged_phenotypes.csv", index=True)


if __name__ == "__main__":
    main()
