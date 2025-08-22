"""
This script is used to prepare GCT files for the GTEx eQTL pipeline.
Specifically for converting BrainSEQ genotypes using TOPMed imputation
as well as the normalized files provided by BrainSEQ.

- Outputs:
  * GCT files for counts and TPM for selected genes and samples
  * A lookup table of sample_ids and brain_ids
  * A list of chromosomes to use
"""

import session_info
import pandas as pd
from pyhere import here
from functools import lru_cache
from sklearn.feature_selection import VarianceThreshold

def to_gct(filename, df):
    description_df = pd.DataFrame({'Description': df.index.values},
                                  index=df.index)
    dfo = pd.concat([description_df, df], axis=1)
    dfo.index.name = 'Names'
    with open(filename, "wt") as out:
        print("#1.2", file=out)
        print(df.shape[0], df.shape[1], sep="\t", file=out)
        dfo.to_csv(out, sep="\t")


def remove_near_constant(df, threshold=1e-3):
    selector = VarianceThreshold(threshold=threshold)
    filtered = selector.fit_transform(df.T.values)
    mask = selector.get_support(indices=True)
    return df.iloc[mask, :]


@lru_cache()
def get_pheno():
    return pd.read_csv("phenotypes.csv", index_col=0)


@lru_cache()
def get_psam():
    psam_file = here("inputs/genotypes/_m/TOPMed_LIBD.psam")
    return pd.read_csv(psam_file, sep="\t", usecols=[0, 1])


@lru_cache()
def load_data():
    pheno_df = get_pheno()
    pheno_df["ids"] = pheno_df.RNum
    pheno_df.set_index("ids", inplace=True)
    norm_df = pd.read_csv("normalized_expression.txt.gz",
                          sep="\t", index_col=0)
    samples = list(set(norm_df.columns)\
                   .intersection(set(pheno_df["RNum"])))
    return pheno_df.loc[samples,:], norm_df.loc[:,samples]


def select_idv(pheno_df, norm_df):
    samples = list(set(pheno_df.loc[norm_df.columns,:].BrNum)\
                   .intersection(set(get_psam()["#FID"])))
    new_psam = get_psam()[(get_psam()["#FID"].isin(samples))]\
        .drop_duplicates(subset="#FID")
    new_psam.to_csv("keepPsam.txt", sep='\t', index=False,
                    header=True)
    return pheno_df.loc[:, ["RNum", "BrNum"]]\
                   .reset_index().set_index("BrNum")\
                   .loc[new_psam["#FID"]].reset_index().set_index("ids")


def main():
    ## Load data
    pheno_df, norm_df = load_data()
    ## Match individuals
    new_pheno = select_idv(pheno_df, norm_df)
    ## Remove near constant features
    norm_df = remove_near_constant(norm_df.loc[:, new_pheno.index])
    ## Save with GCT format
    to_gct("norm.gct", norm_df)
    new_pheno.loc[:, ["RNum", "BrNum"]]\
             .to_csv("sample_id_to_brnum.tsv", sep="\t", index=False)
    ## Save chrom list
    pd.DataFrame({'chr':['chr'+xx for xx in [str(x) for x in range(1,23)]+['X']]})\
      .to_csv('vcf_chr_list.txt', header=False, index=None)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
