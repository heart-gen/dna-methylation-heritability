"""
Fisher's exact enrichment of environmental associations within
AA non-heritable VMR subgroups A, B, and C.

Uses the all_individuals environmental analysis data (logit + DMR
results from environmental-analysis/all_individuals/{tissue}/correlation/).

Background for each Fisher's test: all AA non-heritable VMRs in the tissue
(Groups A + B + C combined).

Foreground: VMRs in the target subgroup that are significantly associated
with the environmental variable (p < 0.05, unadjusted, consistent with
existing enrichment pipeline).

2×2 contingency table:
                     In subgroup   Not in subgroup
  sig association        a               b
  not sig                c               d

OR > 1  →  environmental signal enriched in this subgroup
OR < 1  →  depleted

FDR corrected across all comparisons within each tissue.

Output:
  subgroup_enrichment_analysis.tsv
"""

import numpy as np
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection
import session_info

TISSUE_DIRS = {"Caudate": "caudate", "DLPFC": "dlpfc", "Hippocampus": "hippocampus"}


@lru_cache()
def get_subgroups(tissue):
    fn = here(
        f"environmental-analysis/all_individuals/tissue_compare/"
        f"non_heritable_subgroup/_m/non_heritable_subgroups_{tissue}.tsv"
    )
    df = pd.read_csv(fn, sep="\t")
    df["chr"] = df["chrom"].astype(str)
    if not df["chr"].str.startswith("chr").all():
        df["chr"] = "chr" + df["chr"]
    return df[["chr", "start", "end", "feature_id", "subgroup"]]


@lru_cache()
def get_vmrs(tissue, env, cat=None):
    tdir = TISSUE_DIRS[tissue]
    fn = here(
        f"environmental-analysis/all_individuals/{tdir}/"
        f"correlation/_m/{env.lower()}_logit.csv.gz"
    )
    vmrs = pd.read_csv(fn, sep=",")
    if cat is not None:
        vmrs = vmrs[vmrs["var"] == f"{env}{cat}"]
    vmrs["test"] = "logit"
    vmrs["sig"] = vmrs["p"] < 0.05
    return vmrs


@lru_cache()
def get_dmrs(tissue, env, cat=None):
    tdir = TISSUE_DIRS[tissue]
    fn = here(
        f"environmental-analysis/all_individuals/{tdir}/"
        f"correlation/_m/{env.lower()}_dmr.csv.gz"
    )
    dmrs = pd.read_csv(fn, sep=",")
    if cat is not None:
        dmrs = dmrs[dmrs["level"] == cat]
    dmrs["test"] = "dmr"
    dmrs["sig"] = dmrs["p.value"] < 0.05
    return dmrs


@lru_cache()
def build_merged(tissue, env, cat=None):
    """
    One row per VMR with columns: chr, start, end, logit, dmr, both, subgroup.
    Background = all AA non-heritable VMRs (Groups A + B + C).
    """
    vmrs = get_vmrs(tissue, env, cat)
    dmrs = get_dmrs(tissue, env, cat)

    combined = pd.concat([vmrs, dmrs])
    combined["chr"] = "chr" + combined["chr"].astype(str)

    pivot = combined.pivot_table(
        values="sig", columns="test", fill_value=0,
        index=["chr", "start", "end"]
    ).reset_index()

    for col in ["logit", "dmr"]:
        if col not in pivot.columns:
            pivot[col] = 0

    pivot["both"] = ((pivot["logit"] == 1) & (pivot["dmr"] == 1)).astype(int)

    sg = get_subgroups(tissue)
    merged = sg.merge(pivot, on=["chr", "start", "end"], how="left")

    # VMRs with no association result default to not-sig
    for col in ["logit", "dmr", "both"]:
        merged[col] = merged[col].fillna(0).astype(int)

    return merged


def fisher_subgroup(tissue, test, subgroup, env, cat=None):
    df = build_merged(tissue, env, cat)
    in_grp = df["subgroup"] == subgroup
    sig = df[test.lower()] == 1

    # 2×2: rows = sig/not, cols = in_grp/not_in_grp
    table = [
        [int((~sig & ~in_grp).sum()), int((~sig & in_grp).sum())],
        [int(( sig & ~in_grp).sum()), int(( sig & in_grp).sum())],
    ]
    return fisher_exact(table)


def calculate_enrichment():
    records = []

    env_vars = [
        "smoking", "codeine", "morphine", "cocaine", "ethanol",
        "antipsychotics", "nicotine", "amphetamines",
        "hx_sexual_abuse", "hx_physical_abuse",
    ]
    categorical_vars = {
        "education": ["less_than_hs", "more_than_hs"],
        "marital_status": ["single", "previously_married"],
    }
    subgroups = ["A", "B", "C"]
    tests = ["Logit", "DMR", "Both"]

    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        pvals_all = []
        meta_all = []

        for subgroup in subgroups:
            for test in tests:
                for env in env_vars:
                    try:
                        or_, pval = fisher_subgroup(tissue, test, subgroup, env)
                    except Exception as e:
                        print(f"  SKIP {tissue}/{subgroup}/{test}/{env}: {e}")
                        or_, pval = np.nan, np.nan
                    pvals_all.append(pval)
                    meta_all.append(
                        dict(Tissue=tissue, Subgroup=subgroup, Test=test,
                             Env=env, OR=or_, PValue=pval)
                    )

                for env, cats in categorical_vars.items():
                    for cat in cats:
                        try:
                            or_, pval = fisher_subgroup(tissue, test, subgroup, env, cat)
                        except Exception as e:
                            print(f"  SKIP {tissue}/{subgroup}/{test}/{env}/{cat}: {e}")
                            or_, pval = np.nan, np.nan
                        pvals_all.append(pval)
                        meta_all.append(
                            dict(Tissue=tissue, Subgroup=subgroup, Test=test,
                                 Env=cat, OR=or_, PValue=pval)
                        )

        # FDR across all comparisons within tissue
        valid_mask = [not np.isnan(p) for p in pvals_all]
        fdr_arr = np.full(len(pvals_all), np.nan)
        valid_pvals = [p for p, v in zip(pvals_all, valid_mask) if v]
        if valid_pvals:
            _, fdr_valid = fdrcorrection(valid_pvals)
            vi = 0
            for i, v in enumerate(valid_mask):
                if v:
                    fdr_arr[i] = fdr_valid[vi]
                    vi += 1

        for meta, fdr in zip(meta_all, fdr_arr):
            meta["FDR"] = fdr
            records.append(meta)

    return pd.DataFrame(records)[
        ["Tissue", "Subgroup", "OR", "PValue", "FDR", "Test", "Env"]
    ]


def main():
    df = calculate_enrichment()
    out = here(
        "environmental-analysis/all_individuals/tissue_compare/"
        "non_heritable_subgroup/_m/subgroup_enrichment_analysis.tsv"
    )
    df.to_csv(out, sep="\t", index=False)
    print(f"Saved {len(df)} rows to {out}")
    session_info.show()


if __name__ == "__main__":
    main()
