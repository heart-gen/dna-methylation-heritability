"""
Classify AA (BA) non-heritable VMRs by their EA (WA) heritability status,
using the all_individuals VMR set where both populations were estimated on
the same genomic coordinates.

The existing all_individuals enrichment analysis only uses VMRs where AA and
EA classifications are *concordant* — effectively Group B only. This script
removes that filter to expose all three subgroups:

  Group A — AA non-heritable, EA heritable
             (cryptically heritable candidate; possible LD-driven
              misclassification in AA due to shorter haplotype blocks)
  Group B — AA non-heritable, EA non-heritable
             (jointly non-heritable; strongest environmental-mediation
              candidate — these are the concordant VMRs already studied)
  Group C — AA non-heritable, EA low-prediction or absent after QC
             (EA model had insufficient predictive confidence; treat as
              population-specific to AA)

Outputs (per tissue):
  non_heritable_subgroups_{tissue}.tsv   — one row per VMR with subgroup
  subgroup_counts.tsv                    — summary counts across tissues
"""

import pandas as pd
import numpy as np
from pyhere import here
from functools import lru_cache
import session_info

H2_THRESH = 0.1
R2_THRESH = 0.3

# tissue directory name → display name mapping
TISSUE_DIRS = {"Caudate": "caudate", "DLPFC": "dlpfc", "Hippocampus": "hippocampus"}


def classify(h2, r2):
    if pd.isna(h2) or pd.isna(r2):
        return np.nan
    if r2 > R2_THRESH and h2 >= H2_THRESH:
        return "Heritable"
    elif r2 > R2_THRESH and h2 < H2_THRESH:
        return "Non-heritable"
    else:
        return "Low prediction"


def make_feature_id(df, chrom_col="chrom", start_col="start", end_col="end"):
    chrom = df[chrom_col].astype(str)
    if not chrom.str.startswith("chr").all():
        chrom = "chr" + chrom
    return chrom + "_" + df[start_col].astype(str) + "_" + df[end_col].astype(str)


@lru_cache()
def get_enet(tissue, race):
    """Load and classify elastic-net results for one tissue and one race (AA or EA)."""
    tdir = TISSUE_DIRS[tissue]
    fn = here(
        f"heritability/elastic_net_model/all_individuals/{tdir}/_m/"
        f"{tdir}_summary_elastic-net_{race}.tsv"
    )
    df = pd.read_csv(fn, sep="\t")
    df["feature_id"] = make_feature_id(df)
    df[f"h2_category_{race}"] = df.apply(
        lambda r: classify(r["h2_unscaled"], r["r_squared_cv"]), axis=1
    )
    return df[["feature_id", "chrom", "start", "end",
               "h2_unscaled", "r_squared_cv", f"h2_category_{race}"]].rename(
        columns={
            "h2_unscaled": f"h2_unscaled_{race}",
            "r_squared_cv": f"r_squared_cv_{race}",
        }
    )


def assign_groups(tissue):
    aa = get_enet(tissue, "AA")
    ea = get_enet(tissue, "EA")

    # Outer join on feature_id — keep all VMRs that exist in AA
    merged = aa.merge(ea, on="feature_id", how="left",
                      suffixes=("_AA", "_EA"))

    # Resolve duplicate coordinate columns after merge
    for col in ["chrom", "start", "end"]:
        col_aa = f"{col}_AA" if f"{col}_AA" in merged.columns else col
        if col_aa in merged.columns:
            merged[col] = merged[col_aa]
            merged.drop(columns=[c for c in [f"{col}_AA", f"{col}_EA"]
                                  if c in merged.columns], inplace=True)

    # Restrict to AA non-heritable VMRs
    aa_nh = merged[merged["h2_category_AA"] == "Non-heritable"].copy()

    def _group(row):
        ea_cat = row.get("h2_category_EA", np.nan)
        if pd.isna(ea_cat):
            return "C"  # VMR absent from EA after QC
        elif ea_cat == "Heritable":
            return "A"  # AA non-heritable, EA heritable
        elif ea_cat == "Non-heritable":
            return "B"  # jointly non-heritable
        else:
            return "C"  # EA low-prediction — insufficient confidence

    aa_nh["subgroup"] = aa_nh.apply(_group, axis=1)
    aa_nh["tissue"] = tissue

    print(f"{tissue}: {aa_nh['subgroup'].value_counts().to_dict()}")
    return aa_nh


def main():
    summary_rows = []

    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        df = assign_groups(tissue)
        out = here(
            f"environmental-analysis/all_individuals/tissue_compare/"
            f"non_heritable_subgroup/_m/non_heritable_subgroups_{tissue}.tsv"
        )
        df.to_csv(out, sep="\t", index=False)

        for grp, cnt in df["subgroup"].value_counts().items():
            summary_rows.append({"tissue": tissue, "subgroup": grp, "n_vmrs": int(cnt)})

    summary = pd.DataFrame(summary_rows)
    summary_out = here(
        "environmental-analysis/all_individuals/tissue_compare/"
        "non_heritable_subgroup/_m/subgroup_counts.tsv"
    )
    summary.to_csv(summary_out, sep="\t", index=False)
    print("\nSubgroup counts:\n", summary.to_string(index=False))

    session_info.show()


if __name__ == "__main__":
    main()
