#!/usr/bin/env python3
"""Build TensorQTL covariate matrix for primary CpG meQTL model."""

from __future__ import annotations

import argparse, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import load_paths, load_yaml, norm_region, read_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True, choices=["caudate", "dlpfc", "hippocampus"])
    p.add_argument("--population", default="AA")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    paths = load_paths()
    covars = load_yaml("covariates.yml")
    project = Path(paths["project_root"])
    outdir = project / "meqtl-validation" / "01_cpg_meqtl_mapping" / args.region / "_m" / "prepared"
    outdir.mkdir(parents=True, exist_ok=True)

    inclusion = read_tsv(
        project / "meqtl-validation" / "01_cpg_meqtl_mapping" / args.region / "_m"
        / "preflight" / "sample_inclusion_primary.tsv"
    )
    keep = [r["brnum"] for r in inclusion]

    phen = pd.read_csv(project / paths["phenotype_table"], sep="\t")
    phen["_region"] = phen["region"].map(norm_region)
    phen = phen[(phen["race"] == args.population) & (phen["_region"] == args.region)]
    phen = phen.set_index("brnum").loc[keep]

    cols = covars["primary_meqtl"]["required_phenotype_columns"]
    df = phen[cols].copy()
    # Encode sex / diagnosis as numeric for TensorQTL
    if df["sex"].dtype == object:
        df["sex"] = df["sex"].astype(str).str.upper().map({"M": 1, "MALE": 1, "F": 0, "FEMALE": 0}).astype(float)
    if df["primarydx"].dtype == object:
        df["primarydx"] = df["primarydx"].map({"Control": 0, "Schizo": 1}).astype(float)
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # TensorQTL covariates: samples x covariates
    # (index must equal phenotype BED sample columns)
    if df.isnull().any().any():
        bad = df.columns[df.isnull().any()].tolist()
        raise SystemExit(f"Missing covariate values in columns: {bad}")
    cov = df.copy()
    cov.index.name = "id"
    out = outdir / "covariates.txt"
    cov.to_csv(out, sep="\t")
    cov.reset_index().to_csv(outdir / "covariates_long.tsv", sep="\t", index=False)
    print(f"Wrote {out} shape={cov.shape} (samples x covariates)")


if __name__ == "__main__":
    main()
