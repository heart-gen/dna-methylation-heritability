#!/usr/bin/env python3
"""Aggregate per-CpG coverage to VMR-level mean coverage tables."""

from __future__ import annotations

import argparse, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import PROJECT_ROOT, write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True, choices=["caudate", "dlpfc", "hippocampus"])
    p.add_argument("--population", default="AA")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    prep_root = (
        PROJECT_ROOT
        / "meqtl-validation"
        / "01_cpg_meqtl_mapping"
        / args.region
        / "_m"
        / "prepared"
    )
    prep = prep_root if args.population == "AA" else prep_root / args.population
    maps = sorted(prep.glob("cpg_vmr_map.chr*.tsv"))
    if not maps:
        raise SystemExit(f"No cpg_vmr_map files under {prep}")
    frames = []
    for path in maps:
        try:
            df = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            # Empty chromosomes are valid when the retained VMR universe has
            # no loci on that chromosome. They contribute no denominator rows.
            print(f"INFO: {path.name} is an empty-chromosome sentinel; skip")
            continue
        if "mean_coverage" not in df.columns:
            print(f"WARNING: {path.name} lacks mean_coverage; skip")
            continue
        frames.append(df)
    if not frames:
        raise SystemExit("No coverage-annotated cpg_vmr_map files yet; run step_6_coverage.sh first")
    d = pd.concat(frames, ignore_index=True)
    g = (
        d.groupby("vmr_id", as_index=False)
        .agg(
            n_cpgs=("phenotype_id", "count"),
            mean_cpg_coverage=("mean_coverage", "mean"),
            median_cpg_mean_coverage=("mean_coverage", "median"),
            mean_frac_samples_cov_ge_min=("fraction_samples_cov_ge_min", "mean"),
            min_mean_coverage=("mean_coverage", "min"),
        )
    )
    g["region"] = args.region
    out = prep / "vmr_mean_coverage.tsv"
    g.to_csv(out, sep="\t", index=False)
    write_tsv(
        prep / "vmr_mean_coverage_summary.tsv",
        [{
            "region": args.region,
            "population": args.population,
            "n_vmrs": len(g),
            "n_cpgs": int(d.shape[0]),
            "mean_of_vmr_mean_coverage": float(g["mean_cpg_coverage"].mean()),
            "path": str(out),
        }],
    )
    print(f"Wrote {out} ({len(g):,} VMRs)")


if __name__ == "__main__":
    main()
