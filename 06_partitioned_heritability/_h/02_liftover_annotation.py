#!/usr/bin/env python3
"""06_partitioned_heritability -- lift the annotation BED from hg38 to hg19.

Usage:
    python 02_liftover_annotation.py --run-id sldsc-AA-caudate-YYYYMMDD

Module 01/02 coordinates are hg38; every LDSC reference panel on Quest is hg19.
The liftover therefore happens on the annotation, once, rather than on the
reference panels.

Ported from the legacy region_heritability.py, keeping its liftover helpers and
discarding everything else in it (the r_squared_cv filter, the h2_unscaled
value column, and the quintile split are all banned by AGENTS.md 3).

Unlike the legacy script, a VMR whose interval does not survive liftover is
written to excluded/ with a reason rather than silently dropped: the annotation
denominator has to be reconstructable.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd
import yaml
from pyliftover import LiftOver

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def repo_root() -> Path:
    root = os.environ.get("V2_REPO_ROOT")
    if root:
        return Path(root)
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / ".git").is_dir():
            return parent
    raise SystemExit("Could not locate repository root")


def lift(lo: LiftOver, chrom: str, pos: int):
    """Return the lifted position, or None when the coordinate does not map.

    pyliftover is 0-based; BED is 0-based half-open, so start and end are both
    converted directly.
    """
    try:
        hit = lo.convert_coordinate(chrom, int(pos))
    except Exception:
        return None
    if not hit:
        return None
    return hit[0][0], int(hit[0][1])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "06_partitioned_heritability" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")

    # Prefer the run's config snapshot over the live working tree; see
    # run_config.py for why.
    _snap = run_dir / "code" / "config" / "partitioned_heritability.yml"
    _src = _snap if _snap.exists() else root / "config" / "partitioned_heritability.yml"
    cfg = yaml.safe_load(_src.read_text())
    chain = root / cfg["liftover_chain"]
    if not chain.exists():
        raise SystemExit(f"Liftover chain not found: {chain}")

    src = run_dir / "annotation" / "annotation-hg38.bed"
    df = pd.read_csv(src, sep="\t", header=None,
                     names=["chrom", "start", "end", "score"])
    n_in = len(df)

    lo = LiftOver(str(chain))

    rows, dropped = [], []
    for chrom, start, end, score in df.itertuples(index=False, name=None):
        a = lift(lo, chrom, start)
        b = lift(lo, chrom, end)
        if a is None or b is None:
            dropped.append((chrom, start, end, score, "unmapped"))
            continue
        # A VMR whose two ends land on different contigs has not been lifted as
        # an interval; keeping it would invent a region that exists in neither
        # build.
        if a[0] != b[0]:
            dropped.append((chrom, start, end, score, "split_across_contigs"))
            continue
        new_chrom = a[0]
        lo_pos, hi_pos = sorted((a[1], b[1]))
        if new_chrom not in AUTOSOMES:
            dropped.append((chrom, start, end, score, "non_autosomal_after_lift"))
            continue
        if lo_pos >= hi_pos:
            dropped.append((chrom, start, end, score, "degenerate_after_lift"))
            continue
        rows.append((new_chrom, lo_pos, hi_pos, score))

    out = pd.DataFrame(rows, columns=["chrom", "start", "end", "score"])
    if out.empty:
        raise SystemExit("Liftover produced no intervals")

    out["chrom"] = pd.Categorical(out["chrom"], categories=AUTOSOMES, ordered=True)
    out = out.sort_values(["chrom", "start"]).reset_index(drop=True)

    dst = run_dir / "annotation" / "annotation-hg19.bed"
    out.to_csv(dst, sep="\t", header=False, index=False)

    if dropped:
        pd.DataFrame(dropped,
                     columns=["chrom", "start", "end", "score", "reason"]).to_csv(
            run_dir / "excluded" / "liftover-dropped.tsv", sep="\t", index=False)

    frac = 1.0 - (len(out) / n_in) if n_in else 1.0
    max_missing = float(cfg["gates"]["max_annotation_missing_fraction"])
    stats = pd.DataFrame([{
        "n_hg38": n_in,
        "n_hg19": len(out),
        "n_dropped": len(dropped),
        "fraction_dropped": frac,
        "max_allowed_fraction_dropped": max_missing,
        "score_sd_hg19": float(out["score"].std()),
        "n_distinct_hg19": int(out["score"].nunique()),
    }])
    stats.to_csv(run_dir / "annotation" / "liftover-summary.tsv",
                 sep="\t", index=False)

    print(f"[06] liftover {n_in} -> {len(out)} intervals "
          f"({len(dropped)} dropped, {frac:.4%})")

    # Losing a large share of the annotation changes what the annotation IS, so
    # it fails here rather than propagating into an enrichment estimate.
    if frac > max_missing:
        raise SystemExit(
            f"Liftover dropped {frac:.2%} of intervals, above the configured "
            f"maximum {max_missing:.2%}. The hg19 annotation is not a faithful "
            f"image of the hg38 one.")


if __name__ == "__main__":
    main()
