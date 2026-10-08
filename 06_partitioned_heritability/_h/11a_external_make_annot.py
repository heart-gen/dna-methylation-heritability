#!/usr/bin/env python3
"""External-annotation benchmark: one binary annotation onto reference SNPs.

Usage:
    python 11a_external_make_annot.py --run-id <id> --chrom 22 --annotation NAME
    python 11a_external_make_annot.py --run-id <id> --chrom 22 --annotation NAME --relabel

Default mode writes `ldscores/<NAME>/annot.<chrom>.annot.gz`, a thin-annot file
with ONE column named NAME (1 if the reference SNP lies in an interval of
`annotation/<NAME>.hg19.bed`, else 0), on exactly the reference bim the VMR runs
used, and `annot.<chrom>.common.tsv` with the count, sum and sum of squares of
the annotation over SNPs with MAF > min MAF -- what stage 12 needs for
standardized tau.

`--relabel` runs after `ldsc.py --l2`: for a single-column annotation LDSC writes
the LD-score column as a bare `L2`, which the regression would label `L2_1`.
It is renamed to `<NAME>L2`, so the `.results` row is identified by name, as
stage 05a does for the VMR runs.

The overlap test is the VMR runs' test (03_make_annot.py): the SNP as a 1-bp
BED interval [BP-1, BP) intersected with the merged annotation intervals.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from pybedtools import BedTool


def repo_root() -> Path:
    root = os.environ.get("V2_REPO_ROOT")
    if root:
        return Path(root)
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").is_dir():
            return parent
    raise SystemExit("Could not locate repository root")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--chrom", required=True, type=int)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--relabel", action="store_true")
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "06_partitioned_heritability" / "_m" / "runs" / args.run_id
    out_dir = run_dir / "ldscores" / args.annotation
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / f"annot.{args.chrom}"

    if args.relabel:
        f = Path(f"{prefix}.l2.ldscore.gz")
        with gzip.open(f, "rt") as fh:
            lines = fh.readlines()
        header = lines[0].rstrip("\n").split("\t")
        want = f"{args.annotation}L2"
        if header[-1] == want:
            return
        if header[-1] != "L2":
            raise SystemExit(f"{f}: last column {header[-1]!r}, expected 'L2' or {want!r}")
        header[-1] = want
        lines[0] = "\t".join(header) + "\n"
        with gzip.open(f, "wt") as fh:
            fh.writelines(lines)
        return

    cfg = yaml.safe_load((run_dir / "code" / "config" / "partitioned_heritability.yml").read_text())
    ext = yaml.safe_load((run_dir / "code" / "config" / "sldsc_external_annotations.yml").read_text())
    ref = cfg["ld_references"][cfg["ld_reference_arm"]]
    bim = Path(ref["bim_dir"]) / f"{ref['bim_prefix']}{args.chrom}.bim"
    frq = Path(ref["frq_dir"]) / f"{ref['frq_prefix']}{args.chrom}.frq"
    bed = run_dir / "annotation" / f"{args.annotation}.hg19.bed"
    for p in (bim, frq, bed):
        if not p.exists():
            raise SystemExit(f"missing {p}")

    df_bim = pd.read_csv(bim, sep=r"\s+", usecols=[0, 1, 3], names=["CHR", "SNP", "BP"])
    snps = BedTool([["chr" + str(c), int(bp) - 1, int(bp)]
                    for c, bp in np.array(df_bim[["CHR", "BP"]])])
    ann = BedTool(str(bed)).sort().merge()
    # -c counts overlapping intervals; after merge it is 0 or 1. Order is the
    # bim order because bedtools intersect -c preserves the A file's order.
    counts = [int(r[-1]) for r in snps.intersect(ann, c=True)]
    if len(counts) != len(df_bim):
        raise SystemExit("annotation length != bim length; SNP order would misalign")
    a = (np.asarray(counts) > 0).astype(float)
    with gzip.open(f"{prefix}.annot.gz", "wt") as fh:
        pd.DataFrame({args.annotation: a}).to_csv(fh, sep="\t", index=False)

    fr = pd.read_csv(frq, sep=r"\s+")
    if not (fr["SNP"].astype(str).values == df_bim["SNP"].astype(str).values).all():
        raise SystemExit(f"{frq} is not row-aligned with {bim}")
    common = fr["MAF"].values > float(ext["standardized_tau_min_maf"])
    ac = a[common]
    pd.DataFrame({"chrom": [args.chrom], "n_snps": [len(a)], "n_in_annotation": [int(a.sum())],
                  "n_common": [int(common.sum())], "sum_common": [float(ac.sum())],
                  "sumsq_common": [float((ac ** 2).sum())]}).to_csv(
        f"{prefix}.common.tsv", sep="\t", index=False)
    print(f"[06x] {args.annotation} chr{args.chrom}: {int(a.sum())}/{len(a)} SNPs "
          f"({a.mean():.4%})", file=sys.stderr)


if __name__ == "__main__":
    main()
