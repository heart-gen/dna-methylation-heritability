#!/usr/bin/env python3
"""06_partitioned_heritability -- map the continuous annotation onto reference SNPs.

Usage:
    python 03_make_annot.py --run-id <id> --chrom 22

Ported from the legacy make_annot_continuous.py, reduced to the single
continuous-annotation case. The legacy caller passed five quintile BEDs and its
step_2.sh then stripped the value column with

    awk '{print $1,$2,$3}'

so the pipeline that called itself continuous emitted five binary indicators.
Here there is exactly one annotation column and it holds the score.

Output is a thin-annot file (annotation columns only, no CHR/BP/SNP/CM), which
is what `ldsc.py --l2 --thin-annot` expects.
"""
from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from pybedtools import BedTool

ANNOT_NAME = "LOCAL_SNP_CONTRIBUTION_Z"


def repo_root() -> Path:
    root = os.environ.get("V2_REPO_ROOT")
    if root:
        return Path(root)
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").is_dir():
            return parent
    raise SystemExit("Could not locate repository root")


def read_bim(bimfile: Path) -> pd.DataFrame:
    return pd.read_csv(bimfile, sep=r"\s+", usecols=[0, 1, 2, 3],
                       names=["CHR", "SNP", "CM", "BP"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--chrom", required=True, type=int)
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
    arm = cfg["ld_reference_arm"]
    ref = cfg["ld_references"][arm]
    merge_op = cfg["annotation"]["merge_op"]

    bimfile = Path(ref["bim_dir"]) / f"{ref['bim_prefix']}{args.chrom}.bim"
    if not bimfile.exists():
        raise SystemExit(f"Reference bim not found: {bimfile}")

    bed_path = run_dir / "annotation" / "annotation-hg19.bed"
    if not bed_path.exists():
        raise SystemExit(f"Annotation BED not found: {bed_path} "
                         "(run 02_liftover_annotation.py first)")

    # Overlapping VMRs are merged with the configured operator so a SNP under
    # two intervals gets one value on the score's own scale. Summing would make
    # the annotation depend on VMR density rather than on local genetic control.
    bed = BedTool(str(bed_path)).sort().merge(c=4, o=merge_op)

    df_bim = read_bim(bimfile)
    iter_bim = [["chr" + str(c), int(bp) - 1, int(bp)]
                for c, bp in np.array(df_bim[["CHR", "BP"]])]
    bimbed = BedTool(iter_bim)

    # o="mean" over the merged annotation: a SNP inside no VMR scores 0, which
    # is the correct "no local genetic control measured here" value for a
    # z-scored annotation centred on the VMR universe.
    mapped = bimbed.map(bed, c=4, o="mean")

    values = []
    for row in mapped:
        val = row[3]
        values.append(0.0 if val in (".", "", None) else float(val))

    annot = pd.DataFrame({ANNOT_NAME: values})
    if len(annot) != len(df_bim):
        raise SystemExit(
            f"Annotation length {len(annot)} != bim length {len(df_bim)} on "
            f"chr{args.chrom}; the SNP order would be misaligned.")

    out = run_dir / "ldscores" / f"annot.{args.chrom}.annot.gz"
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt") as f:
        annot.to_csv(f, sep="\t", index=False)

    nonzero = int((annot[ANNOT_NAME] != 0).sum())
    print(f"[06] chr{args.chrom}: {len(annot)} SNPs, {nonzero} in-annotation "
          f"({nonzero / max(len(annot), 1):.2%}), sd {annot[ANNOT_NAME].std():.4f}")


if __name__ == "__main__":
    main()
