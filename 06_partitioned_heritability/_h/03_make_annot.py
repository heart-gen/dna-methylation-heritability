#!/usr/bin/env python3
"""06_partitioned_heritability -- map the two annotations onto reference SNPs.

Usage:
    python 03_make_annot.py --run-id <id> --chrom 22

Ported from the legacy make_annot_continuous.py. The legacy caller passed five
quintile BEDs and its step_2.sh then stripped the value column with

    awk '{print $1,$2,$3}'

so the pipeline that called itself continuous emitted five binary indicators.

Two annotation columns are emitted here, in the order declared by
annotations.ANNOT_COLUMNS: a binary VMR_TESTED membership indicator and the
continuous LOCAL_SNP_CONTRIBUTION_Z score. Both go into the S-LDSC model, so tau
on the score is conditional on membership. `annotations.py` carries the full
argument for why one column cannot ask the module's question; the short version
is that the score is a standardized rank, so zero is a median-ranked VMR, and a
single column cannot distinguish "not a tested VMR" from "a tested VMR with
average local genetic control".

Output is a thin-annot file (annotation columns only, no CHR/BP/SNP/CM), which
is what `ldsc.py --l2 --thin-annot` expects.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotations import (ANNOT_COLUMNS, MEMBERSHIP_ANNOT, SCORE_ANNOT,
                         assert_score_source_column, check_annot_header)


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


def annot_columns(df_bim: pd.DataFrame, merged_bed: BedTool) -> pd.DataFrame:
    """Build the two annotation columns for one chromosome's reference SNPs.

    Both columns come from a SINGLE overlap test, so they cannot disagree about
    which SNPs are in the tested universe. `bedtools map -o mean` reports "."
    for a SNP under no interval, and that -- not the value of the score -- is
    what decides membership. A SNP inside a VMR whose score happens to be 0.0
    therefore gets membership 1 and score 0.0, while a SNP inside no VMR gets
    membership 0 and score 0.0; the one-column construction gave both the same
    0.0 and had no way to tell them apart.

    The score's 0.0 for a non-member SNP is now a placeholder rather than a
    claim: VMR_TESTED = 0 marks the row as carrying no measurement, and tau on
    the score is estimated conditional on that indicator.
    """
    iter_bim = [["chr" + str(c), int(bp) - 1, int(bp)]
                for c, bp in np.array(df_bim[["CHR", "BP"]])]
    bimbed = BedTool(iter_bim)
    mapped = bimbed.map(merged_bed, c=4, o="mean")

    membership = []
    scores = []
    for row in mapped:
        val = row[3]
        if val in (".", "", None):
            membership.append(0.0)
            scores.append(0.0)
        else:
            membership.append(1.0)
            scores.append(float(val))

    annot = pd.DataFrame({MEMBERSHIP_ANNOT: membership, SCORE_ANNOT: scores})
    # Reindex rather than trust the dict literal's order: the column order is
    # positional downstream (.l2.ldscore.gz, .l2.M_5_50), so it is pinned to the
    # declared contract here and nowhere else.
    annot = annot[list(ANNOT_COLUMNS)]
    if len(annot) != len(df_bim):
        raise SystemExit(
            f"Annotation length {len(annot)} != bim length {len(df_bim)}; "
            "the SNP order would be misaligned.")
    return annot


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
    assert_score_source_column(cfg["annotation"]["score_column"])

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
    # Merging is also what makes membership well defined: after it, the
    # intervals are disjoint, so "overlaps something" is a single test.
    bed = BedTool(str(bed_path)).sort().merge(c=4, o=merge_op)

    df_bim = read_bim(bimfile)
    annot = annot_columns(df_bim, bed)
    check_annot_header(annot.columns, "03_make_annot.py")

    out = run_dir / "ldscores" / f"annot.{args.chrom}.annot.gz"
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt") as f:
        annot.to_csv(f, sep="\t", index=False)

    n_member = int(annot[MEMBERSHIP_ANNOT].sum())
    # The count the one-column construction could not see: in-VMR SNPs whose
    # score is exactly zero were indistinguishable from the 99.4% of the genome
    # that is in no VMR at all.
    n_member_zero_score = int(((annot[MEMBERSHIP_ANNOT] == 1.0) &
                               (annot[SCORE_ANNOT] == 0.0)).sum())
    print(f"[06] chr{args.chrom}: {len(annot)} SNPs, {n_member} in-VMR "
          f"({n_member / max(len(annot), 1):.2%}), "
          f"{n_member_zero_score} in-VMR with score exactly 0, "
          f"score sd {annot[SCORE_ANNOT].std():.4f}")


if __name__ == "__main__":
    main()
