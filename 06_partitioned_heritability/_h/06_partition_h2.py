#!/usr/bin/env python3
"""06_partitioned_heritability -- run S-LDSC for one trait against the annotation.

Usage:
    python 06_partition_h2.py --run-id <id> --trait scz

Reports the three metrics named in the legacy interpreting_sldsc_results.md:
enrichment, enrichment p, and the tau coefficient z-score.

For a CONTINUOUS annotation, "Enrichment" from LDSC is the ratio of the share of
heritability to the share of the annotation's total value, and the coefficient
z-score is the test of whether the annotation adds signal over the baselineLD
model. The z-score is the primary statistic here: enrichment on a continuous
annotation is scale-dependent, while tau is conditional on baselineLD and is
what supports "heritability concentrates where local genetic control is high".
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

ANNOT_NAME = "LOCAL_SNP_CONTRIBUTION_Z"


def repo_root() -> Path:
    root = os.environ.get("V2_REPO_ROOT")
    if root:
        return Path(root)
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").is_dir():
            return parent
    raise SystemExit("Could not locate repository root")


def parse_log(log_path: Path) -> dict:
    """Pull total h2 and its SE out of the LDSC log.

    The .results table carries the per-annotation rows but not the model-level
    h2, and a cell whose total h2 is indistinguishable from zero cannot support
    an enrichment claim -- so the gate needs this.
    """
    out = {"total_h2": None, "total_h2_se": None, "lambda_gc": None,
           "mean_chisq": None, "intercept": None}
    if not log_path.exists():
        return out
    for line in log_path.read_text().splitlines():
        s = line.strip()
        if s.startswith("Total Observed scale h2:"):
            val = s.split(":", 1)[1].strip()
            try:
                est, se = val.split("(")
                out["total_h2"] = float(est.strip())
                out["total_h2_se"] = float(se.replace(")", "").strip())
            except ValueError:
                pass
        elif s.startswith("Lambda GC:"):
            try:
                out["lambda_gc"] = float(s.split(":", 1)[1])
            except ValueError:
                pass
        elif s.startswith("Mean Chi^2:"):
            try:
                out["mean_chisq"] = float(s.split(":", 1)[1])
            except ValueError:
                pass
        elif s.startswith("Intercept:"):
            val = s.split(":", 1)[1].strip()
            try:
                out["intercept"] = float(val.split("(")[0].strip())
            except ValueError:
                pass
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--trait", required=True)
    args = ap.parse_args()

    root = repo_root()
    script_dir = Path(__file__).resolve().parent
    run_dir = root / "06_partitioned_heritability" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")

    # Prefer the run's config snapshot over the live working tree; see
    # run_config.py for why.
    _snap = run_dir / "code" / "config" / "partitioned_heritability.yml"
    _src = _snap if _snap.exists() else root / "config" / "partitioned_heritability.yml"
    cfg = yaml.safe_load(_src.read_text())
    ref = cfg["ld_references"][cfg["ld_reference_arm"]]

    sumstats = run_dir / "sumstats" / f"{args.trait}.sumstats.gz"
    if not sumstats.exists():
        raise SystemExit(f"Munged sumstats not found: {sumstats}")

    ld_prefix = run_dir / "ldscores" / "annot."
    if not Path(f"{ld_prefix}1.l2.ldscore.gz").exists():
        raise SystemExit(f"Custom LD scores not found at {ld_prefix}* "
                         "(run 05_compute_ldscores.sh first)")

    out_dir = run_dir / "results" / "sldsc"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / args.trait

    cmd = [sys.executable, str(script_dir / "ldsc_wrapper.py"),
           cfg["ldsc_dir"], "ldsc.py",
           "--h2", str(sumstats),
           "--ref-ld-chr",
           f"{ref['baseline_dir']}/{ref['baseline_prefix']},{ld_prefix}",
           "--w-ld-chr", f"{ref['weights_dir']}/{ref['weights_prefix']}",
           "--frqfile-chr", f"{ref['frq_dir']}/{ref['frq_prefix']}",
           "--overlap-annot", "--thin-annot", "--print-coefficients",
           "--out", str(out_prefix)]

    print(f"[06] S-LDSC {args.trait}", flush=True)
    res = subprocess.run(cmd, cwd=str(run_dir))
    results_file = Path(f"{out_prefix}.results")
    if res.returncode != 0 or not results_file.exists():
        raise SystemExit(f"S-LDSC failed for {args.trait} (exit {res.returncode})")

    df = pd.read_csv(results_file, sep="\t")
    # The custom annotation is appended after the baselineLD columns, so it is
    # the last row; match on name rather than position so a baselineLD version
    # change cannot silently shift which row is read.
    hit = df[df.iloc[:, 0].astype(str).str.startswith(ANNOT_NAME)]
    if hit.empty:
        raise SystemExit(
            f"{args.trait}: no row named {ANNOT_NAME}* in {results_file}. "
            f"Found: {list(df.iloc[:, 0])[-3:]}")
    if len(hit) > 1:
        raise SystemExit(f"{args.trait}: {len(hit)} rows match {ANNOT_NAME}*")
    row = hit.iloc[0]

    rec = {
        "trait": args.trait,
        "annotation": row.iloc[0],
        "prop_snps": row.get("Prop._SNPs"),
        "prop_h2": row.get("Prop._h2"),
        "prop_h2_se": row.get("Prop._h2_std_error"),
        "enrichment": row.get("Enrichment"),
        "enrichment_se": row.get("Enrichment_std_error"),
        "enrichment_p": row.get("Enrichment_p"),
        "tau": row.get("Coefficient"),
        "tau_se": row.get("Coefficient_std_error"),
        "tau_z": row.get("Coefficient_z-score"),
    }
    rec.update(parse_log(Path(f"{out_prefix}.log")))
    pd.DataFrame([rec]).to_csv(out_dir / f"{args.trait}.metrics.tsv",
                               sep="\t", index=False)
    print(f"[06] {args.trait}: enrichment {rec['enrichment']}, "
          f"p {rec['enrichment_p']}, tau_z {rec['tau_z']}")


if __name__ == "__main__":
    main()
