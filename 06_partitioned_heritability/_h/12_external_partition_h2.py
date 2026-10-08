#!/usr/bin/env python3
"""External-annotation benchmark: S-LDSC for one annotation x one trait.

Usage:
    python 12_external_partition_h2.py --run-id <id> --annotation NAME --trait scz

The regression is 06_partition_h2.py's command with one change: the custom
LD scores are this annotation's single column instead of the VMR runs' two.
Same baselineLD, weights, frq files, flags (--overlap-annot --thin-annot
--print-coefficients) and munged sumstats.

Reported, per the PI's list (2026-10-07):
    prop_snps, prop_h2 (+SE), enrichment (+SE, p), tau (+SE, z, two-sided p),
    standardized tau* = tau * sd(a) * M / h2g (Gazal et al. 2017), with sd(a)
    and M over reference SNPs with MAF above the configured floor, and h2g the
    regression's total observed-scale h2. tau*'s SE scales tau's SE by the same
    factor, treating sd(a), M and h2g as fixed.
"""
from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))


def repo_root() -> Path:
    root = os.environ.get("V2_REPO_ROOT")
    if root:
        return Path(root)
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").is_dir():
            return parent
    raise SystemExit("Could not locate repository root")


def parse_log(log_path: Path) -> dict:
    out = {"total_h2": None, "total_h2_se": None, "lambda_gc": None,
           "mean_chisq": None, "intercept": None}
    for line in log_path.read_text().splitlines():
        s = line.strip()
        try:
            if s.startswith("Total Observed scale h2:"):
                est, se = s.split(":", 1)[1].strip().split("(")
                out["total_h2"], out["total_h2_se"] = float(est), float(se.rstrip(")"))
            elif s.startswith("Lambda GC:"):
                out["lambda_gc"] = float(s.split(":", 1)[1])
            elif s.startswith("Mean Chi^2:"):
                out["mean_chisq"] = float(s.split(":", 1)[1])
            elif s.startswith("Intercept:"):
                out["intercept"] = float(s.split(":", 1)[1].strip().split("(")[0])
        except ValueError:
            pass
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--trait", required=True)
    args = ap.parse_args()

    root = repo_root()
    script_dir = Path(__file__).resolve().parent
    run_dir = root / "06_partitioned_heritability" / "_m" / "runs" / args.run_id
    cfg = yaml.safe_load((run_dir / "code" / "config" / "partitioned_heritability.yml").read_text())
    ref = cfg["ld_references"][cfg["ld_reference_arm"]]
    sumstats = run_dir / "sumstats" / f"{args.trait}.sumstats.gz"
    ld_dir = run_dir / "ldscores" / args.annotation
    ld_prefix = f"{ld_dir}/annot."
    for p in [sumstats] + [Path(f"{ld_prefix}{c}.l2.ldscore.gz") for c in range(1, 23)]:
        if not p.exists():
            raise SystemExit(f"missing {p}")

    out_dir = run_dir / "results" / "sldsc" / args.annotation
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / args.trait
    cmd = [sys.executable, str(script_dir / "ldsc_wrapper.py"), cfg["ldsc_dir"], "ldsc.py",
           "--h2", str(sumstats),
           "--ref-ld-chr", f"{ref['baseline_dir']}/{ref['baseline_prefix']},{ld_prefix}",
           "--w-ld-chr", f"{ref['weights_dir']}/{ref['weights_prefix']}",
           "--frqfile-chr", f"{ref['frq_dir']}/{ref['frq_prefix']}",
           "--overlap-annot", "--thin-annot", "--print-coefficients",
           "--out", str(out_prefix)]
    res = subprocess.run(cmd, cwd=str(run_dir))
    results_file = Path(f"{out_prefix}.results")
    if res.returncode != 0 or not results_file.exists():
        raise SystemExit(f"S-LDSC failed for {args.annotation} x {args.trait} (exit {res.returncode})")

    df = pd.read_csv(results_file, sep="\t")
    cat = df.iloc[:, 0].astype(str)
    hit = df[[c in (f"{args.annotation}L2_1", f"{args.annotation}L2") for c in cat]]
    if len(hit) != 1:
        raise SystemExit(f"expected one row named {args.annotation}L2_1, found "
                         f"{list(hit.iloc[:, 0])}; last rows {list(cat)[-3:]}")
    row = hit.iloc[0]
    log = parse_log(Path(f"{out_prefix}.log"))

    common = pd.concat([pd.read_csv(f"{ld_prefix}{c}.common.tsv", sep="\t") for c in range(1, 23)])
    m = float(common["n_common"].sum())
    mean = common["sum_common"].sum() / m
    sd = math.sqrt(max(common["sumsq_common"].sum() / m - mean ** 2, 0.0))
    m_5_50 = sum(float(Path(f"{ld_prefix}{c}.l2.M_5_50").read_text().split()[0])
                 for c in range(1, 23))
    tau, tau_se = float(row["Coefficient"]), float(row["Coefficient_std_error"])
    h2g = log["total_h2"]
    scale = (sd * m / h2g) if h2g and h2g != 0 else float("nan")
    z = tau / tau_se if tau_se > 0 else float("nan")

    rec = {
        "annotation": args.annotation, "trait": args.trait,
        "results_category": row.iloc[0],
        "n_snps_in_annotation_common": int(common["sum_common"].sum()),
        "m_common": int(m), "m_5_50": m_5_50,
        "prop_snps": row["Prop._SNPs"], "prop_h2": row["Prop._h2"],
        "prop_h2_se": row["Prop._h2_std_error"],
        "enrichment": row["Enrichment"], "enrichment_se": row["Enrichment_std_error"],
        "enrichment_p": row["Enrichment_p"],
        "tau": tau, "tau_se": tau_se, "tau_z": z,
        "tau_p_two_sided": math.erfc(abs(z) / math.sqrt(2)) if math.isfinite(z) else float("nan"),
        "annotation_sd_common": sd,
        "tau_star": tau * scale, "tau_star_se": tau_se * abs(scale),
    }
    rec.update(log)
    pd.DataFrame([rec]).to_csv(out_dir / f"{args.trait}.metrics.tsv", sep="\t", index=False)
    print(f"[06x] {args.annotation} x {args.trait}: prop_snps {rec['prop_snps']:.4f}, "
          f"prop_h2 {rec['prop_h2']:.4f}, enrichment {rec['enrichment']:.2f} "
          f"(p {rec['enrichment_p']:.3g}), tau* {rec['tau_star']:.3f} (z {z:.2f})")


if __name__ == "__main__":
    main()
