#!/usr/bin/env python3
"""06_partitioned_heritability -- run S-LDSC for one trait against the model.

Usage:
    python 06_partition_h2.py --run-id <id> --trait scz

Reports the three metrics named in the legacy interpreting_sldsc_results.md --
enrichment, enrichment p, and the tau coefficient z-score -- once per annotation.

The model is baselineLD plus TWO annotations of ours: the binary VMR_TESTED
membership indicator and the continuous LOCAL_SNP_CONTRIBUTION_Z score. tau on
the score is therefore conditional on membership, which is what makes it the
within-VMR gradient the module claims to test; see annotations.py. tau on
VMR_TESTED answers a different question -- are tested VMRs enriched at all --
and stage 07 reports it descriptively, outside the frozen FDR family.

For a CONTINUOUS annotation the tau z-score is the primary statistic:
"Enrichment" from LDSC is the ratio of the heritability share to the share of the
annotation's total value, and for a signed score that denominator is a signed
sum, so the ratio is not interpretable. Each emitted row carries
`enrichment_interpretable` so a reader does not have to remember which is which.
"""
from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotations import (ANNOT_COLUMNS, ANNOT_ROLE, ENRICHMENT_INTERPRETABLE,
                         PRIMARY_ANNOT, check_annot_header, find_results_row)


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


def read_m_5_50(ld_prefix: Path, chroms=range(1, 23)) -> dict:
    """Total M_5_50 per annotation, summed over chromosomes.

    Diagnostic rather than decorative. M_5_50 is the column sum of the
    annotation over MAF 5-50% SNPs, which is the denominator LDSC divides by to
    form Prop._SNPs and hence Enrichment. For VMR_TESTED it is a SNP count; for
    the signed score it is a signed sum, and seeing the two side by side in the
    metrics table is the clearest available statement of why enrichment is
    interpretable for one and not the other.
    """
    totals = [0.0] * len(ANNOT_COLUMNS)
    seen = 0
    for chrom in chroms:
        path = Path(f"{ld_prefix}{chrom}.l2.M_5_50")
        if not path.exists():
            continue
        vals = [float(v) for v in path.read_text().split()]
        if len(vals) != len(ANNOT_COLUMNS):
            raise SystemExit(
                f"{path} has {len(vals)} entries, expected "
                f"{len(ANNOT_COLUMNS)} (one per annotation). The LD scores were "
                "computed from a different annotation set than this stage "
                "expects; recompute stage 05.")
        totals = [t + v for t, v in zip(totals, vals)]
        seen += 1
    if seen == 0:
        return {}
    return dict(zip(ANNOT_COLUMNS, totals))


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

    # Fail closed before launching the regression: --overlap-annot reads these
    # .annot.gz files to build the overlap matrix, and a one-annotation file here
    # would produce a tau that is not conditional on membership -- the exact
    # defect the two-annotation model exists to remove.
    with gzip.open(f"{ld_prefix}1.annot.gz", "rt") as fh:
        check_annot_header(fh.readline().rstrip("\n").split("\t"),
                           f"{ld_prefix}1.annot.gz")

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
    log_fields = parse_log(Path(f"{out_prefix}.log"))
    m_5_50 = read_m_5_50(ld_prefix)

    recs = []
    for name in ANNOT_COLUMNS:
        # One row per annotation, identified by name. Position is never used: our
        # annotations come after baselineLD's, but a baselineLD version change
        # must not be able to shift which row is read.
        row = find_results_row(df, name)
        rec = {
            "trait": args.trait,
            "annotation": row.iloc[0],
            "annotation_name": name,
            "annotation_role": ANNOT_ROLE[name],
            "is_primary_hypothesis": name == PRIMARY_ANNOT,
            "enrichment_interpretable": ENRICHMENT_INTERPRETABLE[name],
            "m_5_50": m_5_50.get(name),
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
        rec.update(log_fields)
        recs.append(rec)

    pd.DataFrame(recs).to_csv(out_dir / f"{args.trait}.metrics.tsv",
                              sep="\t", index=False)
    for rec in recs:
        print(f"[06] {args.trait} [{rec['annotation_role']}]: "
              f"enrichment {rec['enrichment']}, p {rec['enrichment_p']}, "
              f"tau_z {rec['tau_z']}, M_5_50 {rec['m_5_50']}")


if __name__ == "__main__":
    main()
