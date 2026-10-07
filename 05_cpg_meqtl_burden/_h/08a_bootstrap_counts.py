#!/usr/bin/env python3
"""Per-draw VMR burden counts for the slope-inference run (05_cpg_meqtl_burden).

Usage:
    python _h/08a_bootstrap_counts.py --run-id cmb-AA-crossregion-YYYYMMDD

Turns each draw's per-CpG permutation p-values into the count 03_vmr_burden.R
models -- significant tested CpGs per VMR -- by the same three steps the
accepted cells used:
  1. p = pval_beta, falling back to pval_perm (02b_combine_meqtl.R);
  2. Storey q over every tested CpG of the REGION at once, via py_qvalue, the
     method the accepted cells record as fdr_method_used;
  3. significant = q <= mapping.fdr_threshold, counted over the cell's
     tested-cpg-membership.

Draw 0 is not a bootstrap draw. It is the accepted cell's own mapping output
pushed through this same code, and it must reproduce the sealed
n_cpgs_with_sig_meqtl for every VMR exactly. If it does not, this script is not
counting what 03_vmr_burden.R counted, and every draw would inherit the
discrepancy, so it stops.

The denominator stays the accepted cell's: a CpG tested in the cell but left
without a finite p-value in a draw (every cis variant fell under the MAF floor
in that resample) counts as not significant, and the per-draw audit records how
many there were.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import py_qvalue
import yaml

H_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(H_DIR))
from meqtl_covariates import read_manifest, repo_root  # noqa: E402


def to_pvalue(df: pd.DataFrame) -> pd.Series:
    p = df["pval_beta"].where(df["pval_beta"].notna(), df["pval_perm"])
    return p[np.isfinite(p)]


def storey(p: pd.Series) -> tuple[pd.Series, float]:
    res = py_qvalue.qvalue(p.to_numpy(dtype=float), lfdr_out=False)
    return pd.Series(np.asarray(res["qvalues"], dtype=float), index=p.index), float(res["pi0"])


def counts(sig_ids: set, member: pd.DataFrame) -> pd.Series:
    m = member.assign(sig=member["cpg_id"].isin(sig_ids))
    return m.groupby("vmr_id")["sig"].sum().astype(int)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    root = repo_root()
    mod = root / "05_cpg_meqtl_burden" / "_m" / "runs"
    run_dir = mod / args.run_id
    man = read_manifest(run_dir)
    cfg = yaml.safe_load((root / "config" / "meqtl_parameters.yml").read_text())
    thr = float(cfg["mapping"]["fdr_threshold"])
    per_task = int(cfg["cross_region_slope_inference"]["draws_per_task"])
    B = int(man["bootstrap_n"])
    blocks = [(a, min(a + per_task - 1, B)) for a in range(1, B + 1, per_task)]
    chroms = [f"chr{c}" for c in range(1, 23)]

    out_counts, audit = [], []
    for region in man["contrast"].split(","):
        cell = mod / man[f"upstream_cpg_meqtl_burden_{region}"] / "results"
        if (cell / "cpg-meqtl-results.tsv").exists():
            used = pd.read_csv(cell / "cpg-meqtl-results.tsv", sep="\t",
                               usecols=["cpg_id", "fdr_method_used"])
            if set(used["fdr_method_used"]) != {"storey_qvalue_py_qvalue"}:
                raise SystemExit(f"{region}: the cell used {set(used['fdr_method_used'])}, "
                                 "not py_qvalue; this script would count differently")
            tested = set(used["cpg_id"])
        member = pd.read_csv(cell / "tested-cpg-membership.tsv", sep="\t",
                             usecols=["cpg_id", "vmr_id"]).drop_duplicates()
        burden = pd.read_csv(cell / "vmr-meqtl-burden.tsv", sep="\t",
                             usecols=["vmr_id", "n_cpgs_with_sig_meqtl"]).set_index("vmr_id")
        member = member[member["vmr_id"].isin(burden.index)]

        # ---- draw 0: reproduce the sealed counts from the cell's own output
        obs = pd.concat([pd.read_csv(f, sep="\t", index_col=0)
                         for f in sorted((cell / "meqtl").glob("chr*.cis_qtl.tsv.gz"))])
        p0 = to_pvalue(obs)
        q0, pi0 = storey(p0)
        n0 = counts(set(q0.index[q0 <= thr]), member).reindex(burden.index, fill_value=0)
        bad = (n0 != burden["n_cpgs_with_sig_meqtl"]).sum()
        if bad:
            raise SystemExit(f"{region}: draw 0 disagrees with the sealed burden on {bad} "
                             "VMRs; the count is not the one 03_vmr_burden.R modelled")
        print(f"{region}: draw 0 reproduces all {len(n0)} sealed VMR counts (pi0 {pi0:.4f})")
        out_counts.append(pd.DataFrame({"region": region, "draw": 0,
                                        "vmr_id": n0.index, "n_sig": n0.values}))
        audit.append(dict(region=region, draw=0, n_cpgs_tested=len(tested),
                          n_cpgs_with_p=len(p0), n_cpgs_missing=len(tested - set(p0.index)),
                          n_sig_cpgs=int((q0 <= thr).sum()), pi0=pi0, n_samples=np.nan))

        # ---- the draws
        missing_files = []
        frames = []
        for c in chroms:
            for a, b in blocks:
                f = run_dir / "results" / "bootstrap" / region / f"{c}.draws{a:04d}-{b:04d}.tsv.gz"
                if not f.exists():
                    missing_files.append(f.name)
                    continue
                frames.append(pd.read_csv(f, sep="\t"))
        if missing_files:
            raise SystemExit(f"{region}: {len(missing_files)} bootstrap files missing, e.g. "
                             f"{missing_files[:3]}. Unexplained computational failures "
                             "stop a production run (AGENTS.md 9).")
        allp = pd.concat(frames)
        allp = allp[allp["cpg_id"].isin(tested)]
        for d, g in allp.groupby("draw"):
            p = to_pvalue(g.set_index("cpg_id"))
            q, pi0 = storey(p)
            n = counts(set(q.index[q <= thr]), member).reindex(burden.index, fill_value=0)
            out_counts.append(pd.DataFrame({"region": region, "draw": int(d),
                                            "vmr_id": n.index, "n_sig": n.values}))
            audit.append(dict(region=region, draw=int(d), n_cpgs_tested=len(tested),
                              n_cpgs_with_p=len(p),
                              n_cpgs_missing=len(tested - set(p.index)),
                              n_sig_cpgs=int((q <= thr).sum()), pi0=pi0,
                              n_samples=int(g["n_samples"].iloc[0])))
        got = sorted(allp["draw"].unique())
        if got != list(range(1, B + 1)):
            raise SystemExit(f"{region}: draws present {len(got)} of {B}")
        print(f"{region}: {B} draws counted")

    res = run_dir / "results"
    pd.concat(out_counts).to_csv(res / "bootstrap-burden-counts.tsv.gz", sep="\t",
                                 index=False, compression="gzip")
    pd.DataFrame(audit).to_csv(res / "bootstrap-draw-audit.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
