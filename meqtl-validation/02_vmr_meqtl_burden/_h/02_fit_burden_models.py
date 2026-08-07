#!/usr/bin/env python3
"""Primary VMR meQTL-burden ~ continuous predictability models + matching."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True)
    p.add_argument("--burden-tsv", default="")
    p.add_argument("--seed", type=int, default=20260722)
    return p.parse_args()


def _z(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x * 0.0
    return (x - x.mean()) / sd


def fit_models(df: pd.DataFrame) -> list[dict]:
    rows = []
    need = ["proportion_cpgs_with_sig_meqtl", "local_predictability", "n_tested_cpgs",
            "n_cpgs_with_sig_meqtl"]
    d = df.dropna(subset=need).copy()
    d = d[d["n_tested_cpgs"] > 0]
    if d.empty:
        return [{"model": "none", "error": "no complete cases"}]

    d["n_ns"] = d["n_tested_cpgs"] - d["n_cpgs_with_sig_meqtl"]
    endog = d[["n_cpgs_with_sig_meqtl", "n_ns"]]

    specs: list[tuple[str, list[str]]] = [
        ("unadjusted", ["local_predictability"]),
        ("adjusted_minimal", ["local_predictability", "n_tested_cpgs"]),
    ]
    tech = []
    for c in ["average_cpg_coverage", "mean_cpg_variance", "length", "umap_k24_mean",
              "line_l1_frac", "segdup_frac"]:
        if c in d.columns and d[c].notna().sum() >= max(50, int(0.5 * len(d))):
            tech.append(c)
    if tech:
        specs.append(("adjusted_technical", ["local_predictability", "n_tested_cpgs"] + tech))

    for name, cols in specs:
        try:
            use = d.dropna(subset=cols).copy()
            if use.empty:
                rows.append({"model": name, "error": "no complete cases for covariates"})
                continue
            endog_u = use[["n_cpgs_with_sig_meqtl", "n_ns"]]
            exog = use[cols].apply(_z)
            exog = sm.add_constant(exog, has_constant="add")
            res = sm.GLM(endog_u, exog, family=sm.families.Binomial()).fit()
            rows.append({
                "model": name,
                "n_vmrs": len(use),
                "covariates": ",".join(cols),
                "coef_predictability": float(res.params.get("local_predictability", np.nan)),
                "se_predictability": float(res.bse.get("local_predictability", np.nan)),
                "pval_predictability": float(res.pvalues.get("local_predictability", np.nan)),
                "converged": bool(res.converged),
            })
        except Exception as exc:  # noqa: BLE001
            rows.append({"model": name, "error": str(exc)})

    try:
        fit = smf.ols(
            "proportion_cpgs_with_sig_meqtl ~ local_predictability + np.log1p(n_tested_cpgs)",
            data=d,
        ).fit()
        rows.append({
            "model": "ols_proportion_secondary",
            "n_vmrs": len(d),
            "covariates": "local_predictability,log1p(n_tested_cpgs)",
            "coef_predictability": float(fit.params["local_predictability"]),
            "se_predictability": float(fit.bse["local_predictability"]),
            "pval_predictability": float(fit.pvalues["local_predictability"]),
            "r_squared": float(fit.rsquared),
        })
    except Exception as exc:  # noqa: BLE001
        rows.append({"model": "ols_proportion_secondary", "error": str(exc)})
    return rows


def matched_contrast(df: pd.DataFrame, seed: int) -> list[dict]:
    match_cols = ["n_tested_cpgs"]
    for c in ["average_cpg_coverage", "mean_cpg_variance", "length", "umap_k24_mean"]:
        if c in df.columns and df[c].notna().sum() >= 50:
            match_cols.append(c)

    need = ["local_predictability", "proportion_cpgs_with_sig_meqtl"] + match_cols
    d = df.dropna(subset=need).copy()
    if len(d) < 20:
        return [{"analysis": "matched", "error": "too few VMRs"}]

    for c in match_cols:
        d[f"z_{c}"] = _z(d[c]).astype(float)

    q_hi = d["local_predictability"].quantile(0.8)
    q_lo = d["local_predictability"].quantile(0.2)
    hi = d[d["local_predictability"] >= q_hi].copy()
    lo = d[d["local_predictability"] <= q_lo].copy()
    lo = lo.sample(frac=1, random_state=seed)
    used: set = set()
    pairs = []
    zcols = [f"z_{c}" for c in match_cols]
    hi_z = hi[zcols].to_numpy(dtype=float)
    lo_z = lo[zcols].to_numpy(dtype=float)
    lo_idx = lo.index.to_numpy()
    lo_prop = lo["proportion_cpgs_with_sig_meqtl"].to_numpy(dtype=float)
    used_mask = np.zeros(len(lo), dtype=bool)
    for i in range(len(hi)):
        avail = ~used_mask
        if not avail.any():
            break
        diffs = lo_z[avail] - hi_z[i]
        dist = np.sqrt((diffs ** 2).sum(axis=1))
        avail_pos = np.flatnonzero(avail)
        j_local = int(avail_pos[int(np.argmin(dist))])
        used_mask[j_local] = True
        pairs.append((float(hi.iloc[i]["proportion_cpgs_with_sig_meqtl"]), float(lo_prop[j_local])))
    if not pairs:
        return [{"analysis": "matched", "error": "no pairs"}]

    hi_p, lo_p = zip(*pairs)
    hi_p = np.asarray(hi_p, dtype=float)
    lo_p = np.asarray(lo_p, dtype=float)
    diff = hi_p - lo_p
    rng = np.random.default_rng(seed)
    null = []
    pooled = np.concatenate([hi_p, lo_p])
    n = len(hi_p)
    for _ in range(2000):
        rng.shuffle(pooled)
        null.append(pooled[:n].mean() - pooled[n:n + n].mean())
    null = np.asarray(null)
    p_perm = (np.sum(np.abs(null) >= abs(diff.mean())) + 1) / (len(null) + 1)
    return [{
        "analysis": "matched_high_vs_low_predictability",
        "n_pairs": n,
        "mean_proportion_high": float(hi_p.mean()),
        "mean_proportion_low": float(lo_p.mean()),
        "mean_difference": float(diff.mean()),
        "permutation_pvalue": float(p_perm),
        "matching_variables": ",".join(match_cols),
    }]


def main() -> None:
    args = parse_args()
    root = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
    burden_path = Path(args.burden_tsv) if args.burden_tsv else (
        root / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m" / args.region / "vmr_meqtl_burden.tsv.gz"
    )
    if not burden_path.exists():
        raise SystemExit(f"Missing {burden_path}; run 01_aggregate_vmr_burden.py first")
    df = pd.read_csv(burden_path, sep="\t")
    outdir = burden_path.parent
    model_rows = fit_models(df)
    for r in model_rows:
        r["region"] = args.region
    write_tsv(outdir / "burden_model_results.tsv", model_rows)
    match_rows = matched_contrast(df, args.seed)
    for r in match_rows:
        r["region"] = args.region
    write_tsv(outdir / "matched_analysis_results.tsv", match_rows)
    print(f"Wrote model and matched results under {outdir}")


if __name__ == "__main__":
    main()
