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
from _lib.io_utils import ANALYSIS_SCHEMA_VERSION, load_yaml, write_tsv  # noqa: E402
from _lib.stats_utils import (  # noqa: E402
    greedy_nearest_neighbor_pairs,
    paired_randomization_pvalue,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True)
    p.add_argument("--burden-tsv", default="")
    p.add_argument("--seed", type=int, default=20260722)
    p.add_argument("--require-complete-tech-join", action="store_true")
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
    if (d["n_ns"] < 0).any():
        raise ValueError("Significant CpG count exceeds tested CpG denominator")

    if "genomic_annotation" in d.columns:
        ann = d["genomic_annotation"].fillna("").astype(str)
        d["annotation_promoter"] = ann.str.contains("promoter", case=False).astype(int)
        d["annotation_cpg_island"] = ann.str.contains("cpg_island|cpg island", case=False).astype(int)
        d["annotation_gene_body"] = ann.str.contains("exon|intron", case=False).astype(int)

    specs: list[tuple[str, list[str]]] = [
        ("unadjusted", ["local_predictability"]),
        ("adjusted_minimal", ["local_predictability", "n_tested_cpgs"]),
    ]
    tech = []
    for c in [
        "average_cpg_coverage", "mean_cpg_variance", "vmr_mean_methylation",
        "length", "cpg_density", "mean_num_tested_snps_per_cpg",
        "umap_k24_mean", "line_l1_frac",
        "segdup_frac", "blacklist_frac", "annotation_promoter",
        "annotation_cpg_island", "annotation_gene_body",
    ]:
        if c in d.columns and d[c].notna().sum() >= max(50, int(0.5 * len(d))):
            tech.append(c)
    if tech:
        specs.append(("adjusted_prespecified", ["local_predictability", "n_tested_cpgs"] + tech))

    for name, cols in specs:
        try:
            use = d.dropna(subset=cols).copy()
            if use.empty:
                rows.append({"model": name, "error": "no complete cases for covariates"})
                continue
            endog_u = use[["n_cpgs_with_sig_meqtl", "n_ns"]]
            exog = use[cols].apply(_z)
            exog = sm.add_constant(exog, has_constant="add")
            # Quasi-binomial dispersion plus heteroskedasticity-robust covariance
            # prevents correlated CpGs within a VMR from yielding binomial-scale
            # standard errors. The VMR remains the independent row of inference.
            model = sm.GLM(endog_u, exog, family=sm.families.Binomial())
            pilot = model.fit(maxiter=250, tol=1e-8)
            dispersion = max(1.0, float(pilot.pearson_chi2 / max(pilot.df_resid, 1)))
            res = model.fit(scale=dispersion, cov_type="HC3", maxiter=250, tol=1e-8)
            rows.append({
                "model": name,
                "n_vmrs": len(use),
                "covariates": ",".join(cols),
                "coef_predictability": float(res.params.get("local_predictability", np.nan)),
                "se_predictability": float(res.bse.get("local_predictability", np.nan)),
                "pval_predictability": float(res.pvalues.get("local_predictability", np.nan)),
                "converged": bool(res.converged),
                "dispersion_pearson": float(res.scale),
                "covariance": "HC3_quasibinomial",
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


def matched_contrast(df: pd.DataFrame, seed: int) -> tuple[list[dict], pd.DataFrame]:
    d = df.copy()
    if "genomic_annotation" in d.columns:
        ann = d["genomic_annotation"].fillna("").astype(str)
        d["broad_genomic_annotation"] = np.select(
            [ann.str.contains("promoter", case=False), ann.str.contains("cpg_island", case=False), ann.str.contains("exon|intron", case=False)],
            ["promoter", "cpg_island", "gene_body"], default="other",
        )
    match_cols = ["n_tested_cpgs"]
    for c in [
        "average_cpg_coverage", "mean_cpg_variance", "vmr_mean_methylation",
        "length", "cpg_density", "mean_num_tested_snps_per_cpg",
        "umap_k24_mean", "line_l1_frac",
    ]:
        if c in d.columns and d[c].notna().sum() >= max(50, int(0.5 * len(d))):
            match_cols.append(c)
    exact = ["broad_genomic_annotation"] if "broad_genomic_annotation" in d.columns else []
    thresholds = load_yaml("analysis_thresholds.yml")["matching"]
    pairs, balance, meta = greedy_nearest_neighbor_pairs(
        d,
        exposure="local_predictability",
        outcome="proportion_cpgs_with_sig_meqtl",
        numeric_covariates=match_cols,
        exact_covariates=exact,
        caliper_sd=float(thresholds["caliper_sd"]),
        seed=seed,
    )
    if pairs.empty:
        return ([{"analysis": "matched", **meta}], balance)
    diff = pairs["outcome_high"].to_numpy() - pairs["outcome_low"].to_numpy()
    p_perm = paired_randomization_pvalue(diff, seed=seed, n_perm=10000)
    balance_pass = bool(meta["max_abs_smd_after"] <= float(thresholds["balance_smd_max"]))
    return ([{
        "analysis": "matched_high_vs_low_predictability",
        "n_pairs": int(len(pairs)),
        "mean_proportion_high": float(pairs["outcome_high"].mean()),
        "mean_proportion_low": float(pairs["outcome_low"].mean()),
        "mean_difference": float(diff.mean()),
        "permutation_pvalue": float(p_perm),
        "matching_variables": ",".join(match_cols),
        "exact_matching_variables": ",".join(exact),
        "caliper_sd": meta["caliper_sd"],
        "max_abs_smd_after": meta["max_abs_smd_after"],
        "balance_pass": balance_pass,
        "permutation_scheme": "within_pair_label_swap",
    }], balance)


def main() -> None:
    args = parse_args()
    root = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
    burden_path = Path(args.burden_tsv) if args.burden_tsv else (
        root / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m" / args.region / "vmr_meqtl_burden.tsv.gz"
    )
    if not burden_path.exists():
        raise SystemExit(f"Missing {burden_path}; run 01_aggregate_vmr_burden.py first")
    df = pd.read_csv(burden_path, sep="\t")
    if "analysis_schema_version" not in df or not df["analysis_schema_version"].eq(ANALYSIS_SCHEMA_VERSION).all():
        raise SystemExit("Burden table is stale; rerun 01_aggregate_vmr_burden.py with repair schema v2")
    if args.require_complete_tech_join and "tech_join_source" not in df.columns:
        raise SystemExit(
            "Primary AA burden modeling requires the repair-v2 technical join; "
            "run 07_repeat_mappability_sensitivity/_h/step_2_tech_joins.sh first"
        )
    outdir = burden_path.parent
    model_rows = fit_models(df)
    for r in model_rows:
        r["region"] = args.region
        r["analysis_schema_version"] = ANALYSIS_SCHEMA_VERSION
    write_tsv(outdir / "burden_model_results.tsv", model_rows)
    match_rows, balance = matched_contrast(df, args.seed)
    for r in match_rows:
        r["region"] = args.region
        r["analysis_schema_version"] = ANALYSIS_SCHEMA_VERSION
    write_tsv(outdir / "matched_analysis_results.tsv", match_rows)
    if not balance.empty:
        balance.insert(0, "region", args.region)
        balance.to_csv(outdir / "matched_analysis_balance.tsv", sep="\t", index=False)
    print(f"Wrote model and matched results under {outdir}")


if __name__ == "__main__":
    main()
