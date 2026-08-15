#!/usr/bin/env python3
"""Boundary-hit exclusion and overlap-weighted balance sensitivities."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import ANALYSIS_SCHEMA_VERSION, write_tsv  # noqa: E402
from _lib.stats_utils import standardized_mean_difference, zscore  # noqa: E402


def _load_comparison_module():
    path = Path(__file__).with_name("05_fit_predictor_comparison.py")
    spec = importlib.util.spec_from_file_location("predictor_compare_followup", path)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise SystemExit(f"Cannot load {path}")
    spec.loader.exec_module(module)
    return module


COMPARE = _load_comparison_module()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", required=True)
    parser.add_argument("--primary-bridge", required=True)
    parser.add_argument("--aligned-bridge", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def boundary_models(bridge: pd.DataFrame, config: dict, definition: str) -> tuple[list[dict], dict]:
    common, covariates = COMPARE.prepare_common_set(bridge, config)
    boundary = _as_bool(common["h2_upper_boundary_hit"])
    use = common.loc[~boundary].copy()
    use["legacy_z"] = zscore(use["legacy_h2_unscaled"])
    use["calibrated_z"] = zscore(use["calibrated_h2"])
    if len(use) < 100:
        raise SystemExit(f"Too few boundary-excluded VMRs for {definition}: {len(use)}")
    rows = COMPARE.fit_model_comparison(use, covariates)
    for row in rows:
        row["vmr_definition"] = definition
        row["sensitivity"] = "exclude_calibrated_upper_boundary_hits"
    qc = {
        "vmr_definition": definition,
        "n_common": int(len(common)),
        "n_boundary_hits_excluded": int(boundary.sum()),
        "n_boundary_excluded": int(len(use)),
        "selected_covariates": ",".join(covariates),
        "all_models_converged": all(bool(row.get("converged", False)) for row in rows),
    }
    return rows, qc


def _weighted_smd(frame: pd.DataFrame, exposure: pd.Series, weights: pd.Series, column: str) -> float:
    values = pd.to_numeric(frame[column], errors="coerce")
    result = []
    variances = []
    for state in [True, False]:
        mask = exposure.eq(state) & values.notna() & weights.notna()
        x = values[mask].to_numpy(float)
        w = weights[mask].to_numpy(float)
        mean = np.average(x, weights=w)
        variance = np.average((x - mean) ** 2, weights=w)
        result.append(mean)
        variances.append(variance)
    pooled = np.sqrt(np.mean(variances))
    return 0.0 if not np.isfinite(pooled) or pooled == 0 else float((result[0] - result[1]) / pooled)


def overlap_weighted(frame: pd.DataFrame, covariates: list[str]) -> tuple[dict, pd.DataFrame]:
    d = frame.copy()
    exposure = _as_bool(d["calibrated_positive_signal"])
    numeric = [column for column in covariates if not column.startswith("annotation_")]
    x = pd.DataFrame(index=d.index)
    for column in numeric:
        x[column] = zscore(d[column])
    categories = pd.get_dummies(d["broad_genomic_annotation"], prefix="annotation", drop_first=True, dtype=float)
    x = pd.concat([x, categories], axis=1)
    complete = x.notna().all(axis=1) & d["proportion_cpgs_with_sig_meqtl"].notna()
    d = d.loc[complete].copy()
    x = x.loc[complete]
    exposure = exposure.loc[complete]
    design = sm.add_constant(x, has_constant="add")
    try:
        propensity_fit = sm.GLM(exposure.astype(int), design, family=sm.families.Binomial()).fit(maxiter=250)
        propensity = propensity_fit.predict(design)
        method = "binomial_glm"
    except Exception:
        propensity_fit = sm.GLM(exposure.astype(int), design, family=sm.families.Binomial()).fit_regularized(
            alpha=1e-6, L1_wt=0.0, maxiter=1000
        )
        propensity = propensity_fit.predict(design)
        method = "ridge_binomial_glm"
    propensity = pd.Series(np.clip(propensity, 1e-6, 1 - 1e-6), index=d.index)
    weights = pd.Series(np.where(exposure, 1 - propensity, propensity), index=d.index)
    outcome = pd.to_numeric(d["proportion_cpgs_with_sig_meqtl"], errors="coerce")
    outcome_design = sm.add_constant(exposure.astype(int).rename("calibrated_positive_signal"), has_constant="add")
    result = sm.WLS(outcome, outcome_design, weights=weights).fit(cov_type="HC3")

    balance_rows = []
    for column in numeric:
        balance_rows.append({
            "covariate": column,
            "smd_unweighted": standardized_mean_difference(
                d.loc[exposure, column], d.loc[~exposure, column]
            ),
            "smd_overlap_weighted": _weighted_smd(d, exposure, weights, column),
        })
    for category in categories.columns:
        d[category] = categories.loc[d.index, category]
        balance_rows.append({
            "covariate": category,
            "smd_unweighted": standardized_mean_difference(
                d.loc[exposure, category], d.loc[~exposure, category]
            ),
            "smd_overlap_weighted": _weighted_smd(d, exposure, weights, category),
        })
    balance = pd.DataFrame(balance_rows)
    max_smd = float(balance["smd_overlap_weighted"].abs().max()) if len(balance) else 0.0
    treated_weight = weights[exposure]
    control_weight = weights[~exposure]
    return ({
        "analysis": "calibrated_positive_overlap_weighted",
        "n_vmrs": int(len(d)),
        "n_positive": int(exposure.sum()),
        "n_nonpositive": int((~exposure).sum()),
        "effective_n_positive": float(treated_weight.sum() ** 2 / (treated_weight.pow(2).sum())),
        "effective_n_nonpositive": float(control_weight.sum() ** 2 / (control_weight.pow(2).sum())),
        "mean_difference": float(result.params["calibrated_positive_signal"]),
        "se": float(result.bse["calibrated_positive_signal"]),
        "ci_lower": float(result.conf_int().loc["calibrated_positive_signal", 0]),
        "ci_upper": float(result.conf_int().loc["calibrated_positive_signal", 1]),
        "pvalue": float(result.pvalues["calibrated_positive_signal"]),
        "max_abs_smd_after": max_smd,
        "propensity_method": method,
        "weighting_scheme": "overlap_weights",
        "covariance": "HC3_WLS",
    }, balance)


def main() -> None:
    args = parse_args()
    config = COMPARE._load_config(args.config)
    region = args.region.lower()
    primary = pd.read_csv(args.primary_bridge, sep="\t")
    aligned = pd.read_csv(args.aligned_bridge, sep="\t")
    all_rows = []
    qc_rows = []
    for frame, definition in [(primary, "aa_exact_intersection"), (aligned, "all_individual_coordinate_aligned")]:
        rows, qc = boundary_models(frame, config, definition)
        all_rows.extend(rows)
        qc_rows.append(qc)

    primary_common, primary_covariates = COMPARE.prepare_common_set(primary, config)
    weighted, balance = overlap_weighted(primary_common, primary_covariates)
    weighted["vmr_definition"] = "aa_exact_intersection"
    weighted["balance_pass"] = weighted["max_abs_smd_after"] <= float(config["matching"]["balance_smd_max"])

    for row in all_rows:
        row["region"] = region
        row["analysis_schema_version"] = ANALYSIS_SCHEMA_VERSION
    for row in qc_rows:
        row["region"] = region
    weighted["region"] = region
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_rows).to_csv(output / "boundary_excluded_models.tsv", sep="\t", index=False)
    write_tsv(output / "boundary_excluded_qc.tsv", qc_rows)
    write_tsv(output / "overlap_weighted_burden.tsv", [weighted])
    balance.insert(0, "region", region)
    balance.to_csv(output / "overlap_weighted_balance.tsv", sep="\t", index=False)
    print(f"Wrote boundary and overlap-weighted sensitivities under {output}")


if __name__ == "__main__":
    main()

