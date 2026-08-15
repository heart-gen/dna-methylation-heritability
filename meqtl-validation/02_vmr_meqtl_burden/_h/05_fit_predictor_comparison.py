#!/usr/bin/env python3
"""Compare calibrated and legacy predictors against repaired VMR meQTL burden."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import ANALYSIS_SCHEMA_VERSION, load_yaml, write_tsv  # noqa: E402
from _lib.stats_utils import (  # noqa: E402
    greedy_binary_propensity_pairs,
    paired_randomization_pvalue,
    zscore,
)

COMPARISON_SCHEMA_VERSION = 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", required=True)
    parser.add_argument("--bridge", required=True)
    parser.add_argument("--config", default="predictor_comparison.yml")
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def _load_config(value: str) -> dict:
    if "/" not in value:
        return load_yaml(value)
    import yaml

    with Path(value).open() as handle:
        return yaml.safe_load(handle)


def _as_bool(series: pd.Series) -> pd.Series:
    """Parse booleans without treating the string 'False' as truthy."""
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def _annotation_indicators(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "genomic_annotation" in out.columns:
        annotation = out["genomic_annotation"].fillna("").astype(str)
    else:
        annotation = pd.Series("", index=out.index)
    out["annotation_promoter"] = annotation.str.contains("promoter", case=False).astype(int)
    out["annotation_cpg_island"] = annotation.str.contains(
        "cpg_island|cpg island", case=False, regex=True
    ).astype(int)
    out["annotation_gene_body"] = annotation.str.contains(
        "exon|intron", case=False, regex=True
    ).astype(int)
    out["broad_genomic_annotation"] = np.select(
        [
            out["annotation_promoter"].eq(1),
            out["annotation_cpg_island"].eq(1),
            out["annotation_gene_body"].eq(1),
        ],
        ["promoter", "cpg_island", "gene_body"],
        default="other",
    )
    return out


def prepare_common_set(frame: pd.DataFrame, config: dict) -> tuple[pd.DataFrame, list[str]]:
    required = [
        "common_complete_case", "n_tested_cpgs", "n_cpgs_with_sig_meqtl",
        "proportion_cpgs_with_sig_meqtl", "legacy_h2_unscaled", "calibrated_h2",
    ]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise SystemExit(f"Bridge is missing: {', '.join(missing)}")
    d = _annotation_indicators(frame)
    d = d[_as_bool(d["common_complete_case"])].copy()
    d = d[pd.to_numeric(d["n_tested_cpgs"], errors="coerce").gt(0)]
    d["n_nonsignificant_cpgs"] = (
        pd.to_numeric(d["n_tested_cpgs"], errors="coerce")
        - pd.to_numeric(d["n_cpgs_with_sig_meqtl"], errors="coerce")
    )
    if d["n_nonsignificant_cpgs"].lt(0).any():
        raise SystemExit("Significant CpG count exceeds the tested denominator")

    model_cfg = config["model"]
    minimum = float(model_cfg["min_covariate_complete_fraction"])
    candidates = list(model_cfg["candidate_covariates"]) + list(
        model_cfg["annotation_indicators"]
    )
    covariates = [
        column for column in candidates
        if column in d.columns and d[column].notna().mean() >= minimum
    ]
    complete_columns = [
        "n_tested_cpgs", "n_cpgs_with_sig_meqtl", "n_nonsignificant_cpgs",
        "proportion_cpgs_with_sig_meqtl", "legacy_h2_unscaled", "calibrated_h2",
    ] + covariates
    d = d.dropna(subset=complete_columns).copy()
    if len(d) < 100:
        raise SystemExit(f"Too few common complete VMRs after covariate QC: {len(d)}")
    d["legacy_z"] = zscore(d["legacy_h2_unscaled"])
    d["calibrated_z"] = zscore(d["calibrated_h2"])
    return d, covariates


def _fit_quasibinomial(
    frame: pd.DataFrame,
    *,
    model_name: str,
    predictors: list[str],
    covariates: list[str],
    primary_terms: dict[str, str],
) -> list[dict]:
    columns = predictors + covariates
    exog = frame[columns].apply(zscore)
    exog = sm.add_constant(exog, has_constant="add")
    endog = frame[["n_cpgs_with_sig_meqtl", "n_nonsignificant_cpgs"]]
    model = sm.GLM(endog, exog, family=sm.families.Binomial())
    pilot = model.fit(maxiter=250, tol=1e-8)
    dispersion = max(1.0, float(pilot.pearson_chi2 / max(pilot.df_resid, 1)))
    result = model.fit(scale=dispersion, cov_type="HC3", maxiter=250, tol=1e-8)
    rows = []
    for term, family in primary_terms.items():
        rows.append({
            "model": model_name,
            "test_family": family,
            "term": term,
            "estimate_per_sd": float(result.params[term]),
            "se": float(result.bse[term]),
            "ci_lower": float(result.params[term] - 1.96 * result.bse[term]),
            "ci_upper": float(result.params[term] + 1.96 * result.bse[term]),
            "pvalue": float(result.pvalues[term]),
            "n_vmrs": int(len(frame)),
            "covariates": ",".join(covariates),
            "dispersion_pearson": dispersion,
            "covariance": "HC3_quasibinomial",
            "converged": bool(result.converged),
        })
    return rows


def fit_model_comparison(frame: pd.DataFrame, covariates: list[str]) -> list[dict]:
    rows = []
    rows.extend(_fit_quasibinomial(
        frame,
        model_name="calibrated_adjusted_primary",
        predictors=["calibrated_z"],
        covariates=covariates,
        primary_terms={"calibrated_z": "calibrated_primary"},
    ))
    rows.extend(_fit_quasibinomial(
        frame,
        model_name="legacy_adjusted_comparator",
        predictors=["legacy_z"],
        covariates=covariates,
        primary_terms={"legacy_z": "legacy_comparator"},
    ))
    rows.extend(_fit_quasibinomial(
        frame,
        model_name="joint_adjusted_secondary",
        predictors=["calibrated_z", "legacy_z"],
        covariates=covariates,
        primary_terms={
            "calibrated_z": "joint_calibrated_secondary",
            "legacy_z": "joint_legacy_secondary",
        },
    ))
    return rows


def _match_one(
    frame: pd.DataFrame,
    *,
    analysis: str,
    exposure: str,
    covariates: list[str],
    config: dict,
    seed: int,
) -> tuple[dict, pd.DataFrame]:
    numeric = [column for column in covariates if not column.startswith("annotation_")]
    exact = ["broad_genomic_annotation"]
    pairs, balance, meta = greedy_binary_propensity_pairs(
        frame,
        exposure=exposure,
        outcome="proportion_cpgs_with_sig_meqtl",
        numeric_covariates=numeric,
        exact_covariates=exact,
        caliper_sd=float(config["matching"]["caliper_sd"]),
        seed=seed,
    )
    if pairs.empty:
        return ({
            "analysis": analysis,
            "exposure": exposure,
            "n_pairs": 0,
            "mean_burden_exposed": np.nan,
            "mean_burden_unexposed": np.nan,
            "mean_difference": np.nan,
            "permutation_pvalue": np.nan,
            "matching_variables": ",".join(numeric),
            "exact_matching_variables": ",".join(exact),
            "max_abs_smd_after": np.nan,
            "balance_pass": False,
            "permutation_scheme": "within_pair_label_swap",
            **meta,
        }, balance)
    differences = pairs["outcome_exposed"] - pairs["outcome_unexposed"]
    pvalue = paired_randomization_pvalue(
        differences,
        seed=seed,
        n_perm=int(config["matching"]["n_permutations"]),
    )
    balance_pass = bool(
        meta["max_abs_smd_after"] <= float(config["matching"]["balance_smd_max"])
    )
    result = {
        "analysis": analysis,
        "exposure": exposure,
        "n_pairs": int(len(pairs)),
        "mean_burden_exposed": float(pairs["outcome_exposed"].mean()),
        "mean_burden_unexposed": float(pairs["outcome_unexposed"].mean()),
        "mean_difference": float(differences.mean()),
        "permutation_pvalue": float(pvalue),
        "matching_variables": ",".join(numeric),
        "exact_matching_variables": ",".join(exact),
        "max_abs_smd_after": float(meta["max_abs_smd_after"]),
        "balance_pass": balance_pass,
        "permutation_scheme": "within_pair_label_swap",
    }
    balance = balance.copy()
    if not balance.empty:
        balance.insert(0, "analysis", analysis)
    return result, balance


def matched_comparisons(
    frame: pd.DataFrame, covariates: list[str], config: dict
) -> tuple[list[dict], pd.DataFrame]:
    seed = int(config["seed"])
    calibrated = frame.dropna(subset=["calibrated_positive_signal"]).copy()
    calibrated["calibrated_positive_signal"] = _as_bool(
        calibrated["calibrated_positive_signal"]
    )
    calibrated_result, calibrated_balance = _match_one(
        calibrated,
        analysis="calibrated_positive_vs_nonpositive",
        exposure="calibrated_positive_signal",
        covariates=covariates,
        config=config,
        seed=seed,
    )

    legacy = frame[frame["legacy_class"].isin(["high", "low"])].copy()
    legacy["legacy_high"] = _as_bool(legacy["legacy_high"])
    legacy_result, legacy_balance = _match_one(
        legacy,
        analysis="legacy_high_vs_low",
        exposure="legacy_high",
        covariates=covariates,
        config=config,
        seed=seed + 1,
    )
    balance = pd.concat(
        [calibrated_balance, legacy_balance], ignore_index=True, sort=False
    )
    return [calibrated_result, legacy_result], balance


def concordance_summaries(frame: pd.DataFrame) -> tuple[list[dict], list[dict], pd.DataFrame]:
    pairs = [
        ("legacy_h2_unscaled", "calibrated_h2"),
        ("legacy_h2_unscaled", "calibrated_rho2_oof"),
        ("legacy_h2_unscaled", "calibrated_r2_oof"),
        ("legacy_r_squared_cv", "calibrated_rho2_oof"),
        ("legacy_r_squared_cv", "calibrated_r2_oof"),
    ]
    concordance = []
    for left, right in pairs:
        use = frame[[left, right]].apply(pd.to_numeric, errors="coerce").dropna()
        rho, pvalue = spearmanr(use[left], use[right]) if len(use) >= 3 else (np.nan, np.nan)
        concordance.append({
            "legacy_metric": left,
            "calibrated_metric": right,
            "n_vmrs": int(len(use)),
            "spearman_rho": float(rho),
            "pvalue": float(pvalue),
        })

    old = _as_bool(frame["legacy_high"])
    new = _as_bool(frame["calibrated_positive_signal"])
    n11 = int((old & new).sum())
    n10 = int((old & ~new).sum())
    n01 = int((~old & new).sum())
    n00 = int((~old & ~new).sum())
    total = n11 + n10 + n01 + n00
    observed = (n11 + n00) / total if total else np.nan
    old_positive = (n11 + n10) / total if total else np.nan
    new_positive = (n11 + n01) / total if total else np.nan
    expected = (
        old_positive * new_positive + (1 - old_positive) * (1 - new_positive)
        if total else np.nan
    )
    kappa = (observed - expected) / (1 - expected) if total and expected < 1 else np.nan
    union = n11 + n10 + n01
    classification = [{
        "n_vmrs": total,
        "both_positive": n11,
        "legacy_only_positive": n10,
        "calibrated_only_positive": n01,
        "both_nonpositive": n00,
        "agreement": observed,
        "cohen_kappa": kappa,
        "jaccard_positive": n11 / union if union else np.nan,
    }]

    quantile = frame[["vmr_coord_id", "legacy_h2_unscaled", "calibrated_h2"]].copy()
    labels = [f"Q{i}" for i in range(1, 6)]
    quantile["legacy_quintile"] = pd.qcut(
        quantile["legacy_h2_unscaled"].rank(method="first"), 5, labels=labels
    )
    quantile["calibrated_quintile"] = pd.qcut(
        quantile["calibrated_h2"].rank(method="first"), 5, labels=labels
    )
    transition = (
        quantile.groupby(["legacy_quintile", "calibrated_quintile"], observed=False)
        .size()
        .rename("n_vmrs")
        .reset_index()
    )
    return concordance, classification, transition


def main() -> None:
    args = parse_args()
    config = _load_config(args.config)
    bridge = pd.read_csv(args.bridge, sep="\t")
    if (
        "predictor_comparison_schema_version" not in bridge.columns
        or not bridge["predictor_comparison_schema_version"].eq(
            COMPARISON_SCHEMA_VERSION
        ).all()
    ):
        raise SystemExit("Predictor bridge has an incompatible schema")
    common, covariates = prepare_common_set(bridge, config)
    model_rows = fit_model_comparison(common, covariates)
    matched_rows, balance = matched_comparisons(common, covariates, config)
    concordance, classification, transition = concordance_summaries(common)

    region = args.region.lower()
    for rows in [model_rows, matched_rows, concordance, classification]:
        for row in rows:
            row["region"] = region
            row["analysis_schema_version"] = ANALYSIS_SCHEMA_VERSION
            row["predictor_comparison_schema_version"] = COMPARISON_SCHEMA_VERSION
    transition.insert(0, "region", region)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(output_dir / "burden_model_comparison.tsv", model_rows)
    write_tsv(output_dir / "matched_comparison.tsv", matched_rows)
    write_tsv(output_dir / "predictor_concordance.tsv", concordance)
    write_tsv(output_dir / "classification_concordance.tsv", classification)
    transition.to_csv(output_dir / "quintile_transition.tsv", sep="\t", index=False)
    if not balance.empty:
        balance.insert(0, "region", region)
        balance.to_csv(output_dir / "matched_balance.tsv", sep="\t", index=False)
    write_tsv(output_dir / "comparison_model_qc.tsv", [{
        "region": region,
        "n_common_before_model_covariates": int(
            _as_bool(bridge["common_complete_case"]).sum()
        ),
        "n_common_model_complete": int(len(common)),
        "selected_covariates": ",".join(covariates),
        "all_models_converged": all(row.get("converged", False) for row in model_rows),
        "analysis_schema_version": ANALYSIS_SCHEMA_VERSION,
        "predictor_comparison_schema_version": COMPARISON_SCHEMA_VERSION,
    }])
    print(f"Wrote predictor comparison under {output_dir}")


if __name__ == "__main__":
    main()
