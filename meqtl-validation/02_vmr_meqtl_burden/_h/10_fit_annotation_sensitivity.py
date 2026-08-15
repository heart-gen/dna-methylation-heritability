#!/usr/bin/env python3
"""Compare calibrated and legacy predictors across repeat/chromatin sensitivities."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import canonical_vmr_id, write_tsv  # noqa: E402
from _lib.stats_utils import (  # noqa: E402
    greedy_nearest_neighbor_pairs,
    paired_randomization_pvalue,
    zscore,
)

TISSUE = {"caudate": "Caudate", "dlpfc": "DLPFC", "hippocampus": "Hippocampus"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--repressive", required=True)
    parser.add_argument("--technical", required=True)
    parser.add_argument("--calibrated", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seed", type=int, default=20260809)
    parser.add_argument("--n-permutations", type=int, default=10000)
    return parser.parse_args()


def _coord(frame: pd.DataFrame, chrom: str, start: str, end: str) -> pd.Series:
    return pd.Series([
        canonical_vmr_id(c, s, e)
        for c, s, e in zip(frame[chrom], frame[start], frame[end])
    ], index=frame.index)


def load_analysis(args: argparse.Namespace) -> tuple[pd.DataFrame, dict]:
    region = args.region.lower()
    repeat = pd.read_csv(args.repeat, sep="\t")
    repeat = repeat[repeat["tissue"].eq(TISSUE[region]) & repeat["population"].eq("AA")].copy()
    repressive = pd.read_csv(args.repressive, sep="\t")
    repressive = repressive[
        repressive["tissue"].eq(TISSUE[region]) & repressive["population"].eq("AA")
    ].copy()
    for frame in [repeat, repressive]:
        frame["vmr_coord_id"] = _coord(frame, "seqnames", "start", "end")
    if repeat["vmr_coord_id"].duplicated().any() or repressive["vmr_coord_id"].duplicated().any():
        raise SystemExit(f"Duplicate annotation coordinates for {region}")
    keep_repeat = [
        "vmr_coord_id", "h2_unscaled", "num_snps", "vmr_length",
        "in_LINE", "in_L1", "in_any_repeat",
    ]
    keep_repressive = ["vmr_coord_id", "in_H3K9me3", "in_Quies"]
    d = repeat[keep_repeat].merge(repressive[keep_repressive], on="vmr_coord_id", validate="one_to_one")
    d["in_LINE_L1"] = (
        pd.to_numeric(d["in_LINE"], errors="coerce").fillna(0).eq(1)
        | pd.to_numeric(d["in_L1"], errors="coerce").fillna(0).eq(1)
    ).astype(int)
    for feature in ["in_H3K9me3", "in_Quies"]:
        d[feature] = pd.to_numeric(d[feature], errors="coerce").fillna(0).clip(0, 1).astype(int)
    d = d.rename(columns={
        "h2_unscaled": "legacy_h2", "vmr_length": "length",
    })

    calibrated = pd.read_csv(args.calibrated, sep="\t")
    calibrated = calibrated[calibrated["region"].astype(str).str.lower().eq(region)].copy()
    calibrated["vmr_coord_id"] = _coord(calibrated, "chromosome", "start", "end")
    calibrated = calibrated.drop_duplicates("vmr_coord_id", keep="first")
    calibrated["h2_upper_boundary_hit"] = calibrated["h2_upper_boundary_hit"].astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )
    calibrated = calibrated[calibrated["calibration_status"].eq("within_domain")]
    d = d.merge(
        calibrated[["vmr_coord_id", "h2_en_calibrated", "h2_upper_boundary_hit"]].rename(
            columns={"h2_en_calibrated": "calibrated_h2"}
        ), on="vmr_coord_id", how="inner", validate="one_to_one",
    )

    technical = pd.read_csv(args.technical, sep="\t")
    technical["vmr_coord_id"] = _coord(technical, "chrom", "start", "end")
    technical = technical.drop_duplicates("vmr_coord_id", keep="first")
    tech_columns = [column for column in [
        "umap_k24_mean", "high_mappability", "segdup_frac", "overlaps_segdup",
        "snp_prox_frac", "overlaps_snp_prox", "blacklist_frac",
    ] if column in technical.columns]
    d = d.merge(
        technical[["vmr_coord_id"] + tech_columns], on="vmr_coord_id", how="left", validate="one_to_one"
    )
    d["legacy_h2"] = pd.to_numeric(d["legacy_h2"], errors="coerce")
    d["calibrated_h2"] = pd.to_numeric(d["calibrated_h2"], errors="coerce")
    d["length"] = pd.to_numeric(d["length"], errors="coerce")
    d["num_snps"] = pd.to_numeric(d["num_snps"], errors="coerce")
    qc = {
        "region": region,
        "n_annotation_vmrs": int(len(repeat)),
        "n_exact_calibrated_within_domain": int(len(d)),
        "n_upper_boundary_hits": int(d["h2_upper_boundary_hit"].sum()),
        "n_with_umap": int(d.get("umap_k24_mean", pd.Series(dtype=float)).notna().sum()),
        "join_unique": True,
    }
    return d, qc


def fit_logistic(
    frame: pd.DataFrame,
    *,
    outcome: str,
    predictor: str,
    covariates: list[str],
    estimator: str,
    sensitivity: str,
) -> dict:
    columns = [outcome, predictor] + covariates
    use = frame.dropna(subset=columns).copy()
    if len(use) < 50 or use[outcome].nunique() < 2:
        return {
            "feature": outcome, "estimator": estimator, "sensitivity": sensitivity,
            "n_vmrs": int(len(use)), "error": "insufficient_rows_or_outcome_variation",
        }
    exog = pd.DataFrame({predictor: zscore(use[predictor])}, index=use.index)
    for column in covariates:
        exog[column] = zscore(use[column])
    exog = sm.add_constant(exog, has_constant="add")
    try:
        result = sm.GLM(use[outcome].astype(float), exog, family=sm.families.Binomial()).fit(
            cov_type="HC3", maxiter=250
        )
        estimate = float(result.params[predictor])
        se = float(result.bse[predictor])
        return {
            "feature": outcome,
            "estimator": estimator,
            "sensitivity": sensitivity,
            "term": predictor,
            "estimate_log_or_per_sd": estimate,
            "or_per_sd": float(np.exp(estimate)),
            "se": se,
            "ci_lower": float(estimate - 1.96 * se),
            "ci_upper": float(estimate + 1.96 * se),
            "pvalue": float(result.pvalues[predictor]),
            "n_vmrs": int(len(use)),
            "n_feature": int(use[outcome].sum()),
            "covariates": ",".join(covariates),
            "converged": bool(result.converged),
            "error": "",
        }
    except Exception as error:  # noqa: BLE001
        return {
            "feature": outcome, "estimator": estimator, "sensitivity": sensitivity,
            "n_vmrs": int(len(use)), "error": str(error),
        }


def model_sensitivities(d: pd.DataFrame) -> list[dict]:
    rows = []
    features = ["in_LINE_L1", "in_H3K9me3", "in_Quies"]
    estimators = {"calibrated": "calibrated_h2", "legacy": "legacy_h2"}
    candidate_adjustment = ["length", "num_snps", "umap_k24_mean", "segdup_frac", "snp_prox_frac"]
    adjustment = [
        column for column in candidate_adjustment
        if column in d.columns and d[column].notna().sum() >= max(50, int(0.5 * len(d)))
    ]
    subsets: dict[str, tuple[pd.DataFrame, list[str]]] = {
        "original": (d, []),
        "adjusted_prespecified": (d, adjustment),
        "exclude_upper_boundary": (d.loc[~d["h2_upper_boundary_hit"]], adjustment),
    }
    if "umap_k24_mean" in d.columns:
        subsets["high_mappability"] = (d.loc[d["umap_k24_mean"].ge(0.90)], [c for c in adjustment if c != "umap_k24_mean"])
    if "overlaps_snp_prox" in d.columns:
        subsets["snp_proximity_excluded"] = (d.loc[pd.to_numeric(d["overlaps_snp_prox"], errors="coerce").fillna(0).eq(0)], adjustment)
    if "overlaps_segdup" in d.columns:
        subsets["segdup_excluded"] = (d.loc[pd.to_numeric(d["overlaps_segdup"], errors="coerce").fillna(0).eq(0)], adjustment)

    for sensitivity, (subset, covariates) in subsets.items():
        for feature in features:
            # Both estimators use the same complete rows within every sensitivity.
            common = subset.dropna(subset=[feature, "calibrated_h2", "legacy_h2"] + covariates)
            for estimator, predictor in estimators.items():
                rows.append(fit_logistic(
                    common, outcome=feature, predictor=predictor, covariates=covariates,
                    estimator=estimator, sensitivity=sensitivity,
                ))
    return rows


def matched_sensitivities(
    d: pd.DataFrame, *, seed: int, n_permutations: int
) -> tuple[list[dict], pd.DataFrame]:
    match_columns = [
        column for column in ["length", "num_snps", "umap_k24_mean"]
        if column in d.columns and d[column].notna().sum() >= max(50, int(0.4 * len(d)))
    ]
    features = ["in_LINE_L1", "in_H3K9me3", "in_Quies"]
    rows = []
    balances = []
    for offset, (estimator, predictor) in enumerate(
        [("calibrated", "calibrated_h2"), ("legacy", "legacy_h2")]
    ):
        use = d.dropna(subset=[predictor] + match_columns + features).copy()
        pairs, balance, meta = greedy_nearest_neighbor_pairs(
            use,
            exposure=predictor,
            outcome=features[0],
            numeric_covariates=match_columns,
            high_quantile=0.8,
            low_quantile=0.2,
            caliper_sd=0.25,
            seed=seed + offset,
        )
        if not balance.empty:
            balance.insert(0, "estimator", estimator)
            balances.append(balance)
        if pairs.empty:
            for feature in features:
                rows.append({
                    "feature": feature, "estimator": estimator,
                    "sensitivity": "matched_top_vs_bottom_quintile", **meta,
                })
            continue
        high_index = pairs["high_index"].to_numpy()
        low_index = pairs["low_index"].to_numpy()
        for feature in features:
            differences = use.loc[high_index, feature].to_numpy(float) - use.loc[low_index, feature].to_numpy(float)
            rows.append({
                "feature": feature,
                "estimator": estimator,
                "sensitivity": "matched_top_vs_bottom_quintile",
                "n_pairs": int(len(pairs)),
                "mean_feature_high": float(use.loc[high_index, feature].mean()),
                "mean_feature_low": float(use.loc[low_index, feature].mean()),
                "mean_difference": float(differences.mean()),
                "permutation_pvalue": paired_randomization_pvalue(
                    differences, seed=seed + offset, n_perm=n_permutations
                ),
                "matching_variables": ",".join(match_columns),
                "max_abs_smd_after": float(meta["max_abs_smd_after"]),
                "balance_pass": bool(meta["max_abs_smd_after"] <= 0.10),
            })
    return rows, pd.concat(balances, ignore_index=True) if balances else pd.DataFrame()


def main() -> None:
    args = parse_args()
    region = args.region.lower()
    d, qc = load_analysis(args)
    model_rows = model_sensitivities(d)
    matched_rows, balance = matched_sensitivities(
        d, seed=args.seed, n_permutations=args.n_permutations
    )
    for rows in [model_rows, matched_rows]:
        for row in rows:
            row["region"] = region
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(model_rows).to_csv(output / "annotation_models.tsv", sep="\t", index=False)
    pd.DataFrame(matched_rows).to_csv(output / "annotation_matched.tsv", sep="\t", index=False)
    if not balance.empty:
        balance.insert(0, "region", region)
        balance.to_csv(output / "annotation_match_balance.tsv", sep="\t", index=False)
    write_tsv(output / "annotation_join_qc.tsv", [qc])
    print(f"Wrote calibrated annotation sensitivities under {output}")


if __name__ == "__main__":
    main()

