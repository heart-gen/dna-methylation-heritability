#!/usr/bin/env python3
"""Shared reviewer-defensible statistical helpers for meQTL validation."""

from __future__ import annotations

import numpy as np
import pandas as pd


def zscore(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    sd = x.std(ddof=0)
    return x * 0.0 if not np.isfinite(sd) or sd == 0 else (x - x.mean()) / sd


def standardized_mean_difference(high: pd.Series, low: pd.Series) -> float:
    a = pd.to_numeric(high, errors="coerce").dropna().to_numpy(dtype=float)
    b = pd.to_numeric(low, errors="coerce").dropna().to_numpy(dtype=float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
    if not np.isfinite(pooled) or pooled == 0:
        return 0.0
    return float((np.mean(a) - np.mean(b)) / pooled)


def greedy_nearest_neighbor_pairs(
    df: pd.DataFrame,
    *,
    exposure: str,
    outcome: str,
    numeric_covariates: list[str],
    exact_covariates: list[str] | None = None,
    high_quantile: float = 0.8,
    low_quantile: float = 0.2,
    caliper_sd: float = 0.25,
    seed: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Match exposure extremes without replacement and report covariate balance."""
    exact_covariates = exact_covariates or []
    required = [exposure, outcome] + numeric_covariates + exact_covariates
    d = df.dropna(subset=required).copy()
    if len(d) < 20:
        return pd.DataFrame(), pd.DataFrame(), {"error": "too_few_complete_rows"}
    q_hi = d[exposure].quantile(high_quantile)
    q_lo = d[exposure].quantile(low_quantile)
    hi = d[d[exposure] >= q_hi].copy()
    lo = d[d[exposure] <= q_lo].copy()
    if hi.empty or lo.empty:
        return pd.DataFrame(), pd.DataFrame(), {"error": "empty_exposure_extreme"}

    for col in numeric_covariates:
        center = d[col].astype(float).mean()
        scale = d[col].astype(float).std(ddof=0)
        if not np.isfinite(scale) or scale == 0:
            scale = 1.0
        hi[f"__z_{col}"] = (hi[col].astype(float) - center) / scale
        lo[f"__z_{col}"] = (lo[col].astype(float) - center) / scale

    # Estimate a propensity score for membership in the high-predictability
    # extreme. The caliper is expressed in SD units of the propensity logit.
    import statsmodels.api as sm

    combined = pd.concat([hi.assign(__treated=1), lo.assign(__treated=0)], axis=0)
    propensity_cols = [f"__z_{c}" for c in numeric_covariates]
    propensity_x = sm.add_constant(combined[propensity_cols], has_constant="add")
    try:
        propensity_fit = sm.GLM(
            combined["__treated"], propensity_x, family=sm.families.Binomial()
        ).fit()
        combined["__propensity_logit"] = propensity_fit.predict(
            propensity_x, which="linear"
        )
    except Exception:
        # Deterministic fallback: the first standardized covariate remains a
        # transparent scalar matching score.
        combined["__propensity_logit"] = (
            combined[propensity_cols[0]] if propensity_cols else 0.0
        )
    score_sd = float(combined["__propensity_logit"].std(ddof=0))
    if not np.isfinite(score_sd) or score_sd == 0:
        score_sd = 1.0
    hi["__propensity_logit"] = combined.loc[hi.index, "__propensity_logit"]
    lo["__propensity_logit"] = combined.loc[lo.index, "__propensity_logit"]

    rng = np.random.default_rng(seed)
    hi = hi.iloc[rng.permutation(len(hi))].copy()
    available = np.ones(len(lo), dtype=bool)
    pairs: list[dict] = []
    for hi_idx, hrow in hi.iterrows():
        candidate = np.flatnonzero(available)
        if exact_covariates:
            candidate = np.asarray([
                j for j in candidate
                if all(str(lo.iloc[j][c]) == str(hrow[c]) for c in exact_covariates)
            ], dtype=int)
        if not len(candidate):
            continue
        distance = np.abs(
            lo.iloc[candidate]["__propensity_logit"].to_numpy(float)
            - float(hrow["__propensity_logit"])
        ) / score_sd
        best_local = int(np.argmin(distance))
        if np.isfinite(caliper_sd) and distance[best_local] > caliper_sd:
            continue
        j = int(candidate[best_local])
        available[j] = False
        pairs.append({
            "high_index": hi_idx,
            "low_index": lo.index[j],
            "distance": float(distance[best_local]),
            "outcome_high": float(hrow[outcome]),
            "outcome_low": float(lo.iloc[j][outcome]),
        })
    pair_df = pd.DataFrame(pairs)
    if pair_df.empty:
        return pair_df, pd.DataFrame(), {"error": "no_pairs_within_caliper"}

    hi_post = d.loc[pair_df["high_index"]]
    lo_post = d.loc[pair_df["low_index"]]
    balance = []
    for col in numeric_covariates:
        balance.append({
            "covariate": col,
            "smd_before": standardized_mean_difference(hi[col], lo[col]),
            "smd_after": standardized_mean_difference(hi_post[col], lo_post[col]),
        })
    balance_df = pd.DataFrame(balance)
    meta = {
        "n_high_available": int(len(hi)),
        "n_low_available": int(len(lo)),
        "n_pairs": int(len(pair_df)),
        "caliper_sd": float(caliper_sd),
        "match_score": "propensity_logit",
        "max_abs_smd_after": (
            float(balance_df["smd_after"].abs().max()) if len(balance_df) else 0.0
        ),
    }
    return pair_df, balance_df, meta


def paired_randomization_pvalue(differences, *, seed: int, n_perm: int = 10000) -> float:
    """Two-sided matched-pair randomization test using within-pair label swaps."""
    diff = np.asarray(differences, dtype=float)
    diff = diff[np.isfinite(diff)]
    if not len(diff):
        return np.nan
    observed = abs(float(diff.mean()))
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(n_perm):
        signs = rng.choice((-1.0, 1.0), size=len(diff), replace=True)
        exceed += abs(float(np.mean(diff * signs))) >= observed
    return float((exceed + 1) / (n_perm + 1))


def greedy_binary_propensity_pairs(
    df: pd.DataFrame,
    *,
    exposure: str,
    outcome: str,
    numeric_covariates: list[str],
    exact_covariates: list[str] | None = None,
    caliper_sd: float = 0.25,
    seed: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Match binary-exposed rows 1:1 without replacement on propensity logit."""
    import statsmodels.api as sm

    exact_covariates = exact_covariates or []
    required = [exposure, outcome] + numeric_covariates + exact_covariates
    d = df.dropna(subset=required).copy()
    d[exposure] = d[exposure].astype(bool)
    treated = d[d[exposure]].copy()
    control = d[~d[exposure]].copy()
    if len(treated) < 10 or len(control) < 10:
        return pd.DataFrame(), pd.DataFrame(), {
            "error": "too_few_rows_in_binary_exposure_group",
            "n_exposed_available": int(len(treated)),
            "n_unexposed_available": int(len(control)),
        }

    standardized = []
    for column in numeric_covariates:
        center = pd.to_numeric(d[column], errors="coerce").mean()
        scale = pd.to_numeric(d[column], errors="coerce").std(ddof=0)
        if not np.isfinite(scale) or scale == 0:
            scale = 1.0
        name = f"__z_{column}"
        d[name] = (pd.to_numeric(d[column], errors="coerce") - center) / scale
        standardized.append(name)

    propensity_x = sm.add_constant(d[standardized], has_constant="add")
    try:
        propensity_fit = sm.GLM(
            d[exposure].astype(int), propensity_x, family=sm.families.Binomial()
        ).fit(maxiter=250)
        d["__propensity_logit"] = propensity_fit.predict(propensity_x, which="linear")
    except Exception:
        d["__propensity_logit"] = d[standardized[0]] if standardized else 0.0
    score_sd = float(d["__propensity_logit"].std(ddof=0))
    if not np.isfinite(score_sd) or score_sd == 0:
        score_sd = 1.0

    treated = d[d[exposure]].copy()
    control = d[~d[exposure]].copy()
    rng = np.random.default_rng(seed)
    treated = treated.iloc[rng.permutation(len(treated))]
    available = np.ones(len(control), dtype=bool)
    pairs: list[dict] = []
    for treated_idx, treated_row in treated.iterrows():
        candidate = np.flatnonzero(available)
        if exact_covariates:
            candidate = np.asarray([
                j for j in candidate
                if all(
                    str(control.iloc[j][column]) == str(treated_row[column])
                    for column in exact_covariates
                )
            ], dtype=int)
        if not len(candidate):
            continue
        distance = np.abs(
            control.iloc[candidate]["__propensity_logit"].to_numpy(dtype=float)
            - float(treated_row["__propensity_logit"])
        ) / score_sd
        best_local = int(np.argmin(distance))
        if np.isfinite(caliper_sd) and distance[best_local] > caliper_sd:
            continue
        j = int(candidate[best_local])
        available[j] = False
        pairs.append({
            "exposed_index": treated_idx,
            "unexposed_index": control.index[j],
            "distance": float(distance[best_local]),
            "outcome_exposed": float(treated_row[outcome]),
            "outcome_unexposed": float(control.iloc[j][outcome]),
        })
    pair_df = pd.DataFrame(pairs)
    if pair_df.empty:
        return pair_df, pd.DataFrame(), {
            "error": "no_pairs_within_caliper",
            "n_exposed_available": int(len(treated)),
            "n_unexposed_available": int(len(control)),
        }

    exposed_post = d.loc[pair_df["exposed_index"]]
    unexposed_post = d.loc[pair_df["unexposed_index"]]
    balance = pd.DataFrame([
        {
            "covariate": column,
            "smd_before": standardized_mean_difference(treated[column], control[column]),
            "smd_after": standardized_mean_difference(
                exposed_post[column], unexposed_post[column]
            ),
        }
        for column in numeric_covariates
    ])
    meta = {
        "n_exposed_available": int(len(treated)),
        "n_unexposed_available": int(len(control)),
        "n_pairs": int(len(pair_df)),
        "caliper_sd": float(caliper_sd),
        "match_score": "propensity_logit",
        "max_abs_smd_after": (
            float(balance["smd_after"].abs().max()) if len(balance) else 0.0
        ),
    }
    return pair_df, balance, meta
