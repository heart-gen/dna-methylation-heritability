#!/usr/bin/env python3
"""Phase 3: external brain meQTL support ~ continuous VMR predictability.

Primary model:
  VMR any-external-meQTL-support ~ local_predictability + technical/genomic covariates

Resources are analyzed separately (do not pool platforms). Preferred tissue
pairings are primary; cross-region tests are secondary.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import ANALYSIS_SCHEMA_VERSION, load_yaml, write_tsv  # noqa: E402
from _lib.stats_utils import greedy_nearest_neighbor_pairs, paired_randomization_pvalue  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
PHASE3 = PROJECT / "meqtl-validation" / "03_external_meqtl_validation" / "_m"
PHASE2 = PROJECT / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m"

# Preferred tissue pairing for primary interpretation
PREFERRED_REGION = {
    "jaffe_dlpfc_450k_meqtl": "dlpfc",
    "schulz_hippocampus_array_meqtl": "hippocampus",
    "brainseq_wgbs_meqtl_scz_subset": "dlpfc",  # exploratory SCZ-risk subset only
}

EXPLORATORY_RESOURCES = {"brainseq_wgbs_meqtl_scz_subset"}
# Resources with no recoverable tested-negative universe. Jaffe and Schulz are both
# Illumina 450K studies whose universe is supplied by the array manifest in
# `03_harmonize_external_meqtls.py`; BrainSeq is a WGBS SCZ-risk supplement with no
# genome-wide catalog (and is same-cohort, so it is exploratory regardless).
POSITIVE_ONLY_RESOURCES: set[str] = set()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True)
    p.add_argument("--resource-id", required=True)
    p.add_argument("--harmonized-external", default="")
    p.add_argument("--vmr-burden", default="")
    p.add_argument("--seed", type=int, default=20260730)
    return p.parse_args()


def _z(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x * 0.0
    return (x - x.mean()) / sd


def aggregate_vmr_external(ext: pd.DataFrame) -> pd.DataFrame:
    if "vmr_id" not in ext.columns:
        raise SystemExit("harmonized external table must include vmr_id")
    if "external_meqtl_support" not in ext.columns:
        raise SystemExit("harmonized external table must include external_meqtl_support")
    if "external_assayed" not in ext.columns:
        raise SystemExit("harmonized external table must include external_assayed")
    e = ext.copy()
    e = e[pd.to_numeric(e["external_assayed"], errors="coerce").eq(1)].copy()
    e["vmr_id"] = e["vmr_id"].astype(str)
    e["external_meqtl_support"] = (
        pd.to_numeric(e["external_meqtl_support"], errors="coerce").fillna(0).astype(int).clip(0, 1)
    )
    g = e.groupby("vmr_id", sort=False).agg(
        n_cpgs_annotated=("external_meqtl_support", "size"),
        n_cpgs_with_external_support=("external_meqtl_support", "sum"),
        external_meqtl_support=("external_meqtl_support", "max"),
        assay_universe_complete=("assay_universe_complete", "all"),
    ).reset_index()
    g["proportion_cpgs_with_external_support"] = (
        g["n_cpgs_with_external_support"] / g["n_cpgs_annotated"].replace(0, np.nan)
    )
    return g


def fit_support_models(d: pd.DataFrame) -> list[dict]:
    rows = []
    need = ["external_meqtl_support", "local_predictability"]
    base = d.dropna(subset=need).copy()
    if base.empty:
        return [{"model": "none", "error": "no complete cases"}]
    if not base["assay_universe_complete"].fillna(False).astype(bool).all():
        return [{
            "model": "not_estimable",
            "n_vmrs": len(base),
            "error": "external tested-negative assay universe unavailable; positive-control overlap only",
        }]

    specs: list[tuple[str, list[str]]] = [
        ("unadjusted", ["local_predictability"]),
    ]
    specs.append(("adjusted_minimal", ["local_predictability", "n_cpgs_annotated"]))

    tech = []
    for c in [
        "average_cpg_coverage", "mean_cpg_variance", "length",
        "umap_k24_mean", "line_l1_frac", "segdup_frac",
    ]:
        if c in base.columns and base[c].notna().sum() >= max(50, int(0.4 * len(base))):
            tech.append(c)
    if tech:
        specs.append(("adjusted_technical", ["local_predictability", "n_cpgs_annotated"] + tech))

    y_all = base["external_meqtl_support"].astype(float)
    for name, cols in specs:
        try:
            use = base.dropna(subset=cols).copy()
            if use.empty or use["external_meqtl_support"].nunique() < 2:
                rows.append({"model": name, "error": "insufficient outcome variation or cases"})
                continue
            y = use["external_meqtl_support"].astype(float)
            exog = use[cols].apply(_z)
            exog = sm.add_constant(exog, has_constant="add")
            res = sm.GLM(y, exog, family=sm.families.Binomial()).fit(cov_type="HC3")
            rows.append({
                "model": name,
                "n_vmrs": len(use),
                "n_supported_vmrs": int(y.sum()),
                "mean_external_support": float(y.mean()),
                "covariates": ",".join(cols),
                "coef_predictability": float(res.params.get("local_predictability", np.nan)),
                "se_predictability": float(res.bse.get("local_predictability", np.nan)),
                "pval_predictability": float(res.pvalues.get("local_predictability", np.nan)),
                "converged": bool(res.converged),
            })
        except Exception as exc:  # noqa: BLE001
            rows.append({"model": name, "error": str(exc)})

    # Secondary: binomial counts of external-supported CpGs within VMR
    if {"n_cpgs_with_external_support", "n_cpgs_annotated"}.issubset(base.columns):
        try:
            use = base.dropna(subset=["local_predictability", "n_cpgs_annotated"]).copy()
            use = use[use["n_cpgs_annotated"] > 0]
            use["n_ns"] = use["n_cpgs_annotated"] - use["n_cpgs_with_external_support"]
            endog = use[["n_cpgs_with_external_support", "n_ns"]]
            cols = ["local_predictability"]
            cols.append("n_cpgs_annotated")
            exog = sm.add_constant(use[cols].apply(_z), has_constant="add")
            model = sm.GLM(endog, exog, family=sm.families.Binomial())
            pilot = model.fit(maxiter=250)
            dispersion = max(1.0, float(pilot.pearson_chi2 / max(pilot.df_resid, 1)))
            res = model.fit(scale=dispersion, cov_type="HC3", maxiter=250)
            rows.append({
                "model": "binomial_cpg_proportion_secondary",
                "n_vmrs": len(use),
                "n_supported_vmrs": int((use["n_cpgs_with_external_support"] > 0).sum()),
                "mean_external_support": float((use["n_cpgs_with_external_support"] > 0).mean()),
                "covariates": ",".join(cols),
                "coef_predictability": float(res.params.get("local_predictability", np.nan)),
                "se_predictability": float(res.bse.get("local_predictability", np.nan)),
                "pval_predictability": float(res.pvalues.get("local_predictability", np.nan)),
                "converged": bool(res.converged),
            })
        except Exception as exc:  # noqa: BLE001
            rows.append({"model": "binomial_cpg_proportion_secondary", "error": str(exc)})

    rows.append({
        "model": "descriptive",
        "n_vmrs": len(base),
        "n_supported_vmrs": int(y_all.sum()),
        "mean_external_support": float(y_all.mean()),
        "covariates": "",
        "coef_predictability": np.nan,
        "se_predictability": np.nan,
        "pval_predictability": np.nan,
        "converged": True,
    })
    return rows


def matched_contrast(df: pd.DataFrame, seed: int) -> tuple[list[dict], pd.DataFrame]:
    if not df["assay_universe_complete"].fillna(False).astype(bool).all():
        return ([{
            "analysis": "matched",
            "error": "external tested-negative assay universe unavailable; positive-control overlap only",
        }], pd.DataFrame())
    match_cols = []
    for c in ["n_cpgs_annotated", "average_cpg_coverage", "mean_cpg_variance", "length", "umap_k24_mean"]:
        if c in df.columns and df[c].notna().sum() >= 50:
            match_cols.append(c)
    if not match_cols:
        match_cols = ["n_cpgs_annotated"] if "n_cpgs_annotated" in df.columns else []
    need = ["local_predictability", "external_meqtl_support"] + match_cols
    d = df.dropna(subset=need).copy()
    if len(d) < 40 or d["external_meqtl_support"].nunique() < 2:
        return ([{"analysis": "matched", "error": "too few VMRs or no outcome variation"}], pd.DataFrame())
    thresholds = load_yaml("analysis_thresholds.yml")["matching"]
    pairs, balance, meta = greedy_nearest_neighbor_pairs(
        d,
        exposure="local_predictability",
        outcome="external_meqtl_support",
        numeric_covariates=match_cols,
        caliper_sd=float(thresholds["caliper_sd"]),
        seed=seed,
    )
    if pairs.empty:
        return ([{"analysis": "matched", **meta}], balance)
    differences = pairs["outcome_high"] - pairs["outcome_low"]
    p_perm = paired_randomization_pvalue(differences, seed=seed, n_perm=10000)
    return ([{
        "analysis": "matched_high_vs_low_predictability",
        "n_pairs": int(len(pairs)),
        "mean_support_high": float(pairs["outcome_high"].mean()),
        "mean_support_low": float(pairs["outcome_low"].mean()),
        "mean_difference": float(differences.mean()),
        "permutation_pvalue": float(p_perm),
        "matching_variables": ",".join(match_cols),
        "max_abs_smd_after": meta["max_abs_smd_after"],
        "balance_pass": bool(meta["max_abs_smd_after"] <= float(thresholds["balance_smd_max"])),
        "permutation_scheme": "within_pair_label_swap",
    }], balance)


def main() -> None:
    args = parse_args()
    resource = args.resource_id
    region = args.region
    preferred = PREFERRED_REGION.get(resource)
    analysis_role = (
        "exploratory"
        if resource in EXPLORATORY_RESOURCES
        else ("supportive_positive_only" if resource in POSITIVE_ONLY_RESOURCES else
              ("primary" if preferred == region else "secondary_cross_region"))
    )

    harm = Path(args.harmonized_external) if args.harmonized_external else (
        PHASE3 / "harmonized" / f"{resource}.{region}.vmr_support.tsv.gz"
    )
    burden_path = Path(args.vmr_burden) if args.vmr_burden else (
        PHASE2 / region / "vmr_meqtl_burden.tsv.gz"
    )
    if not harm.exists():
        raise SystemExit(f"Missing harmonized external file: {harm}")
    if not burden_path.exists():
        raise SystemExit(f"Missing VMR burden table: {burden_path}")

    ext = pd.read_csv(harm, sep="\t", compression="infer")
    vmr_ext = aggregate_vmr_external(ext)
    burden = pd.read_csv(burden_path, sep="\t", compression="infer")
    if "analysis_schema_version" not in burden or not burden["analysis_schema_version"].eq(ANALYSIS_SCHEMA_VERSION).all():
        raise SystemExit("Stale Phase 2 burden table; regenerate repair schema v2 before external validation")
    burden["vmr_id"] = burden["vmr_id"].astype(str)
    d = burden.merge(vmr_ext, on="vmr_id", how="inner")
    d = d.dropna(subset=["local_predictability"])
    if d.empty:
        raise SystemExit("No overlapping VMRs with predictability and external annotations")

    outdir = PHASE3 / region
    outdir.mkdir(parents=True, exist_ok=True)

    # Persist VMR-level external annotation joined to burden keys
    keep = [
        "vmr_id", "region", "local_predictability", "n_tested_cpgs",
        "n_cpgs_annotated", "n_cpgs_with_external_support",
        "proportion_cpgs_with_external_support", "external_meqtl_support",
        "assay_universe_complete",
        "average_cpg_coverage", "mean_cpg_variance", "length", "umap_k24_mean",
        "line_l1_frac", "segdup_frac",
    ]
    keep = [c for c in keep if c in d.columns]
    d[keep].to_csv(
        outdir / f"vmr_external_support_{resource}.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )

    model_rows = fit_support_models(d)
    for r in model_rows:
        r["region"] = region
        r["resource_id"] = resource
        r["analysis_role"] = analysis_role
        r["preferred_region"] = preferred or ""
        r["assay_universe_complete"] = bool(
            d["assay_universe_complete"].fillna(False).astype(bool).all()
        )
    write_tsv(outdir / f"external_support_model_{resource}.tsv", model_rows)

    match_rows, balance = matched_contrast(d, args.seed)
    for r in match_rows:
        r["region"] = region
        r["resource_id"] = resource
        r["analysis_role"] = analysis_role
    write_tsv(outdir / f"external_matched_{resource}.tsv", match_rows)
    if not balance.empty:
        balance.insert(0, "resource_id", resource)
        balance.insert(0, "region", region)
        balance.to_csv(outdir / f"external_match_balance_{resource}.tsv", sep="\t", index=False)

    print(
        f"Wrote Phase 3 results for {resource} / {region} "
        f"(role={analysis_role}, n_vmrs={len(d)}, "
        f"support_rate={d['external_meqtl_support'].mean():.3f})"
    )


if __name__ == "__main__":
    main()
