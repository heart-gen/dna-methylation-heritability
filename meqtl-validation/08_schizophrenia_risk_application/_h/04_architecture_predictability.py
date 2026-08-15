#!/usr/bin/env python3
"""Analysis 4: do SCZ risk-linked meQTLs favor higher VMR predictability?

Primary: among tested risk-variant–CpG pairs (or VMR-level aggregation),
meQTL support ~ continuous local genetic predictability (+ technical covariates)
and matched permutations.
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
SEED = 20260801


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", default="caudate")
    p.add_argument("--outdir", default="")
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--n-perm", type=int, default=2000)
    return p.parse_args()


def _z(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x * 0.0
    return (x - x.mean()) / sd


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir) if args.outdir else (
        PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m" / args.region
    )
    tests_path = outdir / "risk_variant_cpg_meqtl.tsv.gz"
    if not tests_path.exists():
        raise SystemExit(f"Missing {tests_path}; run 03_test_risk_variant_cpg_meqtl.py first")
    res = pd.read_csv(tests_path, sep="\t")

    burden = pd.read_csv(
        PROJECT / "meqtl-validation/02_vmr_meqtl_burden/_m" / args.region / "vmr_meqtl_burden.tsv.gz",
        sep="\t",
    )
    burden["vmr_id"] = burden["vmr_id"].astype(str)
    if "analysis_schema_version" not in burden or not burden["analysis_schema_version"].eq(ANALYSIS_SCHEMA_VERSION).all():
        raise SystemExit("Stale Phase 2 burden table; regenerate repair schema v2 before Phase 7")

    # VMR-level: any significant risk-variant meQTL among tested CpGs in that VMR
    res["vmr_id"] = res["vmr_id"].astype(str)
    res["sig"] = res["significant_fdr"].astype(bool) if "significant_fdr" in res.columns else res["qval"].le(0.05)
    vmr = res.groupby("vmr_id", as_index=False).agg(
        n_tested_pairs=("phenotype_id", "size"),
        n_sig_pairs=("sig", "sum"),
        min_pval=("pval_nominal", "min"),
        local_predictability=("local_predictability", "first"),
        n_index_snps=("index_snp", "nunique"),
        n_loci=("locus_id", "nunique"),
    )
    vmr["any_sig"] = vmr["n_sig_pairs"] > 0
    vmr = vmr.merge(
        burden[[
            c for c in [
                "vmr_id", "n_tested_cpgs", "proportion_cpgs_with_sig_meqtl",
                "average_cpg_coverage", "mean_cpg_variance", "length",
                "vmr_mean_methylation", "cpg_density", "mean_num_tested_snps_per_cpg",
                "umap_k24_mean", "line_l1_frac", "segdup_frac",
            ] if c in burden.columns
        ]],
        on="vmr_id",
        how="left",
    )

    model_rows = []
    d = vmr.dropna(subset=["local_predictability", "any_sig"]).copy()
    d["any_sig"] = d["any_sig"].astype(int)
    specs = [
        ("unadjusted", ["local_predictability"]),
        ("adjusted_n_pairs", ["local_predictability", "n_tested_pairs"]),
    ]
    tech = [c for c in [
        "average_cpg_coverage", "mean_cpg_variance", "vmr_mean_methylation",
        "length", "cpg_density", "mean_num_tested_snps_per_cpg",
        "umap_k24_mean", "line_l1_frac", "segdup_frac",
    ]
            if c in d.columns and d[c].notna().sum() >= max(30, int(0.4 * len(d)))]
    if tech:
        specs.append(("adjusted_technical", ["local_predictability", "n_tested_pairs"] + tech))

    for name, cols in specs:
        use = d.dropna(subset=cols).copy()
        if use["any_sig"].nunique() < 2 or len(use) < 20:
            model_rows.append({"model": name, "error": "insufficient_variation_or_n", "n_vmrs": len(use)})
            continue
        try:
            exog = use[cols].apply(_z)
            exog = sm.add_constant(exog, has_constant="add")
            fit = sm.GLM(use["any_sig"], exog, family=sm.families.Binomial()).fit(cov_type="HC3")
            model_rows.append({
                "model": name,
                "n_vmrs": len(use),
                "n_vmrs_with_sig": int(use["any_sig"].sum()),
                "covariates": ",".join(cols),
                "coef_predictability": float(fit.params.get("local_predictability", np.nan)),
                "se_predictability": float(fit.bse.get("local_predictability", np.nan)),
                "pval_predictability": float(fit.pvalues.get("local_predictability", np.nan)),
                "converged": bool(fit.converged),
            })
        except Exception as exc:  # noqa: BLE001
            model_rows.append({"model": name, "error": str(exc)})

    # Pair-level logistic: significant_fdr ~ predictability
    pair_rows = []
    pr = res.dropna(subset=["local_predictability", "sig"]).copy()
    pr["sig"] = pr["sig"].astype(int)
    for name, cols in [
        ("pair_unadjusted", ["local_predictability"]),
        ("pair_adjusted_distance", ["local_predictability", "cpg_to_variant_distance", "maf"]),
    ]:
        use_cols = [c for c in cols if c in pr.columns]
        use = pr.dropna(subset=use_cols).copy()
        if use["sig"].nunique() < 2:
            pair_rows.append({"model": name, "error": "no_positive_or_negative_pairs"})
            continue
        try:
            exog = use[use_cols].apply(_z)
            exog = sm.add_constant(exog, has_constant="add")
            fit = sm.GLM(use["sig"], exog, family=sm.families.Binomial()).fit(
                cov_type="cluster", cov_kwds={"groups": use["vmr_id"]}
            )
            pair_rows.append({
                "model": name,
                "n_pairs": len(use),
                "n_sig_pairs": int(use["sig"].sum()),
                "covariates": ",".join(use_cols),
                "coef_predictability": float(fit.params.get("local_predictability", np.nan)),
                "se_predictability": float(fit.bse.get("local_predictability", np.nan)),
                "pval_predictability": float(fit.pvalues.get("local_predictability", np.nan)),
            })
        except Exception as exc:  # noqa: BLE001
            pair_rows.append({"model": name, "error": str(exc)})

    # Matched permutation: high vs low predictability among VMRs tested for SCZ
    perm_rows = []
    match_cols = ["n_tested_pairs"]
    for c in [
        "average_cpg_coverage", "mean_cpg_variance", "vmr_mean_methylation",
        "length", "cpg_density", "mean_num_tested_snps_per_cpg", "umap_k24_mean",
    ]:
        if c in d.columns and d[c].notna().sum() >= 20:
            match_cols.append(c)
    dd = d.dropna(subset=["local_predictability", "any_sig"] + match_cols).copy()
    balance = pd.DataFrame()
    if len(dd) >= 20 and dd["any_sig"].sum() >= 3:
        thresholds = load_yaml("analysis_thresholds.yml")["matching"]
        pairs, balance, meta = greedy_nearest_neighbor_pairs(
            dd,
            exposure="local_predictability",
            outcome="any_sig",
            numeric_covariates=match_cols,
            caliper_sd=float(thresholds["caliper_sd"]),
            seed=args.seed,
        )
        if not pairs.empty:
            differences = pairs["outcome_high"] - pairs["outcome_low"]
            obs = float(differences.mean())
            p_perm = paired_randomization_pvalue(
                differences, seed=args.seed, n_perm=args.n_perm
            )
            perm_rows.append({
                "analysis": "matched_high_vs_low_predictability",
                "n_pairs": int(len(pairs)),
                "mean_any_sig_high": float(pairs["outcome_high"].mean()),
                "mean_any_sig_low": float(pairs["outcome_low"].mean()),
                "mean_difference": obs,
                "permutation_pvalue": float(p_perm),
                "matching_variables": ",".join(match_cols),
                "max_abs_smd_after": meta["max_abs_smd_after"],
                "balance_pass": bool(meta["max_abs_smd_after"] <= float(thresholds["balance_smd_max"])),
                "permutation_scheme": "within_pair_label_swap",
            })
        else:
            perm_rows.append({"analysis": "matched_high_vs_low_predictability", **meta})
    else:
        perm_rows.append({"analysis": "matched_high_vs_low_predictability", "error": "too_few_vmrs_or_signals"})

    for r in model_rows:
        r["region"] = args.region
        r["unit"] = "vmr"
    for r in pair_rows:
        r["region"] = args.region
        r["unit"] = "pair"
    for r in perm_rows:
        r["region"] = args.region

    write_tsv(outdir / "architecture_model_results.tsv", model_rows + pair_rows)
    write_tsv(outdir / "architecture_matched_permutation.tsv", perm_rows)
    if not balance.empty:
        balance.to_csv(outdir / "architecture_match_balance.tsv", sep="\t", index=False)
    vmr.to_csv(outdir / "scz_vmr_level_summary.tsv.gz", sep="\t", index=False, compression="gzip")

    # Decision sketch for §12.12 criterion 1–2
    n_sig_loci = int(res.loc[res["sig"], "locus_id"].nunique()) if res["sig"].any() else 0
    tech = next((r for r in model_rows if r.get("model") == "adjusted_technical" and "coef_predictability" in r), None)
    unadj = next((r for r in model_rows if r.get("model") == "unadjusted" and "coef_predictability" in r), None)
    pick = tech or unadj
    write_tsv(outdir / "architecture_decision_snapshot.tsv", [{
        "region": args.region,
        "n_sig_risk_loci": n_sig_loci,
        "n_sig_pairs": int(res["sig"].sum()),
        "n_sig_vmrs": int(vmr["any_sig"].sum()),
        "primary_model": pick.get("model") if pick else "none",
        "coef_predictability": pick.get("coef_predictability", np.nan) if pick else np.nan,
        "pval_predictability": pick.get("pval_predictability", np.nan) if pick else np.nan,
        "matched_diff": perm_rows[0].get("mean_difference", np.nan),
        "matched_p": perm_rows[0].get("permutation_pvalue", np.nan),
        "criterion_ge1_locus_meqtl": n_sig_loci >= 1,
        "criterion_enrich_higher_predictability": bool(
            pick and pick.get("coef_predictability", 0) > 0 and pick.get("pval_predictability", 1) < 0.05
        ),
    }])
    print(f"Architecture results under {outdir}; sig loci={n_sig_loci}")


if __name__ == "__main__":
    main()
