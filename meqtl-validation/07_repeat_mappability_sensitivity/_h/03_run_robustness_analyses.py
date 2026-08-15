#!/usr/bin/env python3
"""
Phase 6: consolidated repeat / mappability / technical robustness analyses.

For each brain region, test whether high local genetic predictability remains
enriched for LINE/L1, H3K9me3, and quiescent chromatin, and whether the
predictability→meQTL-burden association survives technical restrictions.

Writes one consolidated table plus per-region detail TSVs.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import ANALYSIS_SCHEMA_VERSION, load_yaml, write_tsv  # noqa: E402
from _lib.stats_utils import greedy_nearest_neighbor_pairs, paired_randomization_pvalue  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
PHASE6 = PROJECT / "meqtl-validation" / "07_repeat_mappability_sensitivity" / "_m"
PHASE2 = PROJECT / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m"
ANN = PROJECT / "heritability/elastic_net_model/all_individuals/tissue_comparison/annotation"
SUPPORT = PROJECT / "inputs/supportfiles/_m"
SNP_WIN = SUPPORT / "common_snp_windows_pm150bp.hg38.bed.gz"

TISSUE_MAP = {
    "caudate": "Caudate",
    "dlpfc": "DLPFC",
    "hippocampus": "Hippocampus",
}
SEED = 20260730


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions", nargs="+", default=["caudate", "dlpfc", "hippocampus"])
    p.add_argument("--seed", type=int, default=SEED)
    return p.parse_args()


def _z(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x * 0.0
    return (x - x.mean()) / sd


def coord_id(chrom, start, end) -> str:
    c = str(chrom).replace("chr", "")
    return f"{c}:{int(start)}-{int(end)}"


def bedtools_frac(vmrs: pd.DataFrame, bed_b: Path) -> pd.Series:
    """Coverage fraction of VMR intervals overlapping bed_b."""
    if not bed_b.exists():
        return pd.Series(np.nan, index=vmrs.index)
    with tempfile.TemporaryDirectory(prefix="p6_cov_") as tmp:
        tmp = Path(tmp)
        a = tmp / "a.bed"
        out = tmp / "out.tsv"
        vmrs[["chrom", "start", "end"]].assign(name=[f"v{i}" for i in range(len(vmrs))]).to_csv(
            a, sep="\t", header=False, index=False
        )
        # bedtools coverage -a a -b b
        cmd = ["bedtools", "coverage", "-a", str(a), "-b", str(bed_b)]
        with out.open("w") as handle:
            subprocess.check_call(cmd, stdout=handle)
        cov = pd.read_csv(
            out, sep="\t", header=None,
            names=["chrom", "start", "end", "name", "n_overlaps", "n_bases", "length", "frac"],
        )
        return pd.Series(cov["frac"].to_numpy(dtype=float), index=vmrs.index)


def load_tech(region: str) -> pd.DataFrame:
    path = PHASE6 / region / "vmr_technical_annotations.tsv"
    t = pd.read_csv(path, sep="\t")
    t["coord_id"] = [
        coord_id(c, s, e) for c, s, e in zip(t["chrom"], t["start"], t["end"])
    ]
    if "snp_prox_frac" not in t.columns and SNP_WIN.exists():
        t["snp_prox_frac"] = bedtools_frac(t, SNP_WIN)
        t["overlaps_snp_prox"] = (t["snp_prox_frac"].fillna(0) > 0).astype(int)
        # persist for reuse
        t.to_csv(path, sep="\t", index=False)
    elif "snp_prox_frac" not in t.columns:
        t["snp_prox_frac"] = np.nan
        t["overlaps_snp_prox"] = 0
    return t


def load_feature_vmrs(region: str) -> pd.DataFrame:
    tissue = TISSUE_MAP[region]
    rep = pd.read_csv(ANN / "repeat_elements/_m/vmr_repeat_overlap_AA.tsv", sep="\t")
    rep = rep[rep["tissue"] == tissue].copy()
    reprs = pd.read_csv(ANN / "repressive_chromatin/_m/vmr_repressive_overlap_AA.tsv", sep="\t")
    reprs = reprs[reprs["tissue"] == tissue].copy()
    key_cols = ["seqnames", "start", "end"]
    keep_rep = key_cols + ["h2_unscaled", "vmr_length", "num_snps", "in_LINE", "in_L1", "in_any_repeat"]
    keep_repr = key_cols + ["in_H3K9me3", "in_Quies", "in_Het"]
    d = rep[keep_rep].merge(reprs[keep_repr], on=key_cols, how="inner")
    d["coord_id"] = [coord_id(c, s, e) for c, s, e in zip(d["seqnames"], d["start"], d["end"])]
    d["local_predictability"] = pd.to_numeric(d["h2_unscaled"], errors="coerce")
    d["length"] = pd.to_numeric(d["vmr_length"], errors="coerce")
    d["num_snps"] = pd.to_numeric(d["num_snps"], errors="coerce")
    for c in ["in_LINE", "in_L1", "in_any_repeat", "in_H3K9me3", "in_Quies", "in_Het"]:
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0).astype(int).clip(0, 1)
    # LINE/L1 combined flag preferred for manuscript claim
    d["in_LINE_L1"] = ((d["in_LINE"] == 1) | (d["in_L1"] == 1)).astype(int)
    return d


def join_tech(features: pd.DataFrame, tech: pd.DataFrame) -> pd.DataFrame:
    tcols = [
        "coord_id", "umap_k24_mean", "high_mappability", "segdup_frac",
        "overlaps_segdup", "blacklist_frac", "line_l1_frac",
        "snp_prox_frac", "overlaps_snp_prox",
    ]
    tcols = [c for c in tcols if c in tech.columns]
    return features.merge(tech[tcols], on="coord_id", how="left")


def fit_logistic(df: pd.DataFrame, ycol: str, xcols: list[str]) -> dict:
    use = df.dropna(subset=[ycol, "local_predictability"] + xcols).copy()
    if use.empty or use[ycol].nunique() < 2:
        return {"estimate": np.nan, "pvalue": np.nan, "n": len(use), "error": "insufficient"}
    y = use[ycol].astype(float)
    X = use[["local_predictability"] + xcols].apply(_z)
    X = sm.add_constant(X, has_constant="add")
    try:
        res = sm.GLM(y, X, family=sm.families.Binomial()).fit(cov_type="HC3")
        return {
            "estimate": float(res.params["local_predictability"]),
            "or": float(np.exp(res.params["local_predictability"])),
            "pvalue": float(res.pvalues["local_predictability"]),
            "n": int(len(use)),
            "n_feature": int(y.sum()),
            "error": "",
        }
    except Exception as exc:  # noqa: BLE001
        return {"estimate": np.nan, "pvalue": np.nan, "n": len(use), "error": str(exc)}


def matched_feature_delta(df: pd.DataFrame, ycol: str, match_cols: list[str], seed: int) -> dict:
    need = ["local_predictability", ycol] + match_cols
    d = df.dropna(subset=need).copy()
    if len(d) < 40 or d[ycol].nunique() < 2:
        return {"estimate": np.nan, "pvalue": np.nan, "n": len(d), "error": "too few"}
    thresholds = load_yaml("analysis_thresholds.yml")["matching"]
    pairs, _balance, meta = greedy_nearest_neighbor_pairs(
        d,
        exposure="local_predictability",
        outcome=ycol,
        numeric_covariates=match_cols,
        caliper_sd=float(thresholds["caliper_sd"]),
        seed=seed,
    )
    if pairs.empty:
        return {"estimate": np.nan, "pvalue": np.nan, "n": 0, "error": "no pairs"}
    hi_y = pairs["outcome_high"].to_numpy()
    lo_y_p = pairs["outcome_low"].to_numpy()
    differences = hi_y - lo_y_p
    diff = differences.mean()
    p = paired_randomization_pvalue(differences, seed=seed, n_perm=10000)
    return {
        "estimate": float(diff),
        "pvalue": float(p),
        "n": int(len(pairs)),
        "mean_high": float(hi_y.mean()),
        "mean_low": float(lo_y_p.mean()),
        "max_abs_smd_after": meta["max_abs_smd_after"],
        "error": "",
    }


def fit_burden(df: pd.DataFrame, xcols: list[str]) -> dict:
    need = ["local_predictability", "n_cpgs_with_sig_meqtl", "n_tested_cpgs"] + xcols
    d = df.dropna(subset=[c for c in need if c in df.columns]).copy()
    d = d[d["n_tested_cpgs"] > 0]
    if d.empty:
        return {"estimate": np.nan, "pvalue": np.nan, "n": 0, "error": "empty"}
    d["n_ns"] = d["n_tested_cpgs"] - d["n_cpgs_with_sig_meqtl"]
    cols = ["local_predictability"] + [c for c in xcols if c in d.columns]
    use = d.dropna(subset=cols)
    endog = use[["n_cpgs_with_sig_meqtl", "n_ns"]]
    exog = sm.add_constant(use[cols].apply(_z), has_constant="add")
    try:
        model = sm.GLM(endog, exog, family=sm.families.Binomial())
        pilot = model.fit(maxiter=250)
        dispersion = max(1.0, float(pilot.pearson_chi2 / max(pilot.df_resid, 1)))
        res = model.fit(scale=dispersion, cov_type="HC3", maxiter=250)
        return {
            "estimate": float(res.params["local_predictability"]),
            "pvalue": float(res.pvalues["local_predictability"]),
            "n": int(len(use)),
            "error": "",
        }
    except Exception as exc:  # noqa: BLE001
        return {"estimate": np.nan, "pvalue": np.nan, "n": len(use), "error": str(exc)}


def subset_filters(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    out = {"original": df}
    # adjusted uses same rows; filter subsets:
    if "high_mappability" in df.columns:
        out["high_mappability"] = df[df["high_mappability"].fillna(0).astype(int) == 1]
    if "overlaps_segdup" in df.columns:
        out["segdup_excluded"] = df[df["overlaps_segdup"].fillna(0).astype(int) == 0]
    if "overlaps_snp_prox" in df.columns:
        out["snp_proximity_excluded"] = df[df["overlaps_snp_prox"].fillna(0).astype(int) == 0]
    # combined technical clean set
    mask = pd.Series(True, index=df.index)
    if "high_mappability" in df.columns:
        mask &= df["high_mappability"].fillna(0).astype(int).eq(1)
    if "overlaps_segdup" in df.columns:
        mask &= df["overlaps_segdup"].fillna(0).astype(int).eq(0)
    if "overlaps_snp_prox" in df.columns:
        mask &= df["overlaps_snp_prox"].fillna(0).astype(int).eq(0)
    out["tech_clean"] = df[mask]
    return out


def analyze_feature(region: str, df: pd.DataFrame, feature: str, seed: int) -> dict:
    adj_cols = [c for c in ["length", "num_snps"] if c in df.columns and df[c].notna().sum() > 50]
    match_cols = [c for c in ["length", "num_snps", "umap_k24_mean"] if c in df.columns and df[c].notna().sum() > 50]
    if not match_cols:
        match_cols = ["length"] if "length" in df.columns else []

    subsets = subset_filters(df)
    row = {
        "region": region,
        "analysis": feature,
        "status": "ok",
    }
    # original
    o = fit_logistic(subsets["original"], feature, [])
    row["original_estimate"] = o.get("or", o["estimate"])
    row["original_pvalue"] = o["pvalue"]
    row["original_n"] = o["n"]
    # adjusted
    a = fit_logistic(subsets["original"], feature, adj_cols)
    row["adjusted_estimate"] = a.get("or", a["estimate"])
    row["adjusted_pvalue"] = a["pvalue"]
    # matched (rate difference hi-lo)
    m = matched_feature_delta(subsets["original"], feature, match_cols, seed)
    row["matched_estimate"] = m["estimate"]
    row["matched_pvalue"] = m["pvalue"]
    row["matched_n_pairs"] = m.get("n", np.nan)

    for key, col_est, col_p in [
        ("high_mappability", "high_mappability_estimate", "high_mappability_pvalue"),
        ("snp_proximity_excluded", "snp_proximity_excluded_estimate", "snp_proximity_excluded_pvalue"),
        ("segdup_excluded", "segdup_excluded_estimate", "segdup_excluded_pvalue"),
    ]:
        if key in subsets and len(subsets[key]) >= 50:
            r = fit_logistic(subsets[key], feature, adj_cols)
            row[col_est] = r.get("or", r["estimate"])
            row[col_p] = r["pvalue"]
            row[f"{key}_n"] = r["n"]
        else:
            row[col_est] = np.nan
            row[col_p] = np.nan

    # consistency vs original OR direction (OR>1 means predictability associated with feature)
    ests = [row.get("original_estimate"), row.get("adjusted_estimate"),
            row.get("high_mappability_estimate"), row.get("snp_proximity_excluded_estimate"),
            row.get("segdup_excluded_estimate")]
    # for matched, positive delta means higher feature rate in high-pred
    signs = []
    for e in ests:
        if pd.notna(e):
            signs.append(np.sign(np.log(e)) if e > 0 else 0)
    if pd.notna(row.get("matched_estimate")):
        signs.append(np.sign(row["matched_estimate"]))
    row["direction_consistent"] = bool(len(signs) >= 2 and len(set(signs) - {0}) <= 1)
    pvals = [row.get(k) for k in [
        "original_pvalue", "adjusted_pvalue", "matched_pvalue",
        "high_mappability_pvalue", "snp_proximity_excluded_pvalue", "segdup_excluded_pvalue",
    ] if pd.notna(row.get(k))]
    row["significance_consistent"] = bool(pvals and sum(p < 0.05 for p in pvals) >= max(2, len(pvals) // 2))
    return row


def analyze_burden(region: str, seed: int) -> dict:
    path = PHASE2 / region / "vmr_meqtl_burden.tsv.gz"
    if not path.exists():
        return {"region": region, "analysis": "predictability_meqtl_burden_association",
                "status": "missing_phase2_burden"}
    d = pd.read_csv(path, sep="\t")
    if "analysis_schema_version" not in d or not d["analysis_schema_version"].eq(ANALYSIS_SCHEMA_VERSION).all():
        raise SystemExit(f"{region}: stale Phase 2 burden; regenerate repair schema v2 before Phase 6")
    if "tech_join_source" not in d.columns:
        raise SystemExit(
            f"{region}: repair-v2 technical join missing; run step_2_tech_joins.sh "
            "between Phase 2 aggregation and modeling"
        )
    d = d.dropna(subset=["local_predictability", "n_tested_cpgs", "n_cpgs_with_sig_meqtl"])
    tech = load_tech(region)
    # join tech on vmr_id (task) or coord
    t = tech.copy()
    # Prefer task_id join, then interval_id for remaining gaps.
    feat_cols = [
        c for c in [
            "umap_k24_mean", "high_mappability", "overlaps_segdup",
            "overlaps_snp_prox", "snp_prox_frac", "length", "segdup_frac",
        ] if c in t.columns
    ]
    d["vmr_id"] = d["vmr_id"].astype(str)
    # Repair v2 completes technical joins before fitting Phase 2. Only backfill
    # fields that are absent from the burden table here.
    feat_cols = [c for c in feat_cols if c not in d.columns]
    if "task_id" in t.columns and feat_cols:
        by_task = t.dropna(subset=["task_id"]).copy()
        by_task["vmr_id"] = by_task["task_id"].astype(float).astype(int).astype(str)
        by_task = by_task.drop_duplicates("vmr_id", keep="first")
        merged = d.merge(by_task[["vmr_id"] + feat_cols], on="vmr_id", how="left")
    else:
        merged = d.copy()
        for c in feat_cols:
            if c not in merged.columns:
                merged[c] = np.nan

    if feat_cols and "interval_id" in t.columns:
        for c in feat_cols:
            if c not in merged.columns:
                merged[c] = np.nan
        miss = merged[feat_cols[0]].isna()
        if miss.any():
            by_int = t[["interval_id"] + feat_cols].copy()
            by_int["vmr_id"] = by_int["interval_id"].astype(str)
            by_int = by_int.drop_duplicates("vmr_id", keep="first")
            extra = merged.loc[miss, ["vmr_id"]].merge(
                by_int[["vmr_id"] + feat_cols], on="vmr_id", how="left"
            )
            for c in feat_cols:
                merged.loc[miss, c] = extra[c].to_numpy()

    if "overlaps_snp_prox" not in merged.columns and "snp_prox_frac" in merged.columns:
        merged["overlaps_snp_prox"] = (merged["snp_prox_frac"].fillna(0) > 0).astype(int)

    adj = [c for c in ["n_tested_cpgs", "average_cpg_coverage", "mean_cpg_variance"] if c in merged.columns]
    subsets = subset_filters(merged)
    row = {"region": region, "analysis": "predictability_meqtl_burden_association", "status": "ok"}
    o = fit_burden(subsets["original"], ["n_tested_cpgs"] if "n_tested_cpgs" in merged.columns else [])
    row["original_estimate"] = o["estimate"]
    row["original_pvalue"] = o["pvalue"]
    row["original_n"] = o["n"]
    a = fit_burden(subsets["original"], adj)
    row["adjusted_estimate"] = a["estimate"]
    row["adjusted_pvalue"] = a["pvalue"]
    # matched: reuse Phase 2 style on proportion
    mcols = [c for c in ["n_tested_cpgs", "average_cpg_coverage", "mean_cpg_variance", "length", "umap_k24_mean"]
             if c in merged.columns and merged[c].notna().sum() > 50]
    m = matched_feature_delta(
        subsets["original"].assign(
            _y=subsets["original"]["proportion_cpgs_with_sig_meqtl"]
        ),
        "_y", mcols, seed,
    ) if "proportion_cpgs_with_sig_meqtl" in merged.columns else {"estimate": np.nan, "pvalue": np.nan}
    row["matched_estimate"] = m.get("estimate", np.nan)
    row["matched_pvalue"] = m.get("pvalue", np.nan)

    for key, col_est, col_p in [
        ("high_mappability", "high_mappability_estimate", "high_mappability_pvalue"),
        ("snp_proximity_excluded", "snp_proximity_excluded_estimate", "snp_proximity_excluded_pvalue"),
        ("segdup_excluded", "segdup_excluded_estimate", "segdup_excluded_pvalue"),
    ]:
        if key in subsets and len(subsets[key]) >= 50:
            r = fit_burden(subsets[key], adj)
            row[col_est] = r["estimate"]
            row[col_p] = r["pvalue"]
        else:
            row[col_est] = np.nan
            row[col_p] = np.nan

    ests = [row.get(k) for k in [
        "original_estimate", "adjusted_estimate", "matched_estimate",
        "high_mappability_estimate", "snp_proximity_excluded_estimate", "segdup_excluded_estimate",
    ] if pd.notna(row.get(k))]
    row["direction_consistent"] = bool(ests and all(e > 0 for e in ests))
    pvals = [row.get(k) for k in [
        "original_pvalue", "adjusted_pvalue", "matched_pvalue",
        "high_mappability_pvalue", "snp_proximity_excluded_pvalue", "segdup_excluded_pvalue",
    ] if pd.notna(row.get(k))]
    row["significance_consistent"] = bool(pvals and sum(p < 0.05 for p in pvals) >= max(2, len(pvals) // 2))
    return row


def main() -> None:
    args = parse_args()
    thresholds = load_yaml("analysis_thresholds.yml")
    map_min = float(thresholds.get("repeat_sensitivity", {}).get("high_mappability_min", 0.9))
    all_rows = []

    for region in args.regions:
        print(f"==== {region} ====")
        tech = load_tech(region)
        # ensure high_mappability uses config threshold
        if "umap_k24_mean" in tech.columns:
            tech["high_mappability"] = (tech["umap_k24_mean"] >= map_min).astype(int)
        feat = join_tech(load_feature_vmrs(region), tech)
        feat_path = PHASE6 / region / "vmr_features_with_tech.tsv.gz"
        feat.to_csv(feat_path, sep="\t", index=False, compression="gzip")
        print(f"  features+tech n={len(feat)} wrote {feat_path}")

        for feature, label in [
            ("in_LINE_L1", "LINE_L1_enrichment_vs_predictability"),
            ("in_H3K9me3", "H3K9me3_enrichment_vs_predictability"),
            ("in_Quies", "quiescent_chromatin_enrichment_vs_predictability"),
        ]:
            row = analyze_feature(region, feat, feature, args.seed)
            row["analysis"] = label
            all_rows.append(row)
            print(f"  {label}: OR={row.get('original_estimate')} p={row.get('original_pvalue')}")

        brow = analyze_burden(region, args.seed)
        all_rows.append(brow)
        print(f"  burden: coef={brow.get('original_estimate')} p={brow.get('original_pvalue')}")

        # per-region table
        pd.DataFrame([r for r in all_rows if r.get("region") == region]).to_csv(
            PHASE6 / region / "robustness_results.tsv", sep="\t", index=False
        )

    tab = pd.DataFrame(all_rows)
    # consolidated wide table (one row per analysis×region)
    out = PHASE6 / "consolidated_robustness_table.tsv"
    tab.to_csv(out, sep="\t", index=False)

    # manuscript-facing summary: require direction consistency in ≥2 regions for LINE/H3K9/Quies/burden
    summary = []
    for analysis in tab["analysis"].dropna().unique():
        sub = tab[tab["analysis"] == analysis]
        n_ok = int(sub["direction_consistent"].fillna(False).sum()) if "direction_consistent" in sub else 0
        summary.append({
            "analysis": analysis,
            "n_regions_direction_consistent": n_ok,
            "n_regions": len(sub),
            "passes_phase6_claim": n_ok >= 2,
        })
    write_tsv(PHASE6 / "phase6_claim_summary.tsv", summary)
    write_tsv(
        PHASE6 / "sensitivity_plan.tsv",
        [
            {"step": "restrict_high_mappability", "threshold": f">={map_min}", "status": "applied"},
            {"step": "exclude_segmental_duplications", "threshold": "overlaps_segdup==0", "status": "applied"},
            {"step": "exclude_snp_proximal_vmrs", "threshold": "overlap common SNP ±150bp windows", "status": "applied"},
            {"step": "adjust_length_num_snps_or_coverage_variance", "threshold": "GLM covariates", "status": "applied"},
            {"step": "match_high_vs_low_predictability", "threshold": "NN on z-scored technical vars", "status": "applied"},
            {"step": "repeat_LINE_L1_H3K9me3_Quies_burden", "threshold": "n/a", "status": "applied"},
        ],
    )
    print(f"Wrote {out}")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
