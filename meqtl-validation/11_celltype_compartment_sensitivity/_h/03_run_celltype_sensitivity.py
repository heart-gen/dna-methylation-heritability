#!/usr/bin/env python3
"""Cell-composition sensitivity for LINE/L1 / repressive enrichment and burden.

VMR-level: annotation ~ predictability (+ cellPC_r2 / oligo / dnamCellPC_r2).
Sample-level: Oligo / cellPCs ~ mean methylation of LINE vs non-LINE VMRs.
Appends cell-adjusted rows to Phase 6 consolidated robustness table.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import ANALYSIS_SCHEMA_VERSION, write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
MODULE = PROJECT / "meqtl-validation/11_celltype_compartment_sensitivity/_m"
PHASE6 = PROJECT / "meqtl-validation/07_repeat_mappability_sensitivity/_m"
PHASE2 = PROJECT / "meqtl-validation/02_vmr_meqtl_burden/_m"
ANN = PROJECT / "heritability/elastic_net_model/all_individuals/tissue_comparison/annotation"
SEED = 20260807

TISSUE_MAP = {
    "caudate": "Caudate",
    "dlpfc": "DLPFC",
    "hippocampus": "Hippocampus",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions", nargs="+", default=["caudate", "dlpfc", "hippocampus"])
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


def fit_logistic(df: pd.DataFrame, ycol: str, xcols: list[str]) -> dict:
    need = [ycol, "local_predictability"] + xcols
    use = df.dropna(subset=need).copy()
    if use.empty or use[ycol].nunique() < 2:
        return {"estimate": np.nan, "or": np.nan, "pvalue": np.nan, "n": len(use), "error": "insufficient"}
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
        return {"estimate": np.nan, "or": np.nan, "pvalue": np.nan, "n": len(use), "error": str(exc)}


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


def fit_ols(y: pd.Series, X: pd.DataFrame, predictor: str) -> dict:
    df = pd.concat([y.rename("y"), X], axis=1).dropna()
    if df.shape[0] < 20 or df[predictor].nunique() < 2:
        return {"n": int(df.shape[0]), "beta": np.nan, "pvalue": np.nan, "r2": np.nan}
    res = sm.OLS(
        df["y"].astype(float),
        sm.add_constant(df[X.columns].astype(float), has_constant="add"),
    ).fit()
    return {
        "n": int(res.nobs),
        "beta": float(res.params.get(predictor, np.nan)),
        "pvalue": float(res.pvalues.get(predictor, np.nan)),
        "r2": float(res.rsquared),
    }


def load_feature_vmrs(region: str) -> pd.DataFrame:
    tissue = TISSUE_MAP[region]
    rep = pd.read_csv(ANN / "repeat_elements/_m/vmr_repeat_overlap_AA.tsv", sep="\t")
    rep = rep[rep["tissue"] == tissue].copy()
    reprs = pd.read_csv(ANN / "repressive_chromatin/_m/vmr_repressive_overlap_AA.tsv", sep="\t")
    reprs = reprs[reprs["tissue"] == tissue].copy()
    key_cols = ["seqnames", "start", "end"]
    d = rep[key_cols + ["h2_unscaled", "vmr_length", "num_snps", "in_LINE", "in_L1"]].merge(
        reprs[key_cols + ["in_H3K9me3", "in_Quies"]], on=key_cols, how="inner"
    )
    d["coord_id"] = [coord_id(c, s, e) for c, s, e in zip(d["seqnames"], d["start"], d["end"])]
    d["local_predictability"] = pd.to_numeric(d["h2_unscaled"], errors="coerce")
    d["length"] = pd.to_numeric(d["vmr_length"], errors="coerce")
    d["num_snps"] = pd.to_numeric(d["num_snps"], errors="coerce")
    for c in ["in_LINE", "in_L1", "in_H3K9me3", "in_Quies"]:
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0).astype(int).clip(0, 1)
    d["in_LINE_L1"] = ((d["in_LINE"] == 1) | (d["in_L1"] == 1)).astype(int)
    return d


def analyze_region(region: str) -> tuple[list[dict], list[dict], list[dict]]:
    metrics = pd.read_csv(MODULE / region / "vmr_cell_metrics.tsv.gz", sep="\t")
    metrics = metrics[metrics["coord_id"].astype(str).ne("")].copy()
    feat = load_feature_vmrs(region)
    df = feat.merge(
        metrics[
            [
                c
                for c in [
                    "coord_id",
                    "mean_methylation",
                    "var_methylation",
                    "cellPC_r2",
                    "oligo_r2",
                    "abs_oligo_corr",
                    "oligo_corr",
                    "residual_var_cellPC",
                    "dnamCellPC_r2",
                    "high_mappability",
                    "overlaps_segdup",
                    "overlaps_snp_prox",
                ]
                if c in metrics.columns
            ]
        ],
        on="coord_id",
        how="inner",
    )
    df.to_csv(MODULE / region / "vmr_features_with_cell_metrics.tsv.gz", sep="\t", index=False, compression="gzip")

    base_adj = [c for c in ["length", "num_snps"] if c in df.columns]
    model_specs = [
        ("original", []),
        ("adjusted_technical", base_adj),
        ("cellPC_r2_adj", base_adj + ["cellPC_r2"]),
        ("oligo_r2_adj", base_adj + ["oligo_r2"]),
        ("abs_oligo_corr_adj", base_adj + ["abs_oligo_corr"]),
        ("meth_level_var_adj", base_adj + ["mean_methylation", "var_methylation"]),
        ("cellPC_r2_plus_meth", base_adj + ["cellPC_r2", "mean_methylation", "var_methylation"]),
    ]
    if df["dnamCellPC_r2"].notna().sum() > 100:
        model_specs.append(("dnamCellPC_r2_adj", base_adj + ["dnamCellPC_r2"]))

    enrich_rows = []
    for feature, label in [
        ("in_LINE_L1", "LINE_L1_enrichment_vs_predictability"),
        ("in_H3K9me3", "H3K9me3_enrichment_vs_predictability"),
        ("in_Quies", "quiescent_chromatin_enrichment_vs_predictability"),
    ]:
        for model_name, xcols in model_specs:
            r = fit_logistic(df, feature, xcols)
            enrich_rows.append(
                {
                    "region": region,
                    "analysis": label,
                    "model": model_name,
                    "or": r["or"],
                    "estimate": r["estimate"],
                    "pvalue": r["pvalue"],
                    "n": r["n"],
                    "n_feature": r.get("n_feature", np.nan),
                    "error": r.get("error", ""),
                }
            )
            # high-mappability subset for key models
            if model_name in {"original", "adjusted_technical", "cellPC_r2_adj", "dnamCellPC_r2_adj"}:
                if "high_mappability" in df.columns:
                    sub = df[df["high_mappability"].fillna(0).astype(int) == 1]
                    rr = fit_logistic(sub, feature, xcols)
                    enrich_rows.append(
                        {
                            "region": region,
                            "analysis": label,
                            "model": f"{model_name}__high_mappability",
                            "or": rr["or"],
                            "estimate": rr["estimate"],
                            "pvalue": rr["pvalue"],
                            "n": rr["n"],
                            "n_feature": rr.get("n_feature", np.nan),
                            "error": rr.get("error", ""),
                        }
                    )

    # Burden
    burden_path = PHASE2 / region / "vmr_meqtl_burden.tsv.gz"
    burden_rows = []
    if burden_path.exists():
        b = pd.read_csv(burden_path, sep="\t")
        if "analysis_schema_version" not in b or not b["analysis_schema_version"].eq(ANALYSIS_SCHEMA_VERSION).all():
            raise SystemExit(f"{region}: stale Phase 2 burden; regenerate repair schema v2 before cell-type sensitivity")
        b["vmr_id"] = b["vmr_id"].astype(str)
        # join metrics via task_id / coord
        tech = pd.read_csv(PHASE6 / region / "vmr_technical_annotations.tsv", sep="\t")
        tech["coord_id"] = [
            coord_id(c, s, e) for c, s, e in zip(tech["chrom"], tech["start"], tech["end"])
        ]
        mjoin = metrics.copy()
        # task_id path
        if "task_id" in tech.columns:
            tmap = tech.dropna(subset=["task_id"])[["task_id", "coord_id"]].copy()
            tmap["vmr_id"] = tmap["task_id"].astype(float).astype(int).astype(str)
            bm = b.merge(tmap[["vmr_id", "coord_id"]], on="vmr_id", how="left")
        else:
            bm = b.copy()
            bm["coord_id"] = bm["vmr_id"]
        bm = bm.merge(
            mjoin[
                [
                    c
                    for c in ["coord_id", "cellPC_r2", "oligo_r2", "abs_oligo_corr", "dnamCellPC_r2"]
                    if c in mjoin.columns
                ]
            ],
            on="coord_id",
            how="left",
        )
        tech_adj = [c for c in ["n_tested_cpgs", "average_cpg_coverage", "mean_cpg_variance"] if c in bm.columns]
        for model_name, xcols in [
            ("original", ["n_tested_cpgs"] if "n_tested_cpgs" in bm.columns else []),
            ("adjusted_technical", tech_adj),
            ("cellPC_r2_adj", tech_adj + ["cellPC_r2"]),
            ("oligo_r2_adj", tech_adj + ["oligo_r2"]),
        ]:
            if "dnamCellPC_r2" in xcols or model_name.startswith("dnam"):
                pass
            r = fit_burden(bm, xcols)
            burden_rows.append(
                {
                    "region": region,
                    "analysis": "predictability_meqtl_burden_association",
                    "model": model_name,
                    "estimate": r["estimate"],
                    "or": np.nan,
                    "pvalue": r["pvalue"],
                    "n": r["n"],
                    "error": r.get("error", ""),
                }
            )
        if bm["dnamCellPC_r2"].notna().sum() > 100:
            r = fit_burden(bm, tech_adj + ["dnamCellPC_r2"])
            burden_rows.append(
                {
                    "region": region,
                    "analysis": "predictability_meqtl_burden_association",
                    "model": "dnamCellPC_r2_adj",
                    "estimate": r["estimate"],
                    "or": np.nan,
                    "pvalue": r["pvalue"],
                    "n": r["n"],
                    "error": r.get("error", ""),
                }
            )

    # Sample-level Oligo checks
    meth = pd.read_csv(MODULE / region / "vmr_mean_methylation.tsv.gz", sep="\t", index_col=0)
    scov = pd.read_csv(MODULE / region / "sample_cell_covariates.tsv", sep="\t")
    scov["id"] = scov["id"].astype(str)
    scov = scov.set_index("id")

    # map meth columns to LINE status via metrics/features
    key_line = {}
    for _, r in df.iterrows():
        key_line[str(r["coord_id"])] = int(r["in_LINE_L1"])
    # also task keys from metrics
    met = metrics.copy()
    for _, r in met.iterrows():
        cid = str(r["coord_id"])
        if cid in key_line:
            key_line[str(r["vmr_key"])] = key_line[cid]

    line_cols = [c for c in meth.columns if key_line.get(str(c), 0) == 1]
    non_cols = [c for c in meth.columns if key_line.get(str(c), -1) == 0]
    sample_rows = []
    if line_cols and non_cols:
        y_line = meth[line_cols].mean(axis=1, skipna=True)
        y_non = meth[non_cols].mean(axis=1, skipna=True)
        y_diff = y_line - y_non
        for outcome_name, y in [
            ("mean_meth_LINE_L1_vmrs", y_line),
            ("mean_meth_nonLINE_vmrs", y_non),
            ("mean_meth_LINE_minus_nonLINE", y_diff),
        ]:
            y = y.reindex(scov.index)
            for pred in ["Oligo", "cellPC1", "neuron_agg"]:
                if pred not in scov.columns:
                    continue
                fit = fit_ols(y, scov[[pred]], pred)
                sample_rows.append(
                    {
                        "region": region,
                        "outcome": outcome_name,
                        "predictor": pred,
                        "n_line_vmrs": len(line_cols),
                        "n_nonline_vmrs": len(non_cols),
                        **fit,
                    }
                )
            # multivariable cellPCs
            pcs = [c for c in ["cellPC1", "cellPC2", "cellPC3"] if c in scov.columns]
            if pcs:
                dfp = pd.concat([y.rename("y"), scov[pcs]], axis=1).dropna()
                if len(dfp) >= 20:
                    res = sm.OLS(
                        dfp["y"], sm.add_constant(dfp[pcs].astype(float), has_constant="add")
                    ).fit()
                    # joint F for cellPCs
                    sample_rows.append(
                        {
                            "region": region,
                            "outcome": outcome_name,
                            "predictor": "cellPC1-3_joint",
                            "n_line_vmrs": len(line_cols),
                            "n_nonline_vmrs": len(non_cols),
                            "n": int(res.nobs),
                            "beta": float(res.params.get("cellPC1", np.nan)),
                            "pvalue": float(res.f_pvalue),
                            "r2": float(res.rsquared),
                        }
                    )

    return enrich_rows, burden_rows, sample_rows


def build_consolidated_append(enrich: pd.DataFrame, burden: pd.DataFrame) -> pd.DataFrame:
    """Wide rows compatible with Phase 6 consolidated table (+ cell columns)."""
    rows = []
    analyses = sorted(set(enrich["analysis"]).union(set(burden["analysis"])))
    for region in enrich["region"].unique():
        for analysis in analyses:
            if analysis == "predictability_meqtl_burden_association":
                src = burden[(burden["region"] == region) & (burden["analysis"] == analysis)]
                get = lambda m: src[src["model"] == m]
                val = lambda m, col: (
                    float(get(m).iloc[0][col]) if len(get(m)) and pd.notna(get(m).iloc[0][col]) else np.nan
                )
                row = {
                    "region": region,
                    "analysis": analysis,
                    "status": "ok",
                    "original_estimate": val("original", "estimate"),
                    "original_pvalue": val("original", "pvalue"),
                    "original_n": val("original", "n"),
                    "adjusted_estimate": val("adjusted_technical", "estimate"),
                    "adjusted_pvalue": val("adjusted_technical", "pvalue"),
                    "cellPC_adj_estimate": val("cellPC_r2_adj", "estimate"),
                    "cellPC_adj_pvalue": val("cellPC_r2_adj", "pvalue"),
                    "dnamCellPC_adj_estimate": val("dnamCellPC_r2_adj", "estimate"),
                    "dnamCellPC_adj_pvalue": val("dnamCellPC_r2_adj", "pvalue"),
                    "oligo_r2_adj_estimate": val("oligo_r2_adj", "estimate"),
                    "oligo_r2_adj_pvalue": val("oligo_r2_adj", "pvalue"),
                }
            else:
                src = enrich[(enrich["region"] == region) & (enrich["analysis"] == analysis)]
                get = lambda m: src[src["model"] == m]
                val = lambda m, col: (
                    float(get(m).iloc[0][col]) if len(get(m)) and pd.notna(get(m).iloc[0][col]) else np.nan
                )
                row = {
                    "region": region,
                    "analysis": analysis,
                    "status": "ok",
                    "original_estimate": val("original", "or"),
                    "original_pvalue": val("original", "pvalue"),
                    "original_n": val("original", "n"),
                    "adjusted_estimate": val("adjusted_technical", "or"),
                    "adjusted_pvalue": val("adjusted_technical", "pvalue"),
                    "cellPC_adj_estimate": val("cellPC_r2_adj", "or"),
                    "cellPC_adj_pvalue": val("cellPC_r2_adj", "pvalue"),
                    "dnamCellPC_adj_estimate": val("dnamCellPC_r2_adj", "or"),
                    "dnamCellPC_adj_pvalue": val("dnamCellPC_r2_adj", "pvalue"),
                    "oligo_r2_adj_estimate": val("oligo_r2_adj", "or"),
                    "oligo_r2_adj_pvalue": val("oligo_r2_adj", "pvalue"),
                    "high_mappability_cellPC_estimate": val(
                        "cellPC_r2_adj__high_mappability", "or"
                    ),
                    "high_mappability_cellPC_pvalue": val(
                        "cellPC_r2_adj__high_mappability", "pvalue"
                    ),
                }
            # direction: OR>1 or coef>0
            ests = [
                row.get("original_estimate"),
                row.get("adjusted_estimate"),
                row.get("cellPC_adj_estimate"),
            ]
            if analysis != "predictability_meqtl_burden_association":
                signs = []
                for e in ests:
                    if pd.notna(e) and e > 0:
                        signs.append(np.sign(np.log(e)) if e != 1 else 0)
                row["direction_consistent_with_cellPC"] = bool(
                    len(signs) >= 2 and len(set(signs) - {0}) <= 1 and (not signs or signs[0] > 0)
                )
            else:
                row["direction_consistent_with_cellPC"] = bool(
                    all(pd.isna(e) or e > 0 for e in ests) and any(pd.notna(e) and e > 0 for e in ests)
                )
            # attenuation: cellPC OR / original OR
            o, c = row.get("original_estimate"), row.get("cellPC_adj_estimate")
            if pd.notna(o) and pd.notna(c) and o not in (0, 1) and analysis != "predictability_meqtl_burden_association":
                row["cellPC_attenuation_ratio_logOR"] = float(np.log(c) / np.log(o)) if np.log(o) != 0 else np.nan
            elif pd.notna(o) and pd.notna(c) and analysis == "predictability_meqtl_burden_association":
                row["cellPC_attenuation_ratio_logOR"] = float(c / o) if o != 0 else np.nan
            else:
                row["cellPC_attenuation_ratio_logOR"] = np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def claim_decision(wide: pd.DataFrame, sample_df: pd.DataFrame) -> list[dict]:
    claims = []
    for analysis in [
        "LINE_L1_enrichment_vs_predictability",
        "H3K9me3_enrichment_vs_predictability",
        "quiescent_chromatin_enrichment_vs_predictability",
        "predictability_meqtl_burden_association",
    ]:
        sub = wide[wide["analysis"] == analysis]
        n_dir = int(sub["direction_consistent_with_cellPC"].fillna(False).sum())
        # significance of cellPC_adj in ≥2 regions
        n_sig = int((pd.to_numeric(sub["cellPC_adj_pvalue"], errors="coerce") < 0.05).sum())
        caud = sub[sub["region"] == "caudate"]
        caud_pass = False
        if len(caud):
            r = caud.iloc[0]
            est = r.get("cellPC_adj_estimate")
            p = r.get("cellPC_adj_pvalue")
            if analysis == "predictability_meqtl_burden_association":
                caud_pass = bool(pd.notna(est) and est > 0 and pd.notna(p) and p < 0.05)
            else:
                caud_pass = bool(pd.notna(est) and est > 1 and pd.notna(p) and p < 0.05)
        passes = n_dir >= 2 and n_sig >= 2 and caud_pass
        claims.append(
            {
                "analysis": analysis,
                "n_regions_direction_consistent_cellPC": n_dir,
                "n_regions_cellPC_adj_pvalue_lt_0.05": n_sig,
                "caudate_cellPC_adj_pass": caud_pass,
                "passes_celltype_claim": passes,
            }
        )

    # sample-level: Oligo association with LINE meth should not fully explain architecture
    oligo_line = sample_df[
        (sample_df["outcome"] == "mean_meth_LINE_L1_vmrs") & (sample_df["predictor"] == "Oligo")
    ]
    claims.append(
        {
            "analysis": "sample_level_Oligo_vs_LINE_mean_meth",
            "n_regions_direction_consistent_cellPC": int(
                (pd.to_numeric(oligo_line["pvalue"], errors="coerce") < 0.05).sum()
            ),
            "n_regions_cellPC_adj_pvalue_lt_0.05": int(len(oligo_line)),
            "caudate_cellPC_adj_pass": "",
            "passes_celltype_claim": "",
            "note": "descriptive confounding check; significant Oligo~LINE-meth does not alone refute genomic enrichment",
        }
    )
    return claims


def main() -> None:
    args = parse_args()
    all_enrich, all_burden, all_sample = [], [], []
    for region in args.regions:
        print(f"==== {region} ====")
        e, b, s = analyze_region(region)
        all_enrich.extend(e)
        all_burden.extend(b)
        all_sample.extend(s)
        print(f"  enrich models={len(e)} burden={len(b)} sample={len(s)}")

    enrich = pd.DataFrame(all_enrich)
    burden = pd.DataFrame(all_burden)
    sample_df = pd.DataFrame(all_sample)
    enrich.to_csv(MODULE / "enrichment_celltype_models.tsv", sep="\t", index=False)
    burden.to_csv(MODULE / "burden_celltype_models.tsv", sep="\t", index=False)
    sample_df.to_csv(MODULE / "sample_level_celltype_checks.tsv", sep="\t", index=False)

    wide = build_consolidated_append(enrich, burden)
    wide.to_csv(MODULE / "celltype_robustness_table.tsv", sep="\t", index=False)

    # Merge into Phase 6 consolidated table (add columns; keep original rows)
    p6_path = PHASE6 / "consolidated_robustness_table.tsv"
    p6 = pd.read_csv(p6_path, sep="\t")
    merge_cols = [
        "region",
        "analysis",
        "cellPC_adj_estimate",
        "cellPC_adj_pvalue",
        "dnamCellPC_adj_estimate",
        "dnamCellPC_adj_pvalue",
        "oligo_r2_adj_estimate",
        "oligo_r2_adj_pvalue",
        "high_mappability_cellPC_estimate",
        "high_mappability_cellPC_pvalue",
        "direction_consistent_with_cellPC",
        "cellPC_attenuation_ratio_logOR",
    ]
    # drop prior cell columns if re-run
    drop_prior = [c for c in merge_cols if c not in {"region", "analysis"} and c in p6.columns]
    if drop_prior:
        p6 = p6.drop(columns=drop_prior)
    merged = p6.merge(wide[merge_cols], on=["region", "analysis"], how="left")
    merged.to_csv(p6_path, sep="\t", index=False)
    # also keep a dated copy in module 11
    merged.to_csv(MODULE / "consolidated_robustness_table_with_celltype.tsv", sep="\t", index=False)

    claims = claim_decision(wide, sample_df)
    write_tsv(MODULE / "celltype_claim_summary.tsv", claims)

    # Go/no-go language
    line = next(c for c in claims if c["analysis"] == "LINE_L1_enrichment_vs_predictability")
    if line["passes_celltype_claim"]:
        decision = "keep_main_figure_with_cell_composition_row"
        language = (
            "LINE/L1–repressive enrichment of high-predictability VMRs remains after adjustment "
            "for cell-composition–correlated methylation properties in caudate and ≥1 other region."
        )
    else:
        caud_only = bool(line["caudate_cellPC_adj_pass"]) and int(line["n_regions_direction_consistent_cellPC"]) < 2
        if caud_only:
            decision = "main_text_caudate_others_supplement"
            language = (
                "Cell-adjusted LINE/L1 enrichment is supported in caudate; cross-region consistency "
                "is limited — keep caudate in main text and treat other regions as supplemental."
            )
        else:
            # check attenuation to null
            line_wide = wide[wide["analysis"] == "LINE_L1_enrichment_vs_predictability"]
            n_null = int(
                (
                    (pd.to_numeric(line_wide["cellPC_adj_pvalue"], errors="coerce") >= 0.05)
                    | (pd.to_numeric(line_wide["cellPC_adj_estimate"], errors="coerce") <= 1)
                ).sum()
            )
            if n_null >= 2:
                decision = "downgrade_composition_not_excluded"
                language = (
                    "Bulk compartment correlation; cell composition not excluded as a driver of "
                    "LINE/L1 enrichment after cellPC_r2 adjustment."
                )
            else:
                decision = "mixed_report_with_caveat"
                language = (
                    "Cell-composition adjustment yields mixed regional results; report with explicit "
                    "bulk-cellularity caveat."
                )

    snap = {
        "decision": decision,
        "allowed_language": language,
        "line_l1_passes_celltype_claim": line["passes_celltype_claim"],
        "n_regions_line_direction_cellPC": line["n_regions_direction_consistent_cellPC"],
        "primary_adjustment": "cellPC_r2 from meth ~ MuSiC cellPC1-3",
        "dnam_primary_region": "caudate (M6d / integration gate)",
    }
    write_tsv(MODULE / "celltype_claim_snapshot.tsv", [snap])

    md = MODULE / "CELLTYPE_RESULTS.md"
    md.write_text(
        f"""# Cell-type LINE/L1 sensitivity — results

**Decision:** `{decision}`

{language}

## Claim summary

| Analysis | Regions dir. consistent (cellPC) | Regions p&lt;0.05 | Caudate pass | Claim pass |
|---|---:|---:|---|---|
"""
        + "\n".join(
            f"| {c['analysis']} | {c['n_regions_direction_consistent_cellPC']} | "
            f"{c['n_regions_cellPC_adj_pvalue_lt_0.05']} | {c['caudate_cellPC_adj_pass']} | "
            f"{c['passes_celltype_claim']} |"
            for c in claims
            if c["analysis"] != "sample_level_Oligo_vs_LINE_mean_meth"
        )
        + f"""

## Key outputs

- `{MODULE / 'enrichment_celltype_models.tsv'}`
- `{MODULE / 'sample_level_celltype_checks.tsv'}`
- `{MODULE / 'celltype_robustness_table.tsv'}`
- Phase 6 table updated: `{p6_path}`

## Interpretation rules

- Genomic LINE/L1 overlap is unchanged by deconvolution; adjustment targets methylation variance explained by cell composition as a confounder of the predictability ranking.
- Significant sample-level Oligo ~ LINE-VMR methylation indicates bulk cellularity tracks repeat-rich methylation, but does not by itself refute sequence/compartment enrichment if cellPC-adjusted ORs remain &gt;1.
"""
    )
    print(f"Decision: {decision}")
    print(pd.DataFrame(claims).to_string(index=False))


if __name__ == "__main__":
    main()
