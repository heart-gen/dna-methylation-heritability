#!/usr/bin/env python3
"""Complete LINE/mappability tech joins onto Phase 2 VMR burden tables.

Root cause: burden uses elastic-net task_id keys, while technical annotations are
interval-based and only carry task_id for exact chrom/start/end matches (~27–45%).
This script joins by genomic overlap (best reciprocal overlap) and rewrites
burden tech columns + a join completeness report, then refreshes Phase 6 burden
sensitivity rows in the consolidated robustness table.
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
PHASE2 = PROJECT / "meqtl-validation/02_vmr_meqtl_burden/_m"
PHASE6 = PROJECT / "meqtl-validation/07_repeat_mappability_sensitivity/_m"
TECH_ROOT = PHASE6
JOIN_OUTPUT_ROOT = PHASE6
EN_ROOT = PROJECT / "heritability/elastic_net_model/all_individuals"
SEED = 20260807
BEDTOOLS = Path("/projects/p32505/opt/envs/genomics/bin/bedtools")

TECH_COLS = [
    "umap_k24_mean",
    "high_mappability",
    "line_l1_frac",
    "segdup_frac",
    "overlaps_segdup",
    "blacklist_frac",
    "overlaps_blacklist",
    "snp_prox_frac",
    "overlaps_snp_prox",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions", nargs="+", default=["caudate", "dlpfc", "hippocampus"])
    p.add_argument("--join-only", action="store_true", help="Join annotations but do not update Phase 6 summaries")
    p.add_argument("--min-reciprocal-overlap", type=float, default=0.5)
    p.add_argument("--phase2-root", default=str(PHASE2), help="Burden-table root; supports immutable run directories")
    p.add_argument("--technical-root", default=str(TECH_ROOT), help="Input technical-annotation root")
    p.add_argument("--join-output-root", default=str(JOIN_OUTPUT_ROOT), help="Output root for join tables and completeness report")
    p.add_argument("--elastic-net-root", default=str(EN_ROOT), help="Legacy elastic-net summary root")
    return p.parse_args()


def _z(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x * 0.0
    return (x - x.mean()) / sd


def fit_burden(df: pd.DataFrame, xcols: list[str]) -> dict:
    need = ["local_predictability", "n_cpgs_with_sig_meqtl", "n_tested_cpgs"] + xcols
    d = df.dropna(subset=[c for c in need if c in df.columns]).copy()
    d = d[d["n_tested_cpgs"] > 0]
    if d.empty:
        return {"estimate": np.nan, "pvalue": np.nan, "n": 0}
    d["n_ns"] = d["n_tested_cpgs"] - d["n_cpgs_with_sig_meqtl"]
    cols = ["local_predictability"] + [c for c in xcols if c in d.columns]
    use = d.dropna(subset=cols)
    endog = use[["n_cpgs_with_sig_meqtl", "n_ns"]]
    exog = sm.add_constant(use[cols].apply(_z), has_constant="add")
    model = sm.GLM(endog, exog, family=sm.families.Binomial())
    pilot = model.fit(maxiter=250)
    dispersion = max(1.0, float(pilot.pearson_chi2 / max(pilot.df_resid, 1)))
    res = model.fit(scale=dispersion, cov_type="HC3", maxiter=250)
    return {
        "estimate": float(res.params["local_predictability"]),
        "pvalue": float(res.pvalues["local_predictability"]),
        "n": int(len(use)),
    }


def matched_delta(df: pd.DataFrame, ycol: str, match_cols: list[str], seed: int) -> dict:
    need = ["local_predictability", ycol] + match_cols
    d = df.dropna(subset=need).copy()
    if len(d) < 40:
        return {"estimate": np.nan, "pvalue": np.nan, "n": len(d)}
    thresholds = load_yaml("analysis_thresholds.yml")["matching"]
    pairs, _balance, _meta = greedy_nearest_neighbor_pairs(
        d,
        exposure="local_predictability",
        outcome=ycol,
        numeric_covariates=match_cols,
        caliper_sd=float(thresholds["caliper_sd"]),
        seed=seed,
    )
    if pairs.empty:
        return {"estimate": np.nan, "pvalue": np.nan, "n": 0}
    differences = pairs["outcome_high"] - pairs["outcome_low"]
    p = paired_randomization_pvalue(differences, seed=seed, n_perm=10000)
    return {"estimate": float(differences.mean()), "pvalue": float(p), "n": int(len(pairs))}


def bedtools_best_overlap(
    query: pd.DataFrame, target: pd.DataFrame, min_reciprocal_overlap: float
) -> pd.DataFrame:
    """Return one best-overlap target name per query name."""
    with tempfile.TemporaryDirectory(prefix="tech_join_") as tmp:
        tmp = Path(tmp)
        q = query.copy()
        t = target.copy()
        if "qname" not in q.columns:
            q["qname"] = [f"q{i}" for i in range(len(q))]
        if "tname" not in t.columns:
            t["tname"] = [f"t{i}" for i in range(len(t))]
        qbed = tmp / "q.bed"
        tbed = tmp / "t.bed"
        q[["chrom", "start", "end", "qname"]].to_csv(qbed, sep="\t", header=False, index=False)
        t[["chrom", "start", "end", "tname"]].to_csv(tbed, sep="\t", header=False, index=False)
        out = tmp / "overlap.tsv"
        cmd = [
            str(BEDTOOLS),
            "intersect",
            "-a",
            str(qbed),
            "-b",
            str(tbed),
            "-wo",
        ]
        with out.open("w") as handle:
            subprocess.check_call(cmd, stdout=handle)
        if out.stat().st_size == 0:
            return pd.DataFrame(columns=[
                "qname", "tname", "overlap_bp", "frac_query", "frac_target",
                "min_reciprocal_frac",
            ])
        ov = pd.read_csv(
            out,
            sep="\t",
            header=None,
            names=[
                "q_chrom",
                "q_start",
                "q_end",
                "qname",
                "t_chrom",
                "t_start",
                "t_end",
                "tname",
                "overlap_bp",
            ],
        )
        qlen = (ov["q_end"] - ov["q_start"]).clip(lower=1)
        tlen = (ov["t_end"] - ov["t_start"]).clip(lower=1)
        ov["frac_query"] = ov["overlap_bp"] / qlen
        ov["frac_target"] = ov["overlap_bp"] / tlen
        ov = ov[
            (ov["frac_query"] >= min_reciprocal_overlap)
            & (ov["frac_target"] >= min_reciprocal_overlap)
        ]
        ov["min_reciprocal_frac"] = ov[["frac_query", "frac_target"]].min(axis=1)
        ov = ov.sort_values(
            ["qname", "min_reciprocal_frac", "overlap_bp"],
            ascending=[True, False, False],
        )
        best = ov.drop_duplicates("qname", keep="first")[[
            "qname", "tname", "overlap_bp", "frac_query", "frac_target",
            "min_reciprocal_frac",
        ]]
        return best


def complete_region(region: str, min_reciprocal_overlap: float) -> dict:
    bur_path = PHASE2 / region / "vmr_meqtl_burden.tsv.gz"
    tech_path = TECH_ROOT / region / "vmr_technical_annotations.tsv"
    en_path = EN_ROOT / region / "_m" / f"{region}_summary_elastic-net_AA.tsv"
    bur = pd.read_csv(bur_path, sep="\t")
    if "analysis_schema_version" not in bur or not bur["analysis_schema_version"].eq(ANALYSIS_SCHEMA_VERSION).all():
        raise SystemExit(f"{region}: stale Phase 2 burden; regenerate repair schema v2 before Phase 6")
    tech = pd.read_csv(tech_path, sep="\t")
    en = pd.read_csv(en_path, sep="\t")
    bur["vmr_id"] = bur["vmr_id"].astype(str)
    en["task_id"] = en["task_id"].astype(str)
    en["chrom"] = en["chrom"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")
    tech["chrom"] = tech["chrom"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")

    before = {
        c: float(bur[c].notna().mean()) if c in bur.columns else 0.0 for c in TECH_COLS
    }
    prior_tech = bur[["vmr_id"] + [c for c in TECH_COLS if c in bur.columns]].copy()

    # Drop prior incomplete tech columns then rebuild via multi-key coalesce
    for c in TECH_COLS:
        if c in bur.columns:
            bur = bur.drop(columns=[c])
    for c in [
        "tech_join_source", "tech_join_query_overlap_frac",
        "tech_join_target_overlap_frac", "tech_join_min_reciprocal_frac",
    ]:
        if c in bur.columns:
            bur = bur.drop(columns=[c])

    t = tech.copy()
    t["chrom"] = t["chrom"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")
    t["start"] = t["start"].astype(int)
    t["end"] = t["end"].astype(int)
    t["task_id"] = t["task_id"].apply(
        lambda x: str(int(float(x))) if pd.notna(x) else pd.NA
    )
    if "interval_id" not in t.columns:
        t["interval_id"] = (
            t["chrom"].str.replace("^chr", "", regex=True)
            + ":"
            + t["start"].astype(str)
            + "-"
            + t["end"].astype(str)
        )
    t_cols = [c for c in TECH_COLS if c in t.columns]

    # Strategy 1: direct task_id
    s1 = (
        t.dropna(subset=["task_id"])[["task_id"] + t_cols]
        .drop_duplicates("task_id", keep="first")
        .rename(columns={"task_id": "vmr_id"})
    )
    s1["join_source"] = "task_id"

    # Strategy 2: elastic-net coords -> exact interval_id
    en2 = en.copy()
    en2["interval_id"] = (
        en2["chrom"].str.replace("^chr", "", regex=True)
        + ":"
        + en2["start"].astype(int).astype(str)
        + "-"
        + en2["end"].astype(int).astype(str)
    )
    s2 = bur[["vmr_id"]].merge(en2[["task_id", "interval_id"]], left_on="vmr_id", right_on="task_id", how="inner")
    s2 = s2.merge(
        t[["interval_id"] + t_cols].drop_duplicates("interval_id", keep="first"),
        on="interval_id",
        how="inner",
    )
    s2 = s2[["vmr_id"] + t_cols].drop_duplicates("vmr_id")
    s2["join_source"] = "interval_id"

    # Strategy 3: bedtools best overlap for remaining
    q = bur.merge(
        en[["task_id", "chrom", "start", "end"]],
        left_on="vmr_id",
        right_on="task_id",
        how="left",
    )
    n_with_coords = int(q["start"].notna().sum())
    have = set(s1["vmr_id"]).union(set(s2["vmr_id"]))
    q_ok = q.dropna(subset=["chrom", "start", "end"]).copy()
    q_ok = q_ok[~q_ok["vmr_id"].isin(have)].copy()
    q_ok["start"] = q_ok["start"].astype(int)
    q_ok["end"] = q_ok["end"].astype(int)
    q_ok = q_ok.reset_index(drop=True)
    q_ok["qname"] = [f"q{i}" for i in range(len(q_ok))]
    t = t.reset_index(drop=True)
    t["tname"] = [f"t{i}" for i in range(len(t))]
    if len(q_ok):
        best = bedtools_best_overlap(
            q_ok[["chrom", "start", "end", "qname"]],
            t[["chrom", "start", "end", "tname"]],
            min_reciprocal_overlap,
        )
        linked = q_ok[["vmr_id", "qname"]].merge(best, on="qname", how="left")
        t_tech = t.set_index("tname")[t_cols]
        linked = linked.merge(t_tech, left_on="tname", right_index=True, how="left")
        if "frac_query" in linked.columns:
            linked = linked.sort_values("frac_query", ascending=False).drop_duplicates("vmr_id")
            s3 = linked[[
                "vmr_id", *t_cols, "frac_query", "frac_target", "min_reciprocal_frac"
            ]].copy()
        else:
            linked = linked.drop_duplicates("vmr_id")
            s3 = linked[["vmr_id"] + t_cols].copy()
        s3["join_source"] = "overlap"
    else:
        s3 = pd.DataFrame(columns=["vmr_id"] + t_cols + ["join_source"])

    # Coalesce: prefer task_id > interval_id > overlap
    tech_join = s1.copy()
    for add in (s2, s3):
        if add.empty:
            continue
        new_rows = add[~add["vmr_id"].isin(set(tech_join["vmr_id"]))]
        tech_join = pd.concat([tech_join, new_rows], ignore_index=True, sort=False)

    tech_join = tech_join.rename(columns={"join_source": "tech_join_source"})
    out = bur.merge(tech_join, on="vmr_id", how="left")
    # Fill gaps from prior only when prior increases coverage (protect against regressions)
    if "umap_k24_mean" in prior_tech.columns:
        prior_n = int(prior_tech["umap_k24_mean"].notna().sum())
        new_n = int(out["umap_k24_mean"].notna().sum()) if "umap_k24_mean" in out.columns else 0
        if prior_n > new_n:
            out = out.merge(prior_tech, on="vmr_id", how="left", suffixes=("", "_prior"))
            for c in t_cols:
                prior_c = f"{c}_prior"
                if prior_c in out.columns:
                    out[c] = out[c].where(out[c].notna(), out[prior_c])
                    out = out.drop(columns=[prior_c])
    if "high_mappability" in out.columns and "umap_k24_mean" in out.columns:
        miss_hm = out["high_mappability"].isna() & out["umap_k24_mean"].notna()
        out.loc[miss_hm, "high_mappability"] = (out.loc[miss_hm, "umap_k24_mean"] >= 0.9).astype(int)

    after = {
        c: float(out[c].notna().mean()) if c in out.columns else 0.0 for c in TECH_COLS
    }
    n_joined = int(out["umap_k24_mean"].notna().sum()) if "umap_k24_mean" in out.columns else 0

    # archive previous burden then write
    archive = PHASE2 / region / "vmr_meqtl_burden.before_tech_join_fix.tsv.gz"
    if not archive.exists():
        bur.to_csv(archive, sep="\t", index=False, compression="gzip")
    out_path = bur_path
    # drop helper frac from main burden unless useful
    if "frac_query" in out.columns:
        out = out.rename(columns={
            "frac_query": "tech_join_query_overlap_frac",
            "frac_target": "tech_join_target_overlap_frac",
            "min_reciprocal_frac": "tech_join_min_reciprocal_frac",
        })
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")

    # Also write a dedicated tech-joined annotation for Phase 6
    join_out = JOIN_OUTPUT_ROOT / region / "burden_tech_join.tsv.gz"
    join_out.parent.mkdir(parents=True, exist_ok=True)
    out[
        ["vmr_id"]
        + [c for c in TECH_COLS if c in out.columns]
        + [c for c in [
            "tech_join_source",
            "tech_join_query_overlap_frac", "tech_join_target_overlap_frac",
            "tech_join_min_reciprocal_frac",
        ] if c in out.columns]
    ].to_csv(join_out, sep="\t", index=False, compression="gzip")

    report = {
        "region": region,
        "n_burden": int(len(bur)),
        "n_elastic_net": int(len(en)),
        "n_tech_intervals": int(len(tech)),
        "n_burden_with_en_coords": n_with_coords,
        "n_burden_with_umap_before": int(round(before.get("umap_k24_mean", 0) * len(bur))),
        "frac_umap_before": before.get("umap_k24_mean", 0.0),
        "n_burden_with_umap_after": n_joined,
        "frac_umap_after": after.get("umap_k24_mean", 0.0),
        "frac_line_l1_after": after.get("line_l1_frac", 0.0),
        "frac_snp_prox_after": after.get("overlaps_snp_prox", 0.0),
        "delta_frac_umap": after.get("umap_k24_mean", 0.0) - before.get("umap_k24_mean", 0.0),
        "min_reciprocal_overlap_required": min_reciprocal_overlap,
        "burden_path": str(out_path),
        "archive_path": str(archive),
    }

    # Refresh Phase 6 burden robustness row for this region
    row = analyze_burden_row(region, out)
    return report, row


def analyze_burden_row(region: str, merged: pd.DataFrame) -> dict:
    if "overlaps_snp_prox" not in merged.columns and "snp_prox_frac" in merged.columns:
        merged = merged.copy()
        merged["overlaps_snp_prox"] = (merged["snp_prox_frac"].fillna(0) > 0).astype(int)
    adj = [c for c in ["n_tested_cpgs", "average_cpg_coverage", "mean_cpg_variance"] if c in merged.columns]
    subsets = {"original": merged}
    if "high_mappability" in merged.columns:
        subsets["high_mappability"] = merged[merged["high_mappability"].fillna(0).astype(int) == 1]
    if "overlaps_segdup" in merged.columns:
        subsets["segdup_excluded"] = merged[merged["overlaps_segdup"].fillna(0).astype(int) == 0]
    if "overlaps_snp_prox" in merged.columns:
        subsets["snp_proximity_excluded"] = merged[merged["overlaps_snp_prox"].fillna(0).astype(int) == 0]

    row = {"region": region, "analysis": "predictability_meqtl_burden_association", "status": "ok"}
    o = fit_burden(subsets["original"], ["n_tested_cpgs"] if "n_tested_cpgs" in merged.columns else [])
    row["original_estimate"] = o["estimate"]
    row["original_pvalue"] = o["pvalue"]
    row["original_n"] = o["n"]
    a = fit_burden(subsets["original"], adj)
    row["adjusted_estimate"] = a["estimate"]
    row["adjusted_pvalue"] = a["pvalue"]
    mcols = [
        c
        for c in ["n_tested_cpgs", "average_cpg_coverage", "mean_cpg_variance", "umap_k24_mean"]
        if c in merged.columns and merged[c].notna().sum() > 50
    ]
    if "proportion_cpgs_with_sig_meqtl" in merged.columns and mcols:
        m = matched_delta(
            subsets["original"].assign(_y=subsets["original"]["proportion_cpgs_with_sig_meqtl"]),
            "_y",
            mcols,
            SEED,
        )
    else:
        m = {"estimate": np.nan, "pvalue": np.nan}
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
            row[f"{key}_n"] = r["n"]
        else:
            row[col_est] = np.nan
            row[col_p] = np.nan
    ests = [
        row.get(k)
        for k in [
            "original_estimate",
            "adjusted_estimate",
            "matched_estimate",
            "high_mappability_estimate",
            "snp_proximity_excluded_estimate",
            "segdup_excluded_estimate",
        ]
        if pd.notna(row.get(k))
    ]
    row["direction_consistent"] = bool(ests and all(e > 0 for e in ests))
    return row


def main() -> None:
    global PHASE2, TECH_ROOT, JOIN_OUTPUT_ROOT, EN_ROOT
    args = parse_args()
    PHASE2 = Path(args.phase2_root)
    TECH_ROOT = Path(args.technical_root)
    JOIN_OUTPUT_ROOT = Path(args.join_output_root)
    EN_ROOT = Path(args.elastic_net_root)
    JOIN_OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    reports = []
    burden_rows = []
    for region in args.regions:
        print(f"==== {region} ====")
        report, brow = complete_region(region, args.min_reciprocal_overlap)
        reports.append(report)
        burden_rows.append(brow)
        print(
            f"  umap {report['frac_umap_before']:.3f} -> {report['frac_umap_after']:.3f} "
            f"(n={report['n_burden_with_umap_after']}/{report['n_burden']})"
        )

    write_tsv(JOIN_OUTPUT_ROOT / "tech_join_completeness.tsv", reports)

    if args.join_only:
        print("Join-only mode: Phase 2 burden tables updated; rerun Phase 2 models next")
        return

    if JOIN_OUTPUT_ROOT != TECH_ROOT:
        raise SystemExit(
            "Updating Phase 6 robustness summaries requires join-output-root == technical-root; "
            "use --join-only for immutable comparison runs"
        )

    # Update consolidated robustness table burden rows (preserve other analyses / cell columns)
    cons_path = PHASE6 / "consolidated_robustness_table.tsv"
    cons = pd.read_csv(cons_path, sep="\t")
    new_burden = pd.DataFrame(burden_rows)
    # keep non-burden rows; replace burden rows but retain cell* columns if present
    cell_cols = [c for c in cons.columns if "cell" in c.lower() or c.startswith("oligo_")]
    keep_other = cons[cons["analysis"] != "predictability_meqtl_burden_association"].copy()
    old_burden = cons[cons["analysis"] == "predictability_meqtl_burden_association"].copy()
    if cell_cols and len(old_burden):
        new_burden = new_burden.merge(
            old_burden[["region"] + cell_cols], on="region", how="left"
        )
    # align columns
    all_cols = list(dict.fromkeys(list(keep_other.columns) + list(new_burden.columns)))
    updated = pd.concat([keep_other, new_burden], ignore_index=True, sort=False)
    # stable order: original region/analysis order preference
    region_order = {"caudate": 0, "dlpfc": 1, "hippocampus": 2}
    analysis_order = {
        "LINE_L1_enrichment_vs_predictability": 0,
        "H3K9me3_enrichment_vs_predictability": 1,
        "quiescent_chromatin_enrichment_vs_predictability": 2,
        "predictability_meqtl_burden_association": 3,
    }
    updated["_r"] = updated["region"].map(region_order).fillna(9)
    updated["_a"] = updated["analysis"].map(analysis_order).fillna(9)
    updated = updated.sort_values(["_r", "_a"]).drop(columns=["_r", "_a"])
    for c in all_cols:
        if c not in updated.columns:
            updated[c] = np.nan
    updated.to_csv(cons_path, sep="\t", index=False)

    # per-region robustness burden update
    for region, brow in zip(args.regions, burden_rows):
        rp = PHASE6 / region / "robustness_results.tsv"
        if rp.exists():
            rr = pd.read_csv(rp, sep="\t")
            rr = rr[rr["analysis"] != "predictability_meqtl_burden_association"]
            rr = pd.concat([rr, pd.DataFrame([brow])], ignore_index=True, sort=False)
            rr.to_csv(rp, sep="\t", index=False)

    print(pd.DataFrame(reports)[
        ["region", "frac_umap_before", "frac_umap_after", "delta_frac_umap", "frac_snp_prox_after"]
    ].to_string(index=False))
    print(f"Wrote {PHASE6 / 'tech_join_completeness.tsv'} and updated {cons_path}")


if __name__ == "__main__":
    main()
