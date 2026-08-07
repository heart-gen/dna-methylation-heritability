#!/usr/bin/env python3
"""Aggregate CpG-level cis-meQTL results to VMR-level burden metrics."""

from __future__ import annotations

import argparse
import glob as _glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import load_paths, load_yaml, write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True)
    p.add_argument("--cis-qtl", required=True, help="lead SNP table or cis_qtl.txt.gz")
    p.add_argument("--cpg-vmr-map-glob", default="", help="Glob for cpg_vmr_map.chr*.tsv")
    p.add_argument("--population", default="AA", choices=["AA", "EA"])
    p.add_argument("--outdir", default="", help="Optional output directory override")
    p.add_argument("--fdr", type=float, default=0.05)
    return p.parse_args()


def _aggregate_burden(merged: pd.DataFrame, region: str) -> pd.DataFrame:
    """Vectorized VMR-level meQTL burden metrics."""
    m = merged.copy()
    m["vmr_id"] = m["vmr_id"].astype(str)
    if "sig_meqtl" not in m.columns:
        m["sig_meqtl"] = False
    m["sig_meqtl"] = m["sig_meqtl"].astype("boolean").fillna(False).astype(bool)

    pcol = "pval_beta" if "pval_beta" in m.columns else (
        "pval_nominal" if "pval_nominal" in m.columns else None
    )
    if pcol:
        m["_pval"] = pd.to_numeric(m[pcol], errors="coerce")
    else:
        m["_pval"] = np.nan
    m["_slope_abs"] = pd.to_numeric(m["slope"], errors="coerce").abs() if "slope" in m.columns else np.nan
    if "start_distance" in m.columns:
        m["_dist"] = pd.to_numeric(m["start_distance"], errors="coerce").abs()
    elif "tss_distance" in m.columns:
        m["_dist"] = pd.to_numeric(m["tss_distance"], errors="coerce").abs()
    else:
        m["_dist"] = np.nan
    m["_cov"] = pd.to_numeric(m["mean_coverage"], errors="coerce") if "mean_coverage" in m.columns else np.nan
    m["_var"] = pd.to_numeric(m["variance"], errors="coerce") if "variance" in m.columns else np.nan
    m["_frac"] = pd.to_numeric(m["fraction_nonmissing"], errors="coerce") if "fraction_nonmissing" in m.columns else np.nan

    agg_spec: dict = {
        "n_tested_cpgs": ("phenotype_id", "size"),
        "n_cpgs_with_sig_meqtl": ("sig_meqtl", "sum"),
        "min_cpg_pvalue": ("_pval", "min"),
        "strongest_abs_beta": ("_slope_abs", "max"),
        "median_cpg_lead_snp_distance": ("_dist", "median"),
        "average_cpg_nonmissing_fraction": ("_frac", "mean"),
        "average_cpg_coverage": ("_cov", "mean"),
        "mean_cpg_variance": ("_var", "mean"),
    }
    if "variant_id" in m.columns:
        agg_spec["n_distinct_lead_snps"] = ("variant_id", "nunique")
    agg = m.groupby("vmr_id", sort=False).agg(**agg_spec).reset_index()
    if "n_distinct_lead_snps" not in agg.columns:
        agg["n_distinct_lead_snps"] = 0
    agg["n_cpgs_with_sig_meqtl"] = agg["n_cpgs_with_sig_meqtl"].astype(int)
    agg["proportion_cpgs_with_sig_meqtl"] = (
        agg["n_cpgs_with_sig_meqtl"] / agg["n_tested_cpgs"].replace(0, np.nan)
    )
    agg["region"] = region
    return agg


def main() -> None:
    args = parse_args()
    paths = load_paths()
    meqtl = load_yaml("meqtl_parameters.yml")
    project = Path(paths["project_root"])
    region = args.region
    if args.population == "AA":
        outdir = project / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m" / region
        prepared = (
            project / "meqtl-validation" / "01_cpg_meqtl_mapping" / region / "_m" / "prepared"
        )
        pred_key = "local_predictability_summary_template"
    else:
        outdir = (
            project / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m" / args.population / region
        )
        prepared = (
            project / "meqtl-validation" / "01_cpg_meqtl_mapping" / region / "_m"
            / "prepared" / args.population
        )
        pred_key = "local_predictability_ea_template"
    if args.outdir:
        outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    lead = pd.read_csv(args.cis_qtl, sep="\t", compression="infer")
    if "phenotype_id" not in lead.columns:
        lead = lead.rename(columns={lead.columns[0]: "phenotype_id"})
    qcol = "qval" if "qval" in lead.columns else None
    if qcol:
        lead["sig_meqtl"] = lead[qcol].le(args.fdr).fillna(False).astype(bool)
    else:
        lead["sig_meqtl"] = False

    if args.cpg_vmr_map_glob:
        maps = sorted(Path(p) for p in _glob.glob(args.cpg_vmr_map_glob))
    else:
        maps = sorted(prepared.glob("cpg_vmr_map.chr*.tsv"))
    if not maps:
        raise SystemExit(f"No cpg_vmr_map files found under {prepared}; run Phase 1 preparation first")
    cpg_map = pd.concat([pd.read_csv(p, sep="\t") for p in maps], ignore_index=True)

    merged = cpg_map.merge(lead, on="phenotype_id", how="left")
    burden = _aggregate_burden(merged, region)

    pred_path = project / paths[pred_key].format(region=region)
    pred = pd.read_csv(pred_path, sep="\t")
    score_col = meqtl["predictability_score_column"]
    if score_col not in pred.columns:
        raise SystemExit(f"Missing {score_col} in {pred_path}")
    pred = pred.copy()
    pred["vmr_coord_id"] = (
        pred["chrom"].astype(str).str.replace("chr", "", regex=False)
        + ":"
        + pred["start"].astype(str)
        + "-"
        + pred["end"].astype(str)
    )
    pred["task_id"] = pred["task_id"].astype(str)

    burden = burden.merge(
        pred[["task_id", score_col]].rename(
            columns={"task_id": "vmr_id", score_col: "local_predictability"}
        ),
        on="vmr_id",
        how="left",
    )
    missing = burden["local_predictability"].isna()
    if missing.any():
        extra = burden.loc[missing, ["vmr_id"]].merge(
            pred[["vmr_coord_id", score_col]].rename(
                columns={"vmr_coord_id": "vmr_id", score_col: "local_predictability"}
            ),
            on="vmr_id",
            how="left",
        )
        burden.loc[missing, "local_predictability"] = extra["local_predictability"].to_numpy()

    cov_path = prepared / "vmr_mean_coverage.tsv"
    if cov_path.exists():
        vc = pd.read_csv(cov_path, sep="\t")
        vc["vmr_id"] = vc["vmr_id"].astype(str)
        keep = [c for c in ["vmr_id", "mean_cpg_coverage", "mean_frac_samples_cov_ge_min"] if c in vc.columns]
        burden = burden.merge(vc[keep], on="vmr_id", how="left", suffixes=("", "_vmr"))
        if "mean_cpg_coverage" in burden.columns:
            burden["average_cpg_coverage"] = burden["average_cpg_coverage"].fillna(
                burden["mean_cpg_coverage"]
            )

    tech_path = (
        project / "meqtl-validation" / "07_repeat_mappability_sensitivity" / "_m"
        / region / "vmr_technical_annotations.tsv"
    )
    if tech_path.exists():
        t = pd.read_csv(tech_path, sep="\t")
        feat = [
            c for c in [
                "length", "blacklist_frac", "segdup_frac", "line_l1_frac",
                "umap_k24_mean", "high_mappability", "overlaps_blacklist", "overlaps_segdup",
            ] if c in t.columns
        ]
        if feat and "task_id" in t.columns:
            by_task = t.dropna(subset=["task_id"]).copy()
            by_task["vmr_id"] = by_task["task_id"].astype(float).astype(int).astype(str)
            by_task = by_task.drop_duplicates("vmr_id", keep="first")
            burden = burden.merge(by_task[["vmr_id"] + feat], on="vmr_id", how="left")
        if feat and "interval_id" in t.columns:
            by_int = t.copy()
            by_int["vmr_id"] = by_int["interval_id"].astype(str)
            by_int = by_int.drop_duplicates("vmr_id", keep="first")
            for c in feat:
                if c not in burden.columns:
                    burden[c] = np.nan
            miss = burden[feat[0]].isna()
            if miss.any():
                extra = burden.loc[miss, ["vmr_id"]].merge(
                    by_int[["vmr_id"] + feat], on="vmr_id", how="left"
                )
                for c in feat:
                    burden.loc[miss, c] = extra[c].to_numpy()

    ann_path = project / paths["vmr_annotation_template"].format(region=region)
    if Path(ann_path).exists():
        ann = pd.read_csv(ann_path, sep="\t")
        if "vmr_id" in ann.columns or "task_id" in ann.columns:
            key = "vmr_id" if "vmr_id" in ann.columns else "task_id"
            keep_ann = [
                c for c in ann.columns
                if c == key or c.lower() in {
                    "mappability", "repeat_overlap", "genomic_annotation",
                    "chromatin_annotation", "line_overlap", "l1_overlap",
                }
            ]
            to_merge = ann[keep_ann].rename(columns={key: "vmr_id"}).copy()
            to_merge["vmr_id"] = to_merge["vmr_id"].astype(str)
            to_merge = to_merge.drop_duplicates("vmr_id", keep="first")
            burden = burden.merge(to_merge, on="vmr_id", how="left", suffixes=("", "_ann"))
        elif {"seqnames", "start", "end"}.issubset(ann.columns):
            # Coordinate-keyed annotation table (no task_id)
            ann = ann.copy()
            ann["vmr_id"] = (
                ann["seqnames"].astype(str).str.replace("^chr", "", regex=True)
                + ":"
                + ann["start"].astype(str)
                + "-"
                + ann["end"].astype(str)
            )
            keep = ["vmr_id"]
            for c in ["annot.type", "h2_category", "annot.symbol"]:
                if c in ann.columns:
                    keep.append(c)
            to_merge = ann[keep].drop_duplicates("vmr_id", keep="first")
            burden = burden.merge(to_merge, on="vmr_id", how="left", suffixes=("", "_ann"))

    out = outdir / "vmr_meqtl_burden.tsv.gz"
    burden.to_csv(out, sep="\t", index=False, compression="gzip")
    write_tsv(
        outdir / "aggregation_summary.tsv",
        [{
            "region": region,
            "n_vmrs": len(burden),
            "n_vmrs_with_predictability": int(burden["local_predictability"].notna().sum()),
            "n_vmrs_with_tech": int(burden["umap_k24_mean"].notna().sum()) if "umap_k24_mean" in burden.columns else 0,
            "mean_proportion_sig": float(burden["proportion_cpgs_with_sig_meqtl"].mean()),
            "output": str(out),
        }],
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
