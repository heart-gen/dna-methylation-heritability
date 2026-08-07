#!/usr/bin/env python3
"""Summarize official TensorQTL caudate downsample vs DLPFC/hippocampus M3a."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
OUTDIR = PROJECT / "meqtl-validation/10_downsampling_caudate/_m"
PHASE1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping"
FDR = 0.05


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=str(OUTDIR))
    p.add_argument("--fdr", type=float, default=FDR)
    return p.parse_args()


def region_qc(region: str) -> dict:
    # Prefer locked primary qc; fall back to M3a sensitivity path
    paths = [
        PHASE1 / region / "_m/tensorqtl/qc/meqtl_qc_summary.tsv",
        PHASE1 / region / "_m/covariate_sensitivity/tensorqtl/M3a/qc/meqtl_qc_summary.tsv",
    ]
    for path in paths:
        if path.exists():
            s = pd.read_csv(path, sep="\t").iloc[0]
            return {
                "region": region,
                "n_samples": {"caudate": 153, "dlpfc": 111, "hippocampus": 116}[region],
                "n_tested": int(s["n_phenotypes_tested"]),
                "n_sig_fdr": int(s["n_significant_fdr"]),
                "lambda_gc": float(s["lambda_gc"]),
                "cis_qtl_path": str(s["cis_qtl_path"]),
                "qc_path": str(path),
            }
    raise FileNotFoundError(f"No QC summary for {region}")


def retention_vs_full(full_cis: Path, ds_cis: Path, fdr: float) -> dict:
    """Fraction of full-N FDR-significant CpGs that remain FDR-significant in downsample."""
    if not full_cis.exists() or not ds_cis.exists():
        return {
            "n_full_sig": np.nan,
            "n_full_sig_retained": np.nan,
            "frac_full_sig_retained": np.nan,
        }
    full = pd.read_csv(full_cis, sep="\t", index_col=0, usecols=lambda c: True)
    # Only need qval; re-read lightly if huge
    if "qval" not in full.columns:
        return {
            "n_full_sig": np.nan,
            "n_full_sig_retained": np.nan,
            "frac_full_sig_retained": np.nan,
        }
    ds = pd.read_csv(ds_cis, sep="\t", index_col=0)
    if "qval" not in ds.columns:
        return {
            "n_full_sig": np.nan,
            "n_full_sig_retained": np.nan,
            "frac_full_sig_retained": np.nan,
        }
    full_sig = set(full.index[pd.to_numeric(full["qval"], errors="coerce") <= fdr].astype(str))
    ds_sig = set(ds.index[pd.to_numeric(ds["qval"], errors="coerce") <= fdr].astype(str))
    retained = full_sig & ds_sig
    return {
        "n_full_sig": len(full_sig),
        "n_full_sig_retained": len(retained),
        "frac_full_sig_retained": float(len(retained) / len(full_sig)) if full_sig else np.nan,
    }


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    manifest = pd.read_csv(outdir / "downsample_replicate_manifest.tsv", sep="\t")
    design = pd.read_csv(outdir / "downsample_design_summary.tsv", sep="\t").iloc[0]

    full = {r: region_qc(r) for r in ["caudate", "dlpfc", "hippocampus"]}
    # Resolve full caudate cis path for retention
    full_cis = Path(full["caudate"]["cis_qtl_path"])
    if not full_cis.exists():
        alt = PHASE1 / "caudate/_m/tensorqtl/cpg_meqtl_caudate.cis_qtl.txt.gz"
        full_cis = alt if alt.exists() else full_cis

    rep_rows = []
    for _, mrow in manifest.iterrows():
        rep = int(mrow["replicate"])
        qc_path = Path(mrow["outdir"]) / "qc" / "meqtl_qc_summary.tsv"
        cis_path = Path(mrow["outdir"]) / f"{mrow['prefix']}.cis_qtl.txt.gz"
        if not qc_path.exists():
            print(f"WARNING: missing QC for rep {rep}: {qc_path}")
            continue
        qc = pd.read_csv(qc_path, sep="\t").iloc[0]
        ret = retention_vs_full(full_cis, cis_path, args.fdr)
        rep_rows.append({
            "replicate": rep,
            "n_samples": int(mrow["n_samples"]),
            "n_tested": int(qc["n_phenotypes_tested"]),
            "n_sig_fdr": int(qc["n_significant_fdr"]),
            "lambda_gc": float(qc["lambda_gc"]),
            "median_pval": float(qc["median_pval"]),
            "cis_qtl_path": str(cis_path),
            **ret,
            "method": "tensorqtl_cis_permutation_fdr",
        })

    if not rep_rows:
        raise SystemExit("No replicate QC results found; was step_2 completed?")

    reps = pd.DataFrame(rep_rows)
    reps.to_csv(outdir / "tensorqtl_downsample_replicate_results.tsv", sep="\t", index=False)

    med_sig = float(reps["n_sig_fdr"].median())
    mean_sig = float(reps["n_sig_fdr"].mean())
    q25, q75 = np.percentile(reps["n_sig_fdr"], [25, 75])
    med_retain = float(reps["frac_full_sig_retained"].median())
    med_lambda = float(reps["lambda_gc"].median())

    frac_le_dlpfc = float((reps["n_sig_fdr"] <= full["dlpfc"]["n_sig_fdr"]).mean())
    frac_le_hip = float((reps["n_sig_fdr"] <= full["hippocampus"]["n_sig_fdr"]).mean())
    ratio_dlpfc = med_sig / max(full["dlpfc"]["n_sig_fdr"], 1)
    ratio_hip = med_sig / max(full["hippocampus"]["n_sig_fdr"], 1)
    exceed_both = float(
        ((reps["n_sig_fdr"] > full["dlpfc"]["n_sig_fdr"])
         & (reps["n_sig_fdr"] > full["hippocampus"]["n_sig_fdr"])).mean()
    )
    # Same heuristic as Phase 4 lead-retention, but now method-matched
    not_solely_n = bool(ratio_dlpfc >= 1.1 and ratio_hip >= 1.1 and exceed_both >= 0.7)
    collapses = bool(frac_le_dlpfc >= 0.5 or frac_le_hip >= 0.5)

    comp = [
        {**full["caudate"], "analysis": "full_tensorqtl_M3a"},
        {
            "analysis": "caudate_downsample_tensorqtl_median",
            "region": "caudate",
            "n_samples": int(design["target_n"]),
            "n_tested": float(reps["n_tested"].median()),
            "n_sig_fdr": med_sig,
            "lambda_gc": med_lambda,
            "cis_qtl_path": "",
            "qc_path": "",
        },
        {**full["dlpfc"], "analysis": "full_tensorqtl_M3a"},
        {**full["hippocampus"], "analysis": "full_tensorqtl_M3a"},
    ]
    write_tsv(outdir / "tensorqtl_downsample_vs_regions.tsv", comp)

    claim = [{
        "n_reps_completed": int(len(reps)),
        "n_reps_planned": int(len(manifest)),
        "target_n": int(design["target_n"]),
        "n_shared3": int(design["n_shared_all3"]),
        "method": "tensorqtl_cis_permutation_fdr_M3a",
        "full_caudate_n_sig": full["caudate"]["n_sig_fdr"],
        "downsample_median_n_sig": med_sig,
        "downsample_mean_n_sig": mean_sig,
        "downsample_iqr_n_sig": f"{q25:.1f}-{q75:.1f}",
        "downsample_median_lambda_gc": med_lambda,
        "median_frac_full_sig_retained": med_retain,
        "dlpfc_n_sig": full["dlpfc"]["n_sig_fdr"],
        "hippocampus_n_sig": full["hippocampus"]["n_sig_fdr"],
        "median_ratio_vs_dlpfc": ratio_dlpfc,
        "median_ratio_vs_hippocampus": ratio_hip,
        "frac_reps_sig_le_dlpfc": frac_le_dlpfc,
        "frac_reps_sig_le_hippocampus": frac_le_hip,
        "frac_reps_exceed_both_comparators": exceed_both,
        "criterion_not_solely_sample_size": not_solely_n,
        "flag_collapses_toward_comparator": collapses,
        "fdr_family_matched": True,
        "interpretation": (
            "After official TensorQTL permutation-FDR remapping at N=111, caudate discovery "
            "remains higher than DLPFC/hippocampus; excess is not explained solely by sample size."
            if not_solely_n else
            "After official TensorQTL permutation-FDR remapping at N=111, caudate discovery "
            "advantage attenuates toward DLPFC/hippocampus; do not claim caudate-selective "
            "discovery from full-N counts alone."
        ),
    }]
    write_tsv(outdir / "tensorqtl_downsample_claim_snapshot.tsv", claim)
    print(pd.DataFrame(claim).T.to_string())
    print(f"Wrote summaries under {outdir}")


if __name__ == "__main__":
    main()
