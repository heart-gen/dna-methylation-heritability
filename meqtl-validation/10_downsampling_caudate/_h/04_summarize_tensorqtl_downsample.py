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
                "discovery_rate": int(s["n_significant_fdr"]) / int(s["n_phenotypes_tested"]),
                "lambda_gc": float(s["lambda_gc"]),
                "cis_qtl_path": str(s["cis_qtl_path"]),
                "qc_path": str(path),
            }
    raise FileNotFoundError(f"No QC summary for {region}")


def load_qvalues(path: Path) -> pd.Series:
    """Read a phenotype-indexed q-value vector from a TensorQTL cis table."""
    if not path.exists():
        raise FileNotFoundError(path)
    table = pd.read_csv(path, sep="\t", index_col=0, usecols=lambda c: c == "qval" or c.startswith("phenotype"))
    if "qval" not in table.columns:
        raise ValueError(f"{path} lacks qval")
    qval = pd.to_numeric(table["qval"], errors="coerce")
    qval.index = qval.index.astype(str)
    return qval


def retention_vs_full(full: pd.Series, ds: pd.Series, fdr: float) -> dict:
    """Fraction of full-N FDR-significant CpGs that remain FDR-significant in downsample."""
    full_sig = set(full.index[full <= fdr])
    ds_sig = set(ds.index[ds <= fdr])
    retained = full_sig & ds_sig
    return {
        "n_full_sig": len(full_sig),
        "n_full_sig_retained": len(retained),
        "frac_full_sig_retained": float(len(retained) / len(full_sig)) if full_sig else np.nan,
    }


def common_universe_rates(
    downsample: pd.Series,
    dlpfc: pd.Series,
    hippocampus: pd.Series,
    fdr: float,
) -> dict:
    """Discovery rates on the identical CpG universe testable in all datasets."""
    common = downsample.index.intersection(dlpfc.index).intersection(hippocampus.index)
    if not len(common):
        return {"n_common_cpgs_all3": 0}
    ds_rate = float((downsample.loc[common] <= fdr).mean())
    dlpfc_rate = float((dlpfc.loc[common] <= fdr).mean())
    hip_rate = float((hippocampus.loc[common] <= fdr).mean())
    return {
        "n_common_cpgs_all3": int(len(common)),
        "caudate_downsample_common_discovery_rate": ds_rate,
        "dlpfc_common_discovery_rate": dlpfc_rate,
        "hippocampus_common_discovery_rate": hip_rate,
        "common_rate_ratio_vs_dlpfc": ds_rate / max(dlpfc_rate, 1e-12),
        "common_rate_ratio_vs_hippocampus": ds_rate / max(hip_rate, 1e-12),
        "common_rate_exceeds_both": bool(ds_rate > dlpfc_rate and ds_rate > hip_rate),
    }


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    manifest = pd.read_csv(outdir / "downsample_replicate_manifest.tsv", sep="\t")
    design = pd.read_csv(outdir / "downsample_design_summary.tsv", sep="\t").iloc[0]

    full = {r: region_qc(r) for r in ["caudate", "dlpfc", "hippocampus"]}
    # Resolve full cis tables once for retention and common-universe comparisons.
    full_cis = Path(full["caudate"]["cis_qtl_path"])
    if not full_cis.exists():
        alt = PHASE1 / "caudate/_m/tensorqtl/cpg_meqtl_caudate.cis_qtl.txt.gz"
        full_cis = alt if alt.exists() else full_cis
    full_q = load_qvalues(full_cis)
    comparator_q = {}
    for region in ["dlpfc", "hippocampus"]:
        path = Path(full[region]["cis_qtl_path"])
        if not path.exists():
            path = PHASE1 / region / "_m/tensorqtl" / f"cpg_meqtl_{region}.cis_qtl.txt.gz"
        comparator_q[region] = load_qvalues(path)

    rep_rows = []
    for _, mrow in manifest.iterrows():
        rep = int(mrow["replicate"])
        qc_path = Path(mrow["outdir"]) / "qc" / "meqtl_qc_summary.tsv"
        cis_path = Path(mrow["outdir"]) / f"{mrow['prefix']}.cis_qtl.txt.gz"
        if not qc_path.exists():
            print(f"WARNING: missing QC for rep {rep}: {qc_path}")
            continue
        qc = pd.read_csv(qc_path, sep="\t").iloc[0]
        ds_q = load_qvalues(cis_path)
        ret = retention_vs_full(full_q, ds_q, args.fdr)
        common = common_universe_rates(
            ds_q, comparator_q["dlpfc"], comparator_q["hippocampus"], args.fdr
        )
        rep_rows.append({
            "replicate": rep,
            "n_samples": int(mrow["n_samples"]),
            "n_tested": int(qc["n_phenotypes_tested"]),
            "n_sig_fdr": int(qc["n_significant_fdr"]),
            "discovery_rate": int(qc["n_significant_fdr"]) / int(qc["n_phenotypes_tested"]),
            "lambda_gc": float(qc["lambda_gc"]),
            "median_pval": float(qc["median_pval"]),
            "cis_qtl_path": str(cis_path),
            **ret,
            **common,
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
    med_rate = float(reps["discovery_rate"].median())

    frac_le_dlpfc = float((reps["n_sig_fdr"] <= full["dlpfc"]["n_sig_fdr"]).mean())
    frac_le_hip = float((reps["n_sig_fdr"] <= full["hippocampus"]["n_sig_fdr"]).mean())
    # Primary regional comparison uses the same CpGs in all three datasets.
    ratio_dlpfc = float(reps["common_rate_ratio_vs_dlpfc"].median())
    ratio_hip = float(reps["common_rate_ratio_vs_hippocampus"].median())
    frac_common_exceed_both = float(reps["common_rate_exceeds_both"].mean())
    exceed_both = float(
        ((reps["n_sig_fdr"] > full["dlpfc"]["n_sig_fdr"])
         & (reps["n_sig_fdr"] > full["hippocampus"]["n_sig_fdr"])).mean()
    )
    # Absolute discoveries and rates on different region-specific VMR universes
    # are not comparable. Require an identical-CpG rate advantage that is stable
    # across at least 80% of donor-downsampling replicates.
    not_solely_n = bool(
        ratio_dlpfc >= 1.1 and ratio_hip >= 1.1 and frac_common_exceed_both >= 0.8
    )
    collapses = bool(frac_common_exceed_both < 0.5)

    comp = [
        {**full["caudate"], "analysis": "full_tensorqtl_M3a"},
        {
            "analysis": "caudate_downsample_tensorqtl_median",
            "region": "caudate",
            "n_samples": int(design["target_n"]),
            "n_tested": float(reps["n_tested"].median()),
            "n_sig_fdr": med_sig,
            "lambda_gc": med_lambda,
            "discovery_rate": med_rate,
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
        "downsample_median_discovery_rate": med_rate,
        "downsample_median_n_common_cpgs_all3": float(reps["n_common_cpgs_all3"].median()),
        "downsample_median_common_discovery_rate": float(
            reps["caudate_downsample_common_discovery_rate"].median()
        ),
        "dlpfc_median_common_discovery_rate": float(reps["dlpfc_common_discovery_rate"].median()),
        "hippocampus_median_common_discovery_rate": float(
            reps["hippocampus_common_discovery_rate"].median()
        ),
        "median_frac_full_sig_retained": med_retain,
        "dlpfc_n_sig": full["dlpfc"]["n_sig_fdr"],
        "hippocampus_n_sig": full["hippocampus"]["n_sig_fdr"],
        "dlpfc_discovery_rate": full["dlpfc"]["discovery_rate"],
        "hippocampus_discovery_rate": full["hippocampus"]["discovery_rate"],
        "median_discovery_rate_ratio_vs_dlpfc": ratio_dlpfc,
        "median_discovery_rate_ratio_vs_hippocampus": ratio_hip,
        "frac_reps_common_rate_exceeds_both": frac_common_exceed_both,
        "frac_reps_sig_le_dlpfc": frac_le_dlpfc,
        "frac_reps_sig_le_hippocampus": frac_le_hip,
        "frac_reps_exceed_both_comparators": exceed_both,
        "criterion_not_solely_sample_size": not_solely_n,
        "flag_collapses_toward_comparator": collapses,
        "fdr_family_matched": True,
        "interpretation": (
            "After official TensorQTL permutation-FDR remapping at N=111, caudate's "
            "discovery rate on the identical all-region CpG universe remains at least "
            "10% higher than both comparison regions and exceeds both in at least 80% "
            "of replicates."
            if not_solely_n else
            "The N-matched caudate discovery rate on the identical all-region CpG universe "
            "does not satisfy the prespecified magnitude and stability gate; do not claim "
            "caudate-selective discovery from unequal-universe counts or rates."
        ),
    }]
    write_tsv(outdir / "tensorqtl_downsample_claim_snapshot.tsv", claim)
    print(pd.DataFrame(claim).T.to_string())
    print(f"Wrote summaries under {outdir}")


if __name__ == "__main__":
    main()
