#!/usr/bin/env python3
"""Summarize Phase 4 caudate N-matched downsample vs DLPFC/hippocampus discovery."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
CROSS = PROJECT / "meqtl-validation/04_cross_region_sharing/_m"
PHASE1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping"
PHASE2 = PROJECT / "meqtl-validation/02_vmr_meqtl_burden/_m"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=str(CROSS / "caudate_downsample"))
    return p.parse_args()


def region_discovery(region: str) -> dict:
    qc = pd.read_csv(PHASE1 / region / "_m/tensorqtl/qc/meqtl_qc_summary.tsv", sep="\t").iloc[0]
    burden = PHASE2 / region / "aggregation_summary.tsv"
    n_vmr_sup = np.nan
    if burden.exists():
        b = pd.read_csv(burden, sep="\t")
        if "n_vmrs_with_any_sig_meqtl" in b.columns:
            n_vmr_sup = float(b.iloc[0]["n_vmrs_with_any_sig_meqtl"])
        elif "n_vmrs_meqtl_supported" in b.columns:
            n_vmr_sup = float(b.iloc[0]["n_vmrs_meqtl_supported"])
    n_samp = {"caudate": 153, "dlpfc": 111, "hippocampus": 116}[region]
    return {
        "region": region,
        "n_samples": n_samp,
        "n_sig_fdr": int(qc["n_significant_fdr"]),
        "n_tested": int(qc["n_phenotypes_tested"]),
        "n_sig_vmrs": n_vmr_sup,
        "frac_sig": float(qc["n_significant_fdr"]) / float(qc["n_phenotypes_tested"]),
    }


def update_pending(status: str, reason: str) -> None:
    path = CROSS / "pending_analyses.tsv"
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t")
    mask = df["analysis"] == "caudate_donor_downsample_remap"
    if mask.any():
        df.loc[mask, "status"] = status
        df.loc[mask, "reason"] = reason
        df.to_csv(path, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    reps = pd.read_csv(outdir / "downsample_replicate_results.tsv", sep="\t")
    design = pd.read_csv(outdir / "downsample_design_summary.tsv", sep="\t").iloc[0]

    full = {r: region_discovery(r) for r in ["caudate", "dlpfc", "hippocampus"]}
    med_sig = float(reps["n_sig_fdr"].median())
    mean_sig = float(reps["n_sig_fdr"].mean())
    q25, q75 = np.percentile(reps["n_sig_fdr"], [25, 75])
    med_vmr = float(reps["n_sig_vmrs"].median()) if "n_sig_vmrs" in reps.columns else np.nan
    med_retain = (
        float(reps["frac_full_sig_retained"].median())
        if "frac_full_sig_retained" in reps.columns else np.nan
    )

    frac_le_dlpfc = float((reps["n_sig_fdr"] <= full["dlpfc"]["n_sig_fdr"]).mean())
    frac_le_hip = float((reps["n_sig_fdr"] <= full["hippocampus"]["n_sig_fdr"]).mean())
    ratio_vs_dlpfc = med_sig / max(full["dlpfc"]["n_sig_fdr"], 1)
    ratio_vs_hip = med_sig / max(full["hippocampus"]["n_sig_fdr"], 1)
    exceed_both = float(
        ((reps["n_sig_fdr"] > full["dlpfc"]["n_sig_fdr"])
         & (reps["n_sig_fdr"] > full["hippocampus"]["n_sig_fdr"])).mean()
    )
    # Architecture criterion: after N-match, caudate still exceeds both comparators
    not_solely_n = bool(ratio_vs_dlpfc >= 1.1 and ratio_vs_hip >= 1.1 and exceed_both >= 0.7)
    collapses = bool(frac_le_dlpfc >= 0.5 or frac_le_hip >= 0.5)

    comp_rows = [
        {**full["caudate"], "analysis": "full"},
        {
            "analysis": "caudate_downsample_median",
            "region": "caudate",
            "n_samples": int(design["target_n"]),
            "n_sig_fdr": med_sig,
            "n_tested": float(reps["n_tests"].median()),
            "n_sig_vmrs": med_vmr,
            "frac_sig": float(reps["frac_cpgs_sig"].median()) if "frac_cpgs_sig" in reps else np.nan,
        },
        {**full["dlpfc"], "analysis": "full"},
        {**full["hippocampus"], "analysis": "full"},
    ]
    write_tsv(outdir / "downsample_vs_regions.tsv", comp_rows)

    claim = [{
        "n_reps": int(len(reps)),
        "target_n": int(design["target_n"]),
        "n_shared3": int(design["n_shared_all3"]),
        "method": "lead_snp_retention_M3a_residualized",
        "full_caudate_n_sig": full["caudate"]["n_sig_fdr"],
        "downsample_median_n_sig": med_sig,
        "downsample_mean_n_sig": mean_sig,
        "downsample_iqr_n_sig": f"{q25:.1f}-{q75:.1f}",
        "downsample_median_n_sig_vmrs": med_vmr,
        "median_frac_full_sig_retained": med_retain,
        "dlpfc_n_sig": full["dlpfc"]["n_sig_fdr"],
        "hippocampus_n_sig": full["hippocampus"]["n_sig_fdr"],
        "median_ratio_vs_dlpfc": ratio_vs_dlpfc,
        "median_ratio_vs_hippocampus": ratio_vs_hip,
        "frac_reps_sig_le_dlpfc": frac_le_dlpfc,
        "frac_reps_sig_le_hippocampus": frac_le_hip,
        "frac_reps_exceed_both_comparators": exceed_both,
        "criterion_not_solely_sample_size": not_solely_n,
        "flag_collapses_toward_comparator": collapses,
        "interpretation": (
            "Caudate cis-meQTL discovery remains higher than DLPFC/hippocampus after "
            "N-matched lead-SNP retention downsampling; excess is not explained solely by N."
            if not_solely_n else
            "After N-matched downsampling, caudate discovery advantage attenuates toward "
            "DLPFC/hippocampus; do not claim caudate-selective discovery from full-N counts alone."
        ),
    }]
    write_tsv(outdir / "downsample_claim_snapshot.tsv", claim)

    update_pending(
        status="done",
        reason=(
            f"Lead-SNP retention downsample complete (n_reps={len(reps)}, target_n={int(design['target_n'])}); "
            f"median_n_sig={med_sig:.0f}; not_solely_N={not_solely_n}"
        ),
    )
    print("Downsample claim:")
    print(pd.DataFrame(claim).T.to_string())
    print(f"Wrote summaries under {outdir}")


if __name__ == "__main__":
    main()
