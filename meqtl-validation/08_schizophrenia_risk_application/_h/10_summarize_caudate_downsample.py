#!/usr/bin/env python3
"""Summarize Tier A caudate downsample vs full DLPFC/hippocampus Phase 7 results."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
P7 = PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=str(P7 / "caudate_downsample"))
    return p.parse_args()


def region_summary(region: str) -> dict:
    path = P7 / region / "risk_variant_cpg_tests_summary.tsv"
    s = pd.read_csv(path, sep="\t").iloc[0]
    return {
        "region": region,
        "n_samples": {"caudate": 153, "dlpfc": 111, "hippocampus": 116}[region],
        "n_sig_fdr": int(s["n_sig_fdr"]),
        "n_sig_loci": int(s["n_sig_loci"]),
        "n_sig_vmrs": int(s["n_sig_vmrs"]),
        "n_tests": int(s["n_tests"]),
    }


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    reps = pd.read_csv(outdir / "downsample_replicate_results.tsv", sep="\t")
    design = pd.read_csv(outdir / "downsample_design_summary.tsv", sep="\t").iloc[0]

    full = {r: region_summary(r) for r in ["caudate", "dlpfc", "hippocampus"]}
    med_pairs = float(reps["n_sig_fdr"].median())
    med_loci = float(reps["n_sig_loci"].median())
    mean_pairs = float(reps["n_sig_fdr"].mean())
    q25, q75 = np.percentile(reps["n_sig_fdr"], [25, 75])

    # Fraction of replicates with sig count <= comparator
    frac_le_dlpfc = float((reps["n_sig_fdr"] <= full["dlpfc"]["n_sig_fdr"]).mean())
    frac_le_hip = float((reps["n_sig_fdr"] <= full["hippocampus"]["n_sig_fdr"]).mean())
    frac_le_dlpfc_loci = float((reps["n_sig_loci"] <= full["dlpfc"]["n_sig_loci"]).mean())
    frac_le_hip_loci = float((reps["n_sig_loci"] <= full["hippocampus"]["n_sig_loci"]).mean())

    # Still higher than comparator?
    ratio_vs_dlpfc = med_pairs / max(full["dlpfc"]["n_sig_fdr"], 1)
    ratio_vs_hip = med_pairs / max(full["hippocampus"]["n_sig_fdr"], 1)
    loci_ratio_dlpfc = med_loci / max(full["dlpfc"]["n_sig_loci"], 1)
    loci_ratio_hip = med_loci / max(full["hippocampus"]["n_sig_loci"], 1)

    # Decision heuristics (transparent, not automatic manuscript retain)
    # Pass "not solely N" if median downsampled discoveries remain meaningfully above
    # both comparators (ratio>=1.2) AND majority of replicates exceed both.
    exceed_both = float(
        ((reps["n_sig_fdr"] > full["dlpfc"]["n_sig_fdr"])
         & (reps["n_sig_fdr"] > full["hippocampus"]["n_sig_fdr"])).mean()
    )
    not_solely_n = bool(ratio_vs_dlpfc >= 1.2 and ratio_vs_hip >= 1.2 and exceed_both >= 0.8)
    collapses_to_comparator = bool(frac_le_dlpfc >= 0.5 or frac_le_hip >= 0.5)

    comp_rows = [
        {**full["caudate"], "analysis": "full"},
        {
            "analysis": "caudate_downsample_median",
            "region": "caudate",
            "n_samples": int(design["target_n"]),
            "n_sig_fdr": med_pairs,
            "n_sig_loci": med_loci,
            "n_sig_vmrs": float(reps["n_sig_vmrs"].median()),
            "n_tests": float(reps["n_tests"].median()),
        },
        {**full["dlpfc"], "analysis": "full"},
        {**full["hippocampus"], "analysis": "full"},
    ]
    write_tsv(outdir / "downsample_vs_regions.tsv", comp_rows)

    claim = [{
        "n_reps": int(len(reps)),
        "target_n": int(design["target_n"]),
        "n_shared3": int(design["n_shared_all3"]),
        "full_caudate_n_sig_pairs": full["caudate"]["n_sig_fdr"],
        "full_caudate_n_sig_loci": full["caudate"]["n_sig_loci"],
        "downsample_median_n_sig_pairs": med_pairs,
        "downsample_mean_n_sig_pairs": mean_pairs,
        "downsample_iqr_n_sig_pairs": f"{q25:.1f}-{q75:.1f}",
        "downsample_median_n_sig_loci": med_loci,
        "dlpfc_n_sig_pairs": full["dlpfc"]["n_sig_fdr"],
        "hippocampus_n_sig_pairs": full["hippocampus"]["n_sig_fdr"],
        "median_pair_ratio_vs_dlpfc": ratio_vs_dlpfc,
        "median_pair_ratio_vs_hippocampus": ratio_vs_hip,
        "median_loci_ratio_vs_dlpfc": loci_ratio_dlpfc,
        "median_loci_ratio_vs_hippocampus": loci_ratio_hip,
        "frac_reps_sig_pairs_le_dlpfc": frac_le_dlpfc,
        "frac_reps_sig_pairs_le_hippocampus": frac_le_hip,
        "frac_reps_sig_loci_le_dlpfc": frac_le_dlpfc_loci,
        "frac_reps_sig_loci_le_hippocampus": frac_le_hip_loci,
        "frac_reps_exceed_both_comparators": exceed_both,
        "criterion_not_solely_sample_size": not_solely_n,
        "flag_collapses_toward_comparator": collapses_to_comparator,
        "interpretation": (
            "Caudate SCZ risk-meQTL discovery remains higher than DLPFC/hippocampus after "
            "N-matched downsampling."
            if not_solely_n else
            "After N-matched downsampling, caudate discovery advantage is attenuated and may "
            "be partly/mostly sample-size driven; do not claim caudate selectivity from "
            "full-N significance alone."
        ),
    }]
    write_tsv(outdir / "downsample_claim_snapshot.tsv", claim)

    # Prioritized loci stability
    prior_path = outdir / "downsample_prioritized_locus_results.tsv.gz"
    if prior_path.exists():
        pr = pd.read_csv(prior_path, sep="\t")
        stab = (
            pr.groupby(["index_snp", "locus_id"], as_index=False)
            .agg(
                n_reps=("replicate", "nunique"),
                frac_reps_any_sig=("any_sig", "mean"),
                median_n_sig=("n_sig_fdr", "median"),
                median_min_qval=("min_qval", "median"),
                median_best_beta=("best_beta", "median"),
                median_best_abs_beta=("best_beta", lambda s: float(np.median(np.abs(s)))),
            )
            .sort_values("frac_reps_any_sig", ascending=False)
        )
        stab.to_csv(outdir / "downsample_prioritized_stability.tsv", sep="\t", index=False)

    print("Downsample claim:")
    print(pd.DataFrame(claim).T.to_string())
    print(f"Wrote summaries under {outdir}")


if __name__ == "__main__":
    main()
