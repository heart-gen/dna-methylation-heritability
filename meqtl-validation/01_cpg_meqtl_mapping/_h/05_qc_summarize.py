#!/usr/bin/env python3
"""Summarize cis-meQTL permutation results and write QC tables."""

from __future__ import annotations

import argparse, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True)
    p.add_argument("--cis-qtl", required=True, help="*.cis_qtl.txt.gz from 04_tensorqtl_map.py")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--outdir", default="")
    return p.parse_args()


def genomic_inflation(pvals: np.ndarray) -> float:
    from scipy.stats import chi2

    p = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    if p.size == 0:
        return float("nan")
    chisq = chi2.ppf(1 - p, 1)
    return float(np.median(chisq) / chi2.ppf(0.5, 1))


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir) if args.outdir else Path(args.cis_qtl).parent / "qc"
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.cis_qtl, sep="\t", index_col=0)
    pcol = "pval_beta" if "pval_beta" in df.columns else "pval_nominal"
    if pcol not in df.columns:
        pcol = [c for c in df.columns if c.startswith("pval")][0]

    lam = genomic_inflation(df[pcol].to_numpy(dtype=float))
    n_sig = int((df["qval"] <= args.fdr).sum()) if "qval" in df.columns else 0
    summary = [{
        "region": args.region,
        "n_phenotypes_tested": len(df),
        "n_significant_fdr": n_sig,
        "fdr_threshold": args.fdr,
        "lambda_gc": lam,
        "median_pval": float(df[pcol].median()),
        "cis_qtl_path": str(args.cis_qtl),
    }]
    write_tsv(outdir / "meqtl_qc_summary.tsv", summary)

    # p-value histogram bins
    hist_counts, bin_edges = np.histogram(df[pcol].dropna(), bins=50, range=(0, 1))
    hist_rows = [
        {"bin_start": bin_edges[i], "bin_end": bin_edges[i + 1], "count": int(hist_counts[i])}
        for i in range(len(hist_counts))
    ]
    write_tsv(outdir / "pvalue_histogram.tsv", hist_rows)

    if "qval" in df.columns:
        lead = df.reset_index().rename(columns={"index": "phenotype_id"})
        if "phenotype_id" not in lead.columns and lead.columns[0] != "phenotype_id":
            lead = lead.rename(columns={lead.columns[0]: "phenotype_id"})
        cols = [c for c in [
            "phenotype_id", "num_var", "beta_shape1", "beta_shape2", "true_df",
            "pval_true_df", "variant_id", "start_distance", "end_distance", "tss_distance",
            "ma_samples", "ma_count", "af", "maf", "ref_factor", "pval_nominal", "slope", "slope_se", "pval_perm",
            "pval_beta", "qval", "pval_nominal_threshold",
        ] if c in lead.columns]
        lead[cols].to_csv(outdir / "lead_snp_per_cpg.tsv.gz", sep="\t", index=False)

    print(f"QC summary written to {outdir}; lambda_gc={lam:.3f}; n_sig={n_sig}")


if __name__ == "__main__":
    main()
