#!/usr/bin/env python3
"""
Calibration diagnostics for cis-meQTL permutation results.

Significance is defined by Storey q-value (qval), matching the primary FDR family.
This script does not change significance calls; it diagnoses genomic inflation.

Outputs under {cis_qtl_dir}/qc/calibration/:
  - lambda_by_qval.tsv
  - qq_plot_{pval_col}.png / .pdf
  - calibration_summary.tsv
"""

from __future__ import annotations

import argparse, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402


# Prespecified q-value strata; significance uses qval <= fdr_threshold
DEFAULT_STRATA = [
    ("all_tested", None, None),
    ("nonsignificant_qval_gt_fdr", "gt_fdr", None),
    ("significant_qval_le_fdr", "le_fdr", None),
    ("near_null_qval_gt_0.5", 0.5, None),  # lower bound exclusive via gt
    ("weak_0.05_lt_qval_le_0.5", 0.05, 0.5),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True)
    p.add_argument("--cis-qtl", required=True, help="*.cis_qtl.txt.gz from TensorQTL")
    p.add_argument("--fdr", type=float, default=0.05, help="q-value significance threshold")
    p.add_argument(
        "--pval-cols",
        nargs="+",
        default=["pval_beta", "pval_perm", "pval_nominal"],
        help="P-value columns to calibrate (skip if absent)",
    )
    p.add_argument("--outdir", default="")
    return p.parse_args()


def genomic_inflation(pvals: np.ndarray) -> float:
    from scipy.stats import chi2

    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if p.size == 0:
        return float("nan")
    chisq = chi2.ppf(1.0 - p, 1)
    return float(np.median(chisq) / chi2.ppf(0.5, 1))


def stratum_mask(qval: pd.Series, name: str, fdr: float, lo: float | None, hi: float | None) -> pd.Series:
    q = pd.to_numeric(qval, errors="coerce")
    if name == "all_tested":
        return q.notna()
    if name == "nonsignificant_qval_gt_fdr":
        return q > fdr
    if name == "significant_qval_le_fdr":
        return q <= fdr
    if name == "near_null_qval_gt_0.5":
        return q > 0.5
    if name == "weak_0.05_lt_qval_le_0.5":
        return (q > 0.05) & (q <= 0.5)
    # generic lo/hi if extended later
    mask = q.notna()
    if lo is not None:
        mask &= q > lo
    if hi is not None:
        mask &= q <= hi
    return mask


def qq_arrays(pvals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    p = np.sort(p)
    n = p.size
    if n == 0:
        return np.array([]), np.array([])
    exp = np.arange(1, n + 1) / (n + 1)
    return -np.log10(exp), -np.log10(p)


def plot_qq(
    df: pd.DataFrame,
    pcol: str,
    qcol: str,
    fdr: float,
    out_prefix: Path,
    region: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    panels = [
        ("all tested", df[pcol].to_numpy()),
        (f"qval > {fdr:g} (NS)", df.loc[df[qcol] > fdr, pcol].to_numpy()),
        (f"qval ≤ {fdr:g} (sig)", df.loc[df[qcol] <= fdr, pcol].to_numpy()),
    ]
    for ax, (title, pvals) in zip(axes, panels):
        x, y = qq_arrays(pvals)
        if x.size == 0:
            ax.set_title(title)
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            continue
        lim = max(float(np.nanmax(x)), float(np.nanmax(y)), 1.0)
        # thin for large n
        if x.size > 20000:
            idx = np.linspace(0, x.size - 1, 20000).astype(int)
            ax.scatter(x[idx], y[idx], s=4, c="#2c7fb8", alpha=0.35, linewidths=0)
        else:
            ax.scatter(x, y, s=6, c="#2c7fb8", alpha=0.4, linewidths=0)
        ax.plot([0, lim], [0, lim], color="#333333", lw=1.0, ls="--")
        lam = genomic_inflation(pvals)
        ax.set_xlim(0, lim * 1.02)
        ax.set_ylim(0, lim * 1.02)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"Expected $-\log_{10}(p)$")
        ax.set_ylabel(r"Observed $-\log_{10}(p)$")
        ax.set_title(f"{title}\nλ_GC={lam:.2f}; n={np.isfinite(pvals).sum():,}")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(f"{region} cis-meQTL QQ ({pcol}); significance by qval≤{fdr:g}", fontsize=11)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_prefix}.png", dpi=200)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    cis_path = Path(args.cis_qtl)
    outdir = Path(args.outdir) if args.outdir else cis_path.parent / "qc" / "calibration"
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(cis_path, sep="\t")
    # TensorQTL may use phenotype_id as index column without name
    if "phenotype_id" not in df.columns:
        first = df.columns[0]
        if first.startswith("chr") or first == "Unnamed: 0" or "phenotype" in first.lower():
            df = df.rename(columns={first: "phenotype_id"})
        else:
            df = df.reset_index().rename(columns={"index": "phenotype_id"})

    if "qval" not in df.columns:
        raise SystemExit("cis_qtl table lacks qval; cannot stratify by significance")

    df["qval"] = pd.to_numeric(df["qval"], errors="coerce")
    pcols = [c for c in args.pval_cols if c in df.columns]
    if not pcols:
        raise SystemExit(f"None of requested p-value columns found: {args.pval_cols}")

    lambda_rows: list[dict] = []
    for pcol in pcols:
        df[pcol] = pd.to_numeric(df[pcol], errors="coerce")
        for name, lo, hi in DEFAULT_STRATA:
            mask = stratum_mask(df["qval"], name, args.fdr, lo, hi)
            p = df.loc[mask, pcol].to_numpy(dtype=float)
            lambda_rows.append({
                "region": args.region,
                "pval_column": pcol,
                "stratum": name,
                "fdr_threshold": args.fdr,
                "significance_rule": f"qval <= {args.fdr}",
                "n_cpgs": int(mask.sum()),
                "lambda_gc": genomic_inflation(p),
                "median_pval": float(np.nanmedian(p)) if np.isfinite(p).any() else float("nan"),
                "mean_qval": float(df.loc[mask, "qval"].mean()) if mask.any() else float("nan"),
            })
        plot_qq(
            df, pcol, "qval", args.fdr,
            outdir / f"qq_plot_{pcol}",
            args.region,
        )

    write_tsv(outdir / "lambda_by_qval.tsv", lambda_rows)

    # Compact primary summary emphasizing qval-defined significance
    primary = [r for r in lambda_rows if r["pval_column"] == pcols[0]]
    summary = [{
        "region": args.region,
        "significance_rule": f"qval <= {args.fdr}",
        "n_tested": int(len(df)),
        "n_significant": int((df["qval"] <= args.fdr).sum()),
        "n_nonsignificant": int((df["qval"] > args.fdr).sum()),
        "primary_pval_column": pcols[0],
        "lambda_gc_all": next(r["lambda_gc"] for r in primary if r["stratum"] == "all_tested"),
        "lambda_gc_nonsignificant": next(
            r["lambda_gc"] for r in primary if r["stratum"] == "nonsignificant_qval_gt_fdr"
        ),
        "lambda_gc_significant": next(
            r["lambda_gc"] for r in primary if r["stratum"] == "significant_qval_le_fdr"
        ),
        "lambda_gc_near_null_qval_gt_0.5": next(
            r["lambda_gc"] for r in primary if r["stratum"] == "near_null_qval_gt_0.5"
        ),
        "outdir": str(outdir),
    }]
    write_tsv(outdir / "calibration_summary.tsv", summary)
    print(
        f"{args.region}: calibration written to {outdir} "
        f"(λ_all={summary[0]['lambda_gc_all']:.3f}, "
        f"λ_NS={summary[0]['lambda_gc_nonsignificant']:.3f}, "
        f"λ_nearNull={summary[0]['lambda_gc_near_null_qval_gt_0.5']:.3f})"
    )


if __name__ == "__main__":
    main()
