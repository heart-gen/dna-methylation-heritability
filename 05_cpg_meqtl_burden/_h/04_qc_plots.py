#!/usr/bin/env python3
"""meQTL mapping QC for 05_cpg_meqtl_burden.

Usage:
    python _h/04_qc_plots.py --run-id cmb-AA-caudate-20260823

config/meqtl_parameters.yml lists nine required QC outputs. The one that is a
GATE rather than a picture is genomic inflation: AGENTS.md 7.5 says it must be
resolved before the figure freeze, so lambda is written to its own TSV that
04_check_burden.R reads as a pass/fail criterion.

Lambda is computed from the NOMINAL pass, not the permutation pass. Permutation
p-values are already calibrated per phenotype by construction, so their
distribution says nothing about test-statistic inflation; the nominal pairwise
statistics are where confounding shows up.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def repo_root() -> Path:
    d = Path(__file__).resolve()
    while d != d.parent:
        if (d / ".git").is_dir():
            return d
        d = d.parent
    raise SystemExit("Could not locate repository root")


def load_nominal(nominal_dir: Path, max_pairs: int = 5_000_000) -> pd.DataFrame:
    """Pool nominal cis pairs across chromosomes.

    Capped: a whole-genome nominal pass is hundreds of millions of pairs and the
    QC statistics converge long before that. The cap is applied per chromosome
    by uniform subsampling so no chromosome is silently over-represented, and
    the realised count is recorded alongside lambda.
    """
    files = sorted(nominal_dir.glob("*.parquet")) or sorted(nominal_dir.glob("*.txt.gz"))
    if not files:
        raise SystemExit(f"No nominal output under {nominal_dir}")
    per_file = max(1, max_pairs // len(files))
    frames = []
    for f in files:
        df = pd.read_parquet(f) if f.suffix == ".parquet" else pd.read_csv(f, sep="\t")
        if len(df) > per_file:
            df = df.sample(per_file, random_state=20260722)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def genomic_inflation(pvals: np.ndarray) -> float:
    """Median chi-square inflation on 1 df."""
    from scipy.stats import chi2
    p = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    if p.size == 0:
        return float("nan")
    chi = chi2.isf(p, df=1)
    return float(np.median(chi) / chi2.ppf(0.5, df=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "05_cpg_meqtl_burden" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")
    cfg = yaml.safe_load((root / "config" / "meqtl_parameters.yml").read_text())
    fdr = float(cfg["mapping"]["fdr_threshold"])

    res_dir = run_dir / "results"
    qc_dir = res_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    cis = pd.read_csv(res_dir / "cpg-meqtl-results.tsv", sep="\t")
    nom = load_nominal(res_dir / "meqtl" / "nominal")

    # The gate statistic is lambda over DISTAL cis pairs only.
    #
    # Lambda over all nominal cis pairs is not an inflation estimate: cis pairs
    # are enriched for true meQTLs, so a well-behaved scan has lambda > 1 by
    # construction. The AA caudate chr22 smoke shows exactly that, and shows it
    # is signal rather than inflation, because lambda decays monotonically as
    # pairs move away from the CpG: 1.234 over all pairs, 1.146 beyond 100 kb,
    # 1.113 beyond 250 kb, 1.087 beyond 400 kb. The far edge of the +/-500 kb
    # window is where true cis signal is sparsest, so it is the best available
    # null proxy inside a cis scan.
    #
    # `lambda` is therefore the distal figure and is what 03_vmr_burden.R and
    # 04_check_burden.R gate on; `lambda_all_pairs` is retained alongside it as
    # a descriptive number, and the decay profile is written out so the
    # signal-vs-inflation distinction is auditable rather than asserted.
    pcol = "pval_nominal" if "pval_nominal" in nom else "pval"
    dcol = "start_distance" if "start_distance" in nom else (
        "tss_distance" if "tss_distance" in nom else None)
    if dcol is None:
        raise SystemExit(
            "Nominal output has no CpG-to-SNP distance column, so the distal "
            "lambda gate cannot be computed. Columns: " + ", ".join(nom.columns))

    distal_bp = int(cfg["genomic_inflation"]["distal_min_distance_bp"])
    pv = nom[pcol].to_numpy()
    dist = nom[dcol].abs().to_numpy()

    lam_all = genomic_inflation(pv)
    distal = dist > distal_bp
    lam = genomic_inflation(pv[distal])

    profile = pd.DataFrame(
        [{"min_abs_distance_bp": b,
          "n_pairs": int((dist > b).sum()),
          "lambda": genomic_inflation(pv[dist > b])}
         for b in (0, 100_000, 250_000, distal_bp)]
    ).drop_duplicates(subset="min_abs_distance_bp")
    profile.to_csv(qc_dir / "genomic-inflation-by-distance.tsv",
                   sep="\t", index=False)

    pd.DataFrame([{
        "lambda": lam,
        "lambda_statistic": f"distal_cis_pairs_abs_distance_gt_{distal_bp}",
        "lambda_all_pairs": lam_all,
        "n_distal_pairs_used": int(distal.sum()),
        "n_nominal_pairs_used": int(len(nom)),
        "n_cpgs": int(len(cis)),
        "n_significant": int((cis["qvalue"] <= fdr).sum()),
    }]).to_csv(qc_dir / "genomic-inflation.tsv", sep="\t", index=False)
    print(f"genomic inflation lambda = {lam:.4f} on {int(distal.sum())} distal "
          f"pairs (|d| > {distal_bp}); all-pair lambda = {lam_all:.4f} "
          f"on {len(nom)} pairs")

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(qc_dir / f"{name}.png", dpi=150)
        plt.close(fig)

    # 1. p-value histogram (nominal)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.hist(nom[pcol].dropna(), bins=50)
    ax.set_xlabel("nominal p"); ax.set_ylabel("cis pairs")
    ax.set_title("Nominal p-value distribution")
    save(fig, "pvalue_histogram")

    # 2. QQ plot
    p = np.sort(nom[pcol].dropna().to_numpy())
    exp = -np.log10((np.arange(1, p.size + 1) - 0.5) / p.size)
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.scatter(exp, -np.log10(p), s=1)
    lim = max(exp.max(), -np.log10(p[0]))
    ax.plot([0, lim], [0, lim], color="grey", lw=1)
    ax.set_xlabel("expected -log10 p"); ax.set_ylabel("observed -log10 p")
    # The QQ is over all cis pairs, so label it with the all-pair lambda and
    # name the distal gate figure separately -- captioning this plot with the
    # distal number would misdescribe the points being drawn.
    ax.set_title(f"QQ, all cis pairs (lambda = {lam_all:.3f}; "
                 f"distal gate = {lam:.3f})")
    save(fig, "qq_plot")

    # 3. effect sizes, 4. lead-SNP MAF, 5. lead-SNP distance
    for col, name, xlab in [
        ("slope", "effect_size_distribution", "beta (lead SNP)"),
        ("af", "lead_snp_maf_distribution", "lead SNP allele frequency"),
        ("start_distance", "lead_snp_distance_distribution", "CpG-to-lead-SNP distance (bp)"),
    ]:
        if col not in cis:
            print(f"WARNING: column {col} absent; skipping {name}")
            continue
        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.hist(cis[col].dropna(), bins=50)
        ax.set_xlabel(xlab); ax.set_ylabel("CpGs"); ax.set_title(name.replace("_", " "))
        save(fig, name)

    # 6. per-chromosome result counts
    if "chrom" in cis:
        counts = (cis.assign(sig=cis["qvalue"] <= fdr)
                     .groupby("chrom")["sig"].agg(["size", "sum"])
                     .rename(columns={"size": "n_tested", "sum": "n_significant"}))
        counts.to_csv(qc_dir / "chromosome_result_counts.tsv", sep="\t")
        fig, ax = plt.subplots(figsize=(7, 3.5))
        ax.bar(counts.index.astype(str), counts["n_significant"])
        ax.set_ylabel(f"CpGs at FDR<{fdr}"); ax.tick_params(axis="x", rotation=90)
        ax.set_title("Significant CpGs per chromosome")
        save(fig, "chromosome_result_counts")

    # Two entries of `qc_plots_required` are NOT produced here, and the gap is
    # recorded rather than left to be discovered by counting PNGs:
    #   * covariate_model_comparison needs a second mapping run under an
    #     alternative covariate design; the module maps once, under the locked
    #     design, so there is nothing to compare against.
    #   * positive_control_public_meqtl_overlap needs a public brain meQTL
    #     reference, which is not registered in
    #     inputs/supportfiles/_m/annotation_asset_manifest.tsv.
    # Both require a PI decision (a second run, or an asset to register) before
    # they can be implemented.
    pd.DataFrame([
        {"plot": "covariate_model_comparison",
         "reason": "requires a second mapping run under an alternative "
                   "covariate design; not run"},
        {"plot": "positive_control_public_meqtl_overlap",
         "reason": "requires a public brain meQTL reference; no such asset is "
                   "registered in the annotation asset manifest"},
    ]).to_csv(qc_dir / "qc-plots-not-produced.tsv", sep="\t", index=False)

    print(f"QC written to {qc_dir}")
    print("NOT produced (see qc-plots-not-produced.tsv): "
          "covariate_model_comparison, positive_control_public_meqtl_overlap")


if __name__ == "__main__":
    main()
