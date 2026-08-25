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


# --------------------------------------------------------- positive control
# `positive_control_public_meqtl_overlap` asks a QC question, not a validation
# question: does this scan recover cis-meQTLs that are already known? A
# positive control does NOT require an independent cohort -- overlap with the
# discovery donors would if anything make it a better control, because a low
# recovery rate then indicts the scan rather than the reference. That is why
# the Phase 3 independence exclusions do not carry over here.
#
# Two rules are inherited from the Phase 3 repair (PHASE3_DIAGNOSIS.md) and
# must not be relaxed:
#   1. The denominator is the resource's ASSAYED universe, never its results
#      table. A WGBS CpG absent from a 450K catalog is unobserved, not a
#      negative. The harmonized tables carry the full 485,441-probe universe
#      with an `external_assayed` flag, so the inner join below IS the
#      restriction to assayed CpGs.
#   2. Matching is an exact (chrom, pos_1based) join in hg38, the same
#      convention 03_harmonize_external_meqtls.py used, including its
#      "keep the supported / best-FDR row per position" dedup.
#
# The reported statistic is the recovery rate among externally supported CpGs,
# contrasted against the rate among assayed-but-unsupported ones. The contrast
# is what makes it interpretable: a scan that calls everything significant
# would score a perfect recovery rate and is caught by the second column.


def _resources_from_manifest(root: Path) -> list[dict]:
    """Public meQTL references registered for this QC plot.

    Selection is by the manifest's notes field rather than a hard-coded list,
    so registering a further resource does not require editing this script.
    """
    man_path = root / "inputs" / "supportfiles" / "_m" / "annotation_asset_manifest.tsv"
    if not man_path.is_file():
        return []
    man = pd.read_csv(man_path, sep="\t")
    hits = man[man["notes"].fillna("").str.contains(
        "positive_control_public_meqtl_overlap", regex=False)]
    out = []
    for _, r in hits.iterrows():
        p = Path(str(r["path"]))
        if not p.is_file():
            print(f"WARNING: registered asset missing on disk, skipping: {p}")
            continue
        out.append({"filename": str(r["filename"]), "path": p})
    return out


def positive_control(cis: pd.DataFrame, root: Path, qc_dir: Path,
                     fdr: float, region: str, save) -> bool:
    """Overlap the scan against registered public brain meQTL catalogs.

    Returns True if at least one resource yielded a usable comparison.
    """
    resources = _resources_from_manifest(root)
    if not resources:
        print("positive control: no public meQTL resource registered in the "
              "annotation asset manifest; not produced")
        return False

    # cpg_id is 'chrN:pos' in hg38, 1-based -- verified against the Jaffe
    # universe on the chr22 smoke, where offset 0 matches 172 CpGs and offsets
    # -1/+1 match none.
    parts = cis["cpg_id"].astype(str).str.rsplit(":", n=1, expand=True)
    scan = cis.assign(chrom=parts[0], pos_1based=pd.to_numeric(parts[1], errors="coerce"))
    scan = scan.dropna(subset=["pos_1based"])
    scan["pos_1based"] = scan["pos_1based"].astype(int)
    scan["scan_significant"] = (scan["qvalue"] <= fdr).astype(int)

    rows = []
    for res in resources:
        ref = pd.read_csv(res["path"], sep="\t")
        needed = {"chrom", "pos_1based", "external_assayed", "external_meqtl_support"}
        if not needed.issubset(ref.columns):
            print(f"WARNING: {res['filename']} lacks {sorted(needed - set(ref.columns))}; "
                  "skipping")
            continue
        resource_id = (str(ref["resource_id"].iloc[0])
                       if "resource_id" in ref.columns and len(ref) else res["filename"])
        ref_tissue = (str(ref["tissue_region"].iloc[0])
                      if "tissue_region" in ref.columns and len(ref) else "unknown")

        assay = ref.loc[ref["external_assayed"] == 1].copy()
        assay["chrom"] = assay["chrom"].astype(str)
        assay["pos_1based"] = assay["pos_1based"].astype(int)
        if "external_fdr" in assay.columns:
            assay = assay.sort_values(["external_meqtl_support", "external_fdr"],
                                      ascending=[False, True])
        assay = assay.drop_duplicates(subset=["chrom", "pos_1based"], keep="first")

        m = scan.merge(assay[["chrom", "pos_1based", "external_meqtl_support"]],
                       on=["chrom", "pos_1based"], how="inner")
        if m.empty:
            print(f"WARNING: {resource_id} shares no assayed CpG with the scan; skipping")
            continue
        m["external_meqtl_support"] = m["external_meqtl_support"].astype(int)

        sup = m["external_meqtl_support"] == 1
        n_sup, n_unsup = int(sup.sum()), int((~sup).sum())
        rec_sup = float(m.loc[sup, "scan_significant"].mean()) if n_sup else float("nan")
        rec_unsup = float(m.loc[~sup, "scan_significant"].mean()) if n_unsup else float("nan")

        table = [[int(m.loc[sup, "scan_significant"].sum()), n_sup - int(m.loc[sup, "scan_significant"].sum())],
                 [int(m.loc[~sup, "scan_significant"].sum()), n_unsup - int(m.loc[~sup, "scan_significant"].sum())]]
        try:
            from scipy.stats import fisher_exact
            odds, pval = fisher_exact(table)
        except Exception:
            odds, pval = float("nan"), float("nan")

        # Tissue matching is descriptive, not a filter: a cross-tissue control
        # is weaker but still informative, and caudate has no matched resource
        # at all, so dropping unmatched pairs would leave that cell with none.
        matched = ref_tissue.lower() == str(region).lower()
        rows.append({
            "resource_id": resource_id,
            "resource_tissue": ref_tissue,
            "run_region": region,
            "tissue_match": "matched" if matched else "cross_tissue",
            "n_assayed_shared": len(m),
            "n_scan_cpgs": len(scan),
            "frac_scan_cpgs_assayed": len(m) / len(scan) if len(scan) else float("nan"),
            "n_external_supported": n_sup,
            "n_external_unsupported": n_unsup,
            "recovery_rate_supported": rec_sup,
            "recovery_rate_unsupported": rec_unsup,
            "odds_ratio": float(odds),
            "fisher_p": float(pval),
            "fdr_threshold": fdr,
        })
        print(f"positive control {resource_id} ({ref_tissue} vs {region}, "
              f"{'matched' if matched else 'cross-tissue'}): "
              f"{len(m)} shared assayed CpGs, recovery {rec_sup:.3f} supported "
              f"vs {rec_unsup:.3f} unsupported, OR={odds:.2f}, p={pval:.3g}")

    if not rows:
        return False
    out = pd.DataFrame(rows)
    out.to_csv(qc_dir / "positive_control_public_meqtl_overlap.tsv",
               sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(1.9 * len(out) + 2.6, 3.8))
    x = np.arange(len(out))
    w = 0.38
    ax.bar(x - w / 2, out["recovery_rate_supported"], w, label="external meQTL support")
    ax.bar(x + w / 2, out["recovery_rate_unsupported"], w, label="assayed, no support")
    for i, r in out.reset_index(drop=True).iterrows():
        ax.text(i - w / 2, r["recovery_rate_supported"], f"n={r['n_external_supported']}",
                ha="center", va="bottom", fontsize=7)
        ax.text(i + w / 2, r["recovery_rate_unsupported"], f"n={r['n_external_unsupported']}",
                ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.resource_id}\n({r.resource_tissue}, {r.tissue_match})"
                        for r in out.itertuples()], fontsize=7)
    ax.set_ylabel(f"fraction significant in this scan (FDR<{fdr})")
    ax.set_ylim(0, 1)
    ax.set_title(f"Positive control: recovery of known brain meQTLs ({region})")
    ax.legend(fontsize=7)
    save(fig, "positive_control_public_meqtl_overlap")
    return True


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

    # 7. positive control against registered public brain meQTL catalogs.
    region = str(cis["region"].iloc[0]) if "region" in cis and len(cis) else "unknown"
    produced_positive_control = positive_control(cis, root, qc_dir, fdr, region, save)

    # Whatever `qc_plots_required` still cannot be produced is recorded, rather
    # than left to be discovered by counting PNGs.
    #
    # covariate_model_comparison needs a second mapping run under an
    # alternative covariate design; the module maps once, under the locked
    # design, so there is nothing to compare against, and which alternative
    # design to use is a PI decision.
    not_produced = [
        {"plot": "covariate_model_comparison",
         "reason": "requires a second mapping run under an alternative "
                   "covariate design; not run"},
    ]
    if not produced_positive_control:
        not_produced.append(
            {"plot": "positive_control_public_meqtl_overlap",
             "reason": "no usable public brain meQTL reference: none registered "
                       "in the annotation asset manifest, or none shares an "
                       "assayed CpG with this scan"})
    pd.DataFrame(not_produced).to_csv(
        qc_dir / "qc-plots-not-produced.tsv", sep="\t", index=False)

    print(f"QC written to {qc_dir}")
    print("NOT produced (see qc-plots-not-produced.tsv): "
          + ", ".join(r["plot"] for r in not_produced))


if __name__ == "__main__":
    main()
