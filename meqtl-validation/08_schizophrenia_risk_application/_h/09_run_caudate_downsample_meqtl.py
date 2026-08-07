#!/usr/bin/env python3
"""Tier A: retest caudate SCZ risk–CpG pairs on downsampled donor sets.

Loads phenotypes/dosages once; for each replicate, subsets samples, residualizes
on M3a covariates, and applies BH-FDR over the same pair family as the full
caudate Analysis 3 table.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
P7 = PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m"
PHASE1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping/caudate/_m"
SEED = 20260801


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=str(P7 / "caudate_downsample"))
    p.add_argument("--pairs-tsv", default=str(P7 / "caudate/risk_variant_cpg_meqtl.tsv.gz"))
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--maf", type=float, default=0.05)
    p.add_argument("--max-reps", type=int, default=0, help="0 = all manifest replicates")
    return p.parse_args()


def load_covariates(path: Path, sample_order: list[str]) -> pd.DataFrame:
    cov = pd.read_csv(path, sep="\t", index_col=0)
    sample_set = set(sample_order)
    if set(cov.index.astype(str)).intersection(sample_set):
        cov.index = cov.index.astype(str)
    elif set(cov.columns.astype(str)).intersection(sample_set):
        cov = cov.T
        cov.index = cov.index.astype(str)
    else:
        raise SystemExit("No overlapping covariate/phenotype samples")
    return cov.loc[sample_order].apply(pd.to_numeric, errors="coerce")


def residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    n = X.shape[0]
    Xd = np.column_stack([np.ones(n), X])
    beta, _, _, _ = np.linalg.lstsq(Xd, y, rcond=None)
    return y - Xd @ beta


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    q_sorted = ranked * n / (np.arange(1, n + 1))
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)
    out = np.empty(n)
    out[order] = q_sorted
    q[ok] = out
    return q


def associate_pairs(
    Y: np.ndarray,
    G: np.ndarray,
    X: np.ndarray,
    pair_var_idx: np.ndarray,
    pair_pheno_idx: np.ndarray,
    maf_min: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized residualized association for pre-indexed pairs.

    Y: n×p phenotypes; G: n×v dosages; X: n×k covariates.
    Returns beta, se, pval, maf_ok mask for each pair.
    """
    n = Y.shape[0]
    Y_res = residualize(Y, X)
    G_res = residualize(G, X)
    # per-variant MAF on raw dosages (before residualize)
    maf = np.nanmean(G, axis=0) / 2.0
    maf = np.minimum(maf, 1.0 - maf)
    maf_ok_var = np.isfinite(maf) & (maf >= maf_min)

    # Center
    Yc = Y_res - Y_res.mean(axis=0, keepdims=True)
    Gc = G_res - G_res.mean(axis=0, keepdims=True)
    g = Gc[:, pair_var_idx]
    y = Yc[:, pair_pheno_idx]
    gss = np.sum(g * g, axis=0)
    yss = np.sum(y * y, axis=0)
    cov_gy = np.sum(g * y, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        slope = np.where(gss > 0, cov_gy / gss, np.nan)
        r_xy = np.where((gss > 0) & (yss > 0), cov_gy / np.sqrt(gss * yss), np.nan)
    df = n - 2 - X.shape[1]
    if df < 5:
        df = n - 2
    with np.errstate(invalid="ignore", divide="ignore"):
        tstat = r_xy * np.sqrt(df / np.clip(1 - r_xy ** 2, 1e-12, None))
    pval = 2 * stats.t.sf(np.abs(tstat), df)
    resid_ss = yss - (cov_gy ** 2) / np.clip(gss, 1e-12, None)
    se = np.sqrt(np.clip(resid_ss / (df * np.clip(gss, 1e-12, None)), 0, None))
    maf_ok = maf_ok_var[pair_var_idx] & np.isfinite(pval) & (gss > 1e-12)
    return slope, se, pval, maf_ok


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    rep_dir = outdir / "replicate_summaries"
    rep_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(outdir / "downsample_replicate_manifest.tsv", sep="\t")
    if args.max_reps > 0:
        manifest = manifest.head(args.max_reps).copy()

    pairs = pd.read_csv(args.pairs_tsv, sep="\t", compression="infer")
    for c in ["index_snp", "phenotype_id", "variant_id", "vmr_id"]:
        pairs[c] = pairs[c].astype(str)
    pairs["locus_id"] = pairs["locus_id"]
    n_pairs = len(pairs)
    print(f"Loaded {n_pairs} pairs from {args.pairs_tsv}")

    # Prioritized loci for tracking
    prior_path = P7 / "prioritized/prioritized_loci.tsv"
    prior_snps = set()
    prior_loci = set()
    if prior_path.exists():
        pr = pd.read_csv(prior_path, sep="\t")
        prior_snps = set(pr["index_snp"].astype(str))
        prior_loci = set(pr["locus_id"].dropna().astype(str))

    from tensorqtl import pgen
    import tensorqtl

    pheno_bed = PHASE1 / "prepared/cpg_phenotypes.all_autosomes.bed.gz"
    cov_path = PHASE1 / "prepared/covariates.txt"
    geno_prefix = str(PHASE1 / "genotypes/meqtl_AA")

    phenotype_df, _ = tensorqtl.read_phenotype_bed(str(pheno_bed))
    phenotype_df.columns = phenotype_df.columns.astype(str)
    need_pheno = sorted(set(pairs["phenotype_id"]).intersection(phenotype_df.index))
    phenotype_df = phenotype_df.loc[need_pheno]

    probe = pgen.PgenReader(geno_prefix)
    geno_ids = set(map(str, probe.sample_ids))
    samples_req = [s for s in phenotype_df.columns if s in geno_ids]
    cov_all = load_covariates(cov_path, samples_req)
    samples_req = [s for s in samples_req if s in cov_all.index]
    pgr = pgen.PgenReader(geno_prefix, select_samples=samples_req)
    samples_full = list(map(str, pgr.sample_ids))
    phenotype_df = phenotype_df[samples_full]
    cov_all = cov_all.loc[samples_full]
    if cov_all.isnull().any().any():
        cov_all = cov_all.dropna(axis=1, how="any")

    # Impute phenotype missingness once on full matrix
    Y_full = phenotype_df.to_numpy(dtype=float).T  # n × p
    col_mean = np.nanmean(Y_full, axis=0)
    inds = np.where(np.isnan(Y_full))
    if inds[0].size:
        Y_full[inds] = np.take(col_mean, inds[1])
    pheno_ids = list(phenotype_df.index)
    pheno_index = {p: i for i, p in enumerate(pheno_ids)}

    # Load dosages for unique variants
    uniq_vars = sorted(pairs["variant_id"].unique())
    panel = set(map(str, pgr.variant_ids))
    uniq_vars = [v for v in uniq_vars if v in panel]
    var_index = {v: i for i, v in enumerate(uniq_vars)}
    G_full = np.column_stack([
        np.asarray(pgr.read_dosages(v), dtype=float).ravel() for v in uniq_vars
    ])
    if np.isnan(G_full).any():
        # column-wise mean impute
        for j in range(G_full.shape[1]):
            m = np.nanmean(G_full[:, j])
            G_full[np.isnan(G_full[:, j]), j] = m

    # Filter pairs to those with loadable variant + phenotype
    pairs = pairs[
        pairs["variant_id"].isin(var_index) & pairs["phenotype_id"].isin(pheno_index)
    ].reset_index(drop=True)
    pair_var_idx = pairs["variant_id"].map(var_index).to_numpy(dtype=int)
    pair_pheno_idx = pairs["phenotype_id"].map(pheno_index).to_numpy(dtype=int)
    print(f"Retesting {len(pairs)} pairs; phenotypes={len(pheno_ids)}; variants={len(uniq_vars)}")

    sample_pos = {s: i for i, s in enumerate(samples_full)}
    X_full = cov_all.to_numpy(dtype=float)

    rep_rows = []
    prior_rows = []
    # Keep full-N reference counts (from original table before filtering)
    pairs_ref = pd.read_csv(args.pairs_tsv, sep="\t", compression="infer")
    if "significant_fdr" in pairs_ref.columns:
        sigmask = pairs_ref["significant_fdr"].astype(bool)
        full_sig = int(sigmask.sum())
        full_loci = int(pairs_ref.loc[sigmask, "locus_id"].nunique())
    else:
        full_sig = np.nan
        full_loci = np.nan

    for _, mrow in manifest.iterrows():
        rep = int(mrow["replicate"])
        samp = Path(mrow["path"]).read_text().strip().splitlines()
        idx = np.array([sample_pos[s] for s in samp if s in sample_pos], dtype=int)
        if len(idx) != int(mrow["n_samples"]):
            print(f"WARNING rep {rep}: expected {mrow['n_samples']} got {len(idx)}")
        Y = Y_full[idx]
        G = G_full[idx]
        X = X_full[idx]
        beta, se, pval, ok = associate_pairs(Y, G, X, pair_var_idx, pair_pheno_idx, args.maf)
        qval = np.full(len(pairs), np.nan)
        qval[ok] = bh_fdr(pval[ok])
        sig = ok & np.isfinite(qval) & (qval <= args.fdr)

        out = pairs[["index_snp", "phenotype_id", "variant_id", "locus_id", "vmr_id"]].copy()
        out["beta"] = beta
        out["se"] = se
        out["pval_nominal"] = pval
        out["qval"] = qval
        out["significant_fdr"] = sig
        out["maf_ok"] = ok
        out["n"] = len(idx)
        out["replicate"] = rep

        n_sig = int(sig.sum())
        n_loci = int(out.loc[sig, "locus_id"].nunique()) if sig.any() else 0
        n_vmr = int(out.loc[sig, "vmr_id"].nunique()) if sig.any() else 0
        n_snp = int(out.loc[sig, "index_snp"].nunique()) if sig.any() else 0
        rep_rows.append({
            "replicate": rep,
            "n_samples": len(idx),
            "n_tests": int(ok.sum()),
            "n_sig_fdr": n_sig,
            "n_sig_loci": n_loci,
            "n_sig_vmrs": n_vmr,
            "n_sig_index_snps": n_snp,
            "n_sig_full_caudate": full_sig,
            "n_sig_loci_full_caudate": full_loci,
        })

        # Prioritized locus tracking
        for snp in sorted(prior_snps):
            sub = out[out["index_snp"] == snp]
            if sub.empty:
                continue
            ssig = sub["significant_fdr"]
            i = sub["pval_nominal"].idxmin()
            prior_rows.append({
                "replicate": rep,
                "index_snp": snp,
                "locus_id": sub.at[i, "locus_id"],
                "n_tests": int(sub["maf_ok"].sum()),
                "n_sig_fdr": int(ssig.sum()),
                "min_pval": float(sub["pval_nominal"].min()),
                "min_qval": float(sub["qval"].min()) if sub["qval"].notna().any() else np.nan,
                "best_beta": float(sub.at[i, "beta"]),
                "best_se": float(sub.at[i, "se"]),
                "any_sig": bool(ssig.any()),
            })

        # Save compact sig hits only
        hits = out.loc[sig, [
            "replicate", "index_snp", "locus_id", "phenotype_id", "variant_id",
            "vmr_id", "beta", "se", "pval_nominal", "qval", "n",
        ]]
        hits.to_csv(rep_dir / f"sig_hits_rep{rep:03d}.tsv.gz", sep="\t", index=False, compression="gzip")
        print(f"rep {rep:03d}: n={len(idx)} tests={int(ok.sum())} sig={n_sig} loci={n_loci}")

    write_tsv(outdir / "downsample_replicate_results.tsv", rep_rows)
    if prior_rows:
        pd.DataFrame(prior_rows).to_csv(
            outdir / "downsample_prioritized_locus_results.tsv.gz",
            sep="\t", index=False, compression="gzip",
        )
    write_tsv(outdir / "downsample_run_summary.tsv", [{
        "n_reps": len(rep_rows),
        "n_pairs_input": n_pairs,
        "n_pairs_tested": len(pairs),
        "median_n_sig_fdr": float(np.median([r["n_sig_fdr"] for r in rep_rows])),
        "mean_n_sig_fdr": float(np.mean([r["n_sig_fdr"] for r in rep_rows])),
        "median_n_sig_loci": float(np.median([r["n_sig_loci"] for r in rep_rows])),
        "mean_n_sig_loci": float(np.mean([r["n_sig_loci"] for r in rep_rows])),
        "n_sig_full_caudate": full_sig,
        "fdr": args.fdr,
        "maf": args.maf,
    }])
    print(f"Done. Median sig pairs={np.median([r['n_sig_fdr'] for r in rep_rows]):.1f}")


if __name__ == "__main__":
    main()
