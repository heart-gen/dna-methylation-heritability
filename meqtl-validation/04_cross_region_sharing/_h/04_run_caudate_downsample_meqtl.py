#!/usr/bin/env python3
"""Phase 4 Experiment 3: retest caudate cis lead SNP–CpG pairs on N-matched donors.

Loads phenotypes/dosages once; for each downsample replicate, residualizes on
locked M3a covariates and applies BH-FDR over the same CpG family as the full
caudate Phase 1 lead table.

This is lead-SNP retention under donor downsampling (tractable architecture
sensitivity), not full cis remapping of each replicate. Documented as such in
outputs and the readiness memo.
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
PHASE1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping/caudate/_m"
OUTDIR = PROJECT / "meqtl-validation/04_cross_region_sharing/_m/caudate_downsample"
FDR = 0.05
MAF_MIN = 0.05


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=str(OUTDIR))
    p.add_argument(
        "--leads-tsv",
        default=str(PHASE1 / "tensorqtl/qc/lead_snp_per_cpg.tsv.gz"),
    )
    p.add_argument("--fdr", type=float, default=FDR)
    p.add_argument("--maf", type=float, default=MAF_MIN)
    p.add_argument("--max-reps", type=int, default=0, help="0 = all manifest replicates")
    p.add_argument(
        "--sig-only",
        action="store_true",
        help="Retest only full-N FDR-significant CpGs (faster; FDR among that subset)",
    )
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
    n = Y.shape[0]
    Y_res = residualize(Y, X)
    G_res = residualize(G, X)
    maf = np.nanmean(G, axis=0) / 2.0
    maf = np.minimum(maf, 1.0 - maf)
    maf_ok_var = np.isfinite(maf) & (maf >= maf_min)

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


def load_cpg_vmr_map() -> pd.DataFrame:
    maps = sorted((PHASE1 / "prepared").glob("cpg_vmr_map.chr*.tsv"))
    if not maps:
        return pd.DataFrame(columns=["phenotype_id", "vmr_id"])
    frames = [pd.read_csv(p, sep="\t") for p in maps]
    m = pd.concat(frames, ignore_index=True)
    m["phenotype_id"] = m["phenotype_id"].astype(str)
    m["vmr_id"] = m["vmr_id"].astype(str)
    return m[["phenotype_id", "vmr_id"]].drop_duplicates("phenotype_id")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    rep_dir = outdir / "replicate_summaries"
    rep_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(outdir / "downsample_replicate_manifest.tsv", sep="\t")
    if args.max_reps > 0:
        manifest = manifest.head(args.max_reps).copy()

    leads = pd.read_csv(args.leads_tsv, sep="\t", compression="infer")
    leads["phenotype_id"] = leads["phenotype_id"].astype(str)
    leads["variant_id"] = leads["variant_id"].astype(str)
    leads["qval"] = pd.to_numeric(leads["qval"], errors="coerce")
    leads["sig_full"] = leads["qval"].le(args.fdr).fillna(False)
    n_sig_full = int(leads["sig_full"].sum())
    n_tested_full = len(leads)
    if args.sig_only:
        leads = leads.loc[leads["sig_full"]].reset_index(drop=True)
    print(
        f"Loaded {len(leads)} lead pairs "
        f"(full-N sig={n_sig_full}/{n_tested_full}; sig_only={args.sig_only})"
    )

    cpg_vmr = load_cpg_vmr_map()
    leads = leads.merge(cpg_vmr, on="phenotype_id", how="left")

    from tensorqtl import pgen
    import tensorqtl

    pheno_bed = PHASE1 / "prepared/cpg_phenotypes.all_autosomes.bed.gz"
    cov_path = PHASE1 / "prepared/covariates.txt"
    geno_prefix = str(PHASE1 / "genotypes/meqtl_AA")

    phenotype_df, _ = tensorqtl.read_phenotype_bed(str(pheno_bed))
    phenotype_df.columns = phenotype_df.columns.astype(str)
    need_pheno = sorted(set(leads["phenotype_id"]).intersection(phenotype_df.index.astype(str)))
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

    Y_full = phenotype_df.to_numpy(dtype=float).T
    col_mean = np.nanmean(Y_full, axis=0)
    inds = np.where(np.isnan(Y_full))
    if inds[0].size:
        Y_full[inds] = np.take(col_mean, inds[1])
    pheno_ids = list(phenotype_df.index.astype(str))
    pheno_index = {p: i for i, p in enumerate(pheno_ids)}

    uniq_vars = sorted(leads["variant_id"].unique())
    panel = set(map(str, pgr.variant_ids))
    uniq_vars = [v for v in uniq_vars if v in panel]
    var_index = {v: i for i, v in enumerate(uniq_vars)}
    print(f"Loading dosages for {len(uniq_vars)} unique lead variants...")
    G_full = np.column_stack([
        np.asarray(pgr.read_dosages(v), dtype=float).ravel() for v in uniq_vars
    ])
    if np.isnan(G_full).any():
        for j in range(G_full.shape[1]):
            m = np.nanmean(G_full[:, j])
            G_full[np.isnan(G_full[:, j]), j] = m

    leads = leads[
        leads["variant_id"].isin(var_index) & leads["phenotype_id"].isin(pheno_index)
    ].reset_index(drop=True)
    pair_var_idx = leads["variant_id"].map(var_index).to_numpy(dtype=int)
    pair_pheno_idx = leads["phenotype_id"].map(pheno_index).to_numpy(dtype=int)
    print(
        f"Retesting {len(leads)} pairs; phenotypes={len(pheno_ids)}; "
        f"variants={len(uniq_vars)}; donors_full={len(samples_full)}"
    )

    sample_pos = {s: i for i, s in enumerate(samples_full)}
    X_full = cov_all.to_numpy(dtype=float)

    rep_rows = []
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
        qval = np.full(len(leads), np.nan)
        qval[ok] = bh_fdr(pval[ok])
        sig = ok & np.isfinite(qval) & (qval <= args.fdr)

        n_sig = int(sig.sum())
        n_vmr = int(leads.loc[sig, "vmr_id"].nunique()) if sig.any() and "vmr_id" in leads.columns else 0
        n_retain_full_sig = int(
            (sig & leads["sig_full"].to_numpy()).sum()
        ) if "sig_full" in leads.columns else np.nan
        frac_retain = (
            float(n_retain_full_sig / max(int(leads["sig_full"].sum()), 1))
            if "sig_full" in leads.columns else np.nan
        )
        rep_rows.append({
            "replicate": rep,
            "n_samples": len(idx),
            "n_tests": int(ok.sum()),
            "n_sig_fdr": n_sig,
            "n_sig_vmrs": n_vmr,
            "n_sig_full_caudate": n_sig_full,
            "n_full_sig_retained": n_retain_full_sig,
            "frac_full_sig_retained": frac_retain,
            "frac_cpgs_sig": float(n_sig / max(int(ok.sum()), 1)),
        })

        # Compact: save only counts per replicate (full hit tables are huge)
        print(
            f"rep {rep:03d}: n={len(idx)} tests={int(ok.sum())} "
            f"sig={n_sig} retain_full_sig={n_retain_full_sig}"
        )

    write_tsv(outdir / "downsample_replicate_results.tsv", rep_rows)
    write_tsv(outdir / "downsample_run_summary.tsv", [{
        "n_reps": len(rep_rows),
        "n_pairs_input": n_tested_full if not args.sig_only else n_sig_full,
        "n_pairs_tested": len(leads),
        "sig_only": bool(args.sig_only),
        "method": "lead_snp_retention_M3a_residualized",
        "median_n_sig_fdr": float(np.median([r["n_sig_fdr"] for r in rep_rows])),
        "mean_n_sig_fdr": float(np.mean([r["n_sig_fdr"] for r in rep_rows])),
        "median_n_sig_vmrs": float(np.median([r["n_sig_vmrs"] for r in rep_rows])),
        "median_frac_full_sig_retained": float(
            np.median([r["frac_full_sig_retained"] for r in rep_rows])
        ),
        "n_sig_full_caudate": n_sig_full,
        "fdr": args.fdr,
        "maf": args.maf,
    }])
    print(f"Done. Median sig CpGs={np.median([r['n_sig_fdr'] for r in rep_rows]):.1f}")


if __name__ == "__main__":
    main()
