#!/usr/bin/env python3
"""
Estimate region-specific methylation latent factors (PCA on M0 residuals).

PEER is not installed in the genomics env; this implements the planned PEER-style
latent set as PCA on the CpG phenotype matrix after regressing baseline
covariates (age, sex, diagnosis, snpPC1–5). Factors are written for k=5/10/15.

Outputs under:
  {region}/_m/covariate_sensitivity/latent/
    latent_factors_k{K}.tsv   (samples × factors)
    latent_estimation_summary.tsv
"""

from __future__ import annotations

import argparse, sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import load_paths, write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
MAX_K = 15
DEFAULT_SUBSAMPLE = 50000
SEED = 20260730


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True, choices=["caudate", "dlpfc", "hippocampus"])
    p.add_argument("--subsample-cpgs", type=int, default=DEFAULT_SUBSAMPLE)
    p.add_argument("--seed", type=int, default=SEED)
    return p.parse_args()


def load_phenotype_matrix(bed_path: Path, sample_ids: list[str]) -> pd.DataFrame:
    """Load phenotype BED/BED.gz as CpGs × samples (float)."""
    # tensorqtl BED: #chr start end phenotype_id sample...
    header = pd.read_csv(bed_path, sep="\t", nrows=0, compression="infer")
    cols = list(header.columns)
    id_col = "phenotype_id" if "phenotype_id" in cols else cols[3]
    sample_cols = [c for c in sample_ids if c in cols]
    missing = [s for s in sample_ids if s not in cols]
    if missing:
        raise SystemExit(f"Phenotype BED missing samples: {missing[:5]}")
    df = pd.read_csv(
        bed_path,
        sep="\t",
        compression="infer",
        usecols=[id_col] + sample_cols,
    )
    df = df.set_index(id_col)
    return df[sample_ids].astype(float)


def residualize(Y: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Residualize each phenotype (column of Y.T style: samples × features) on cov."""
    # cov: n_samples × n_cov; Y: n_cpgs × n_samples → work samples × cpgs
    X = np.column_stack([np.ones(cov.shape[0]), cov])
    beta, _, _, _ = np.linalg.lstsq(X, Y.T, rcond=None)
    fitted = X @ beta
    return (Y.T - fitted).T


def main() -> None:
    args = parse_args()
    paths = load_paths()
    region = args.region
    base = PROJECT / "meqtl-validation" / "01_cpg_meqtl_mapping" / region / "_m"
    prep = base / "prepared"
    outdir = base / "covariate_sensitivity" / "latent"
    outdir.mkdir(parents=True, exist_ok=True)

    cov0 = pd.read_csv(prep / "covariates.txt", sep="\t", index_col=0)
    cov0.index = cov0.index.astype(str)
    sample_ids = cov0.index.tolist()

    bed = prep / "cpg_phenotypes.all_autosomes.bed.gz"
    if not bed.exists():
        bed = prep / "cpg_phenotypes.all_autosomes.bed"
    if not bed.exists():
        raise SystemExit(f"Missing phenotype BED under {prep}")

    print(f"Loading phenotypes from {bed}")
    phen = load_phenotype_matrix(bed, sample_ids)
    n_cpg = phen.shape[0]
    rng = np.random.default_rng(args.seed)
    if args.subsample_cpgs and n_cpg > args.subsample_cpgs:
        idx = rng.choice(n_cpg, size=args.subsample_cpgs, replace=False)
        phen = phen.iloc[idx]
        print(f"Subsampled {args.subsample_cpgs} / {n_cpg} CpGs for latent estimation")

    # Drop zero-variance CpGs
    sd = phen.std(axis=1, ddof=0)
    phen = phen.loc[sd > 0]
    Y = phen.to_numpy(dtype=float)
    # Impute any remaining NaNs with CpG mean
    if np.isnan(Y).any():
        col_means = np.nanmean(Y, axis=1)
        inds = np.where(np.isnan(Y))
        Y[inds] = np.take(col_means, inds[0])

    cov = cov0.to_numpy(dtype=float)
    Y_res = residualize(Y, cov)
    # Standardize features (CpGs) across samples before PCA on samples×CpGs
    Z = Y_res.T  # samples × cpgs
    Z = (Z - Z.mean(axis=0)) / np.clip(Z.std(axis=0, ddof=0), 1e-8, None)

    pca = PCA(n_components=MAX_K, random_state=args.seed)
    scores = pca.fit_transform(Z)
    factor_names = [f"methPC{i}" for i in range(1, MAX_K + 1)]
    factors = pd.DataFrame(scores, index=sample_ids, columns=factor_names)
    factors.index.name = "id"

    for k in (5, 10, 15):
        out = outdir / f"latent_factors_k{k}.tsv"
        factors.iloc[:, :k].to_csv(out, sep="\t")
        print(f"Wrote {out}")

    write_tsv(
        outdir / "latent_estimation_summary.tsv",
        [{
            "region": region,
            "method": "PCA_on_M0_residuals",
            "n_samples": len(sample_ids),
            "n_cpgs_used": int(phen.shape[0]),
            "n_cpgs_available": n_cpg,
            "max_k": MAX_K,
            "seed": args.seed,
            "variance_explained_k5": float(pca.explained_variance_ratio_[:5].sum()),
            "variance_explained_k10": float(pca.explained_variance_ratio_[:10].sum()),
            "variance_explained_k15": float(pca.explained_variance_ratio_[:15].sum()),
            "note": "PEER package unavailable; PCA residual factors used as planned PEER-style latents",
        }],
    )
    print(f"Latent estimation complete for {region}")


if __name__ == "__main__":
    main()
