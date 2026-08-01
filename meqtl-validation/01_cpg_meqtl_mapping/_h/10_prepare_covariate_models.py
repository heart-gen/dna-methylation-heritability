#!/usr/bin/env python3
"""
Build covariate matrices for Phase 1 sensitivity models M0–M5.

Writes under {region}/_m/covariate_sensitivity/:
  covariates_M0.txt … covariates_M5.txt
  covariate_model_manifest.tsv
"""

from __future__ import annotations

import argparse, sys
from pathlib import Path

import pandas as pd
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import load_paths, norm_region, write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
BASE_COLS = ["agedeath", "sex", "primarydx", "snpPC1", "snpPC2", "snpPC3", "snpPC4", "snpPC5"]
SEED = 20260730


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True, choices=["caudate", "dlpfc", "hippocampus"])
    p.add_argument("--population", default="AA")
    p.add_argument("--peer-k-for-m4", type=int, default=0,
                   help="If >0, also write M4 = M1 + PEER/k (after pilot choice)")
    return p.parse_args()


def median_impute(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    return x.fillna(x.median())


def load_cell_pcs(region: str, sample_ids: list[str], n_pcs: int = 3) -> pd.DataFrame | None:
    path = PROJECT / "inputs" / "cell_proportions" / "_m" / f"music-proportions-{region}.tsv"
    if not path.exists():
        return None
    long = pd.read_csv(path, sep="\t")
    wide = long.pivot_table(index="sample_id", columns="cell_type", values="proportion", aggfunc="first")
    wide.index = wide.index.astype(str)
    wide = wide.reindex(sample_ids)
    if wide.isna().all(axis=None):
        return None
    wide = wide.fillna(wide.median(numeric_only=True))
    # Drop one cell type to reduce composition collinearity, then PCA
    if wide.shape[1] > 1:
        wide = wide.drop(columns=wide.columns[-1])
    mat = wide.to_numpy(dtype=float)
    n_comp = min(n_pcs, mat.shape[1], mat.shape[0] - 1)
    if n_comp < 1:
        return None
    pcs = PCA(n_components=n_comp, random_state=SEED).fit_transform(mat)
    return pd.DataFrame(pcs, index=sample_ids, columns=[f"cellPC{i}" for i in range(1, n_comp + 1)])


def write_cov(df: pd.DataFrame, path: Path) -> None:
    out = df.copy()
    out.index.name = "id"
    if out.isnull().any().any():
        bad = out.columns[out.isnull().any()].tolist()
        raise SystemExit(f"Missing values in {path.name}: {bad}")
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t")


def main() -> None:
    args = parse_args()
    paths = load_paths()
    region = args.region
    base = PROJECT / "meqtl-validation" / "01_cpg_meqtl_mapping" / region / "_m"
    prep = base / "prepared"
    outdir = base / "covariate_sensitivity"
    outdir.mkdir(parents=True, exist_ok=True)
    latent_dir = outdir / "latent"

    inclusion = pd.read_csv(base / "preflight" / "sample_inclusion_primary.tsv", sep="\t")
    keep = inclusion["brnum"].astype(str).tolist()

    phen = pd.read_csv(PROJECT / paths["phenotype_table"], sep="\t")
    phen["brnum"] = phen["brnum"].astype(str)
    phen["_region"] = phen["region"].map(norm_region)
    phen = phen[(phen["race"] == args.population) & (phen["_region"] == region)]
    phen = phen.set_index("brnum").loc[keep]

    # Encode like primary prepare step
    df = phen.copy()
    if df["sex"].dtype == object:
        df["sex"] = df["sex"].astype(str).str.upper().map({"M": 1, "MALE": 1, "F": 0, "FEMALE": 0}).astype(float)
    if df["primarydx"].dtype == object:
        df["primarydx"] = df["primarydx"].map({"Control": 0, "Schizo": 1}).astype(float)

    m0 = df[BASE_COLS].apply(pd.to_numeric, errors="coerce")
    # Prefer existing prepared covariates for exact M0 match
    existing = prep / "covariates.txt"
    if existing.exists():
        m0 = pd.read_csv(existing, sep="\t", index_col=0)
        m0.index = m0.index.astype(str)
        m0 = m0.loc[keep]

    write_cov(m0, outdir / "covariates_M0.txt")

    # M1: + pmi + ph (median impute if needed)
    m1 = m0.copy()
    m1["pmi"] = median_impute(df["pmi"])
    m1["ph"] = median_impute(df["ph"])
    write_cov(m1, outdir / "covariates_M1.txt")

    # M2: + snpPC6–10
    m2 = m0.copy()
    for c in ["snpPC6", "snpPC7", "snpPC8", "snpPC9", "snpPC10"]:
        m2[c] = pd.to_numeric(df[c], errors="coerce")
    if m2.isnull().any().any():
        for c in ["snpPC6", "snpPC7", "snpPC8", "snpPC9", "snpPC10"]:
            m2[c] = median_impute(m2[c])
    write_cov(m2, outdir / "covariates_M2.txt")

    # M3a–c: M0 + latent k
    manifest = []
    for k, mid in [(5, "M3a"), (10, "M3b"), (15, "M3c")]:
        lat_path = latent_dir / f"latent_factors_k{k}.tsv"
        if not lat_path.exists():
            raise SystemExit(f"Missing {lat_path}; run 09_estimate_latent_factors.py first")
        lat = pd.read_csv(lat_path, sep="\t", index_col=0)
        lat.index = lat.index.astype(str)
        mx = m0.join(lat.loc[keep], how="left")
        write_cov(mx, outdir / f"covariates_{mid}.txt")
        manifest.append({"model_id": mid, "n_covariates": mx.shape[1], "description": f"M0 + methPC1–{k}"})

    # M5: cell composition PCs
    cell = load_cell_pcs(region, keep, n_pcs=3)
    if cell is not None:
        m5 = m0.join(cell, how="left")
        write_cov(m5, outdir / "covariates_M5.txt")
        manifest.append({"model_id": "M5", "n_covariates": m5.shape[1], "description": "M0 + cellPC1–3"})
    else:
        print("WARNING: cell proportions unavailable; skipping M5")

    # Optional M4 after pilot
    if args.peer_k_for_m4 in (5, 10, 15):
        lat = pd.read_csv(latent_dir / f"latent_factors_k{args.peer_k_for_m4}.tsv", sep="\t", index_col=0)
        lat.index = lat.index.astype(str)
        m4 = m1.join(lat.loc[keep], how="left")
        write_cov(m4, outdir / "covariates_M4.txt")
        manifest.append({
            "model_id": "M4",
            "n_covariates": m4.shape[1],
            "description": f"M1 + methPC1–{args.peer_k_for_m4}",
        })

    manifest = [
        {"model_id": "M0", "n_covariates": m0.shape[1], "description": "baseline age+sex+dx+snpPC1–5"},
        {"model_id": "M1", "n_covariates": m1.shape[1], "description": "M0 + pmi + ph"},
        {"model_id": "M2", "n_covariates": m2.shape[1], "description": "M0 + snpPC6–10"},
    ] + manifest
    for row in manifest:
        row["region"] = region
        row["n_samples"] = len(keep)
        row["path"] = str(outdir / f"covariates_{row['model_id']}.txt")
    write_tsv(outdir / "covariate_model_manifest.tsv", manifest)
    print(f"Wrote covariate models for {region} under {outdir}")


if __name__ == "__main__":
    main()
