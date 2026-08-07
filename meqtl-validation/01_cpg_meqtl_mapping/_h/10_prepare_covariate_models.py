#!/usr/bin/env python3
"""
Build covariate matrices for Phase 1 sensitivity models M0–M6d.

Writes under {region}/_m/covariate_sensitivity/:
  covariates_M0.txt … covariates_M6d.txt (when its prespecified gate passes)
  covariate_model_manifest.tsv
"""

from __future__ import annotations

import argparse, sys
from pathlib import Path

import numpy as np
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


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def load_dnam_cell_pcs(
    region: str,
    sample_ids: list[str],
    outdir: Path,
    n_pcs: int = 3,
) -> tuple[pd.DataFrame | None, str]:
    """Create CLR-PCA covariates only after the locked cross-modality QC gate."""
    validation_path = (
        PROJECT / "inputs" / "cell_proportions" / "_m"
        / f"dnam-scmd-validation-{region}.tsv"
    )
    prop_path = (
        PROJECT / "inputs" / "cell_proportions" / "_m"
        / f"dnam-scmd-proportions-{region}.tsv"
    )
    if not validation_path.exists() or not prop_path.exists():
        return None, "DNAm proportion or validation file is missing"

    validation = pd.read_csv(validation_path, sep="\t")
    if len(validation) != 1 or not truthy(validation.iloc[0].get("integration_pass", False)):
        return None, "region did not pass the prespecified DNAm integration gate"

    long = pd.read_csv(prop_path, sep="\t")
    required = {"sample_id", "cell_type", "proportion"}
    missing = required.difference(long.columns)
    if missing:
        raise SystemExit(f"Missing columns in {prop_path}: {sorted(missing)}")
    wide = long.pivot_table(
        index="sample_id", columns="cell_type", values="proportion", aggfunc="first"
    )
    wide.index = wide.index.astype(str)
    retained_ids = [sample_id for sample_id in sample_ids if sample_id in wide.index]
    absent = [sample_id for sample_id in sample_ids if sample_id not in wide.index]
    if not retained_ids:
        return None, "no model samples have DNAm proportions"
    wide = wide.reindex(retained_ids).apply(pd.to_numeric, errors="coerce")
    if wide.isna().any().any():
        return None, "DNAm fractions contain missing values for model samples"
    if (wide < 0).any().any() or (wide.sum(axis=1) <= 0).any():
        return None, "DNAm fractions are negative or have a non-positive row sum"

    positive = wide.to_numpy(dtype=float)
    positive_values = positive[positive > 0]
    if positive_values.size == 0:
        return None, "DNAm composition matrix has no positive fractions"
    zero_replacement = float(positive_values.min() / 2.0)
    positive[positive <= 0] = zero_replacement
    closed = positive / positive.sum(axis=1, keepdims=True)
    logged = np.log(closed)
    clr = logged - logged.mean(axis=1, keepdims=True)

    n_comp = min(n_pcs, clr.shape[1] - 1, clr.shape[0] - 1)
    if n_comp < 1:
        return None, "too few samples or cell classes for DNAm CLR-PCA"
    pca = PCA(n_components=n_comp, random_state=SEED)
    scores = pca.fit_transform(clr)
    columns = [f"dnamCellPC{i}" for i in range(1, n_comp + 1)]
    pcs = pd.DataFrame(scores, index=retained_ids, columns=columns)

    diagnostics = outdir / "dnam_cell_pca"
    diagnostics.mkdir(parents=True, exist_ok=True)
    clr_table = pd.DataFrame(clr, index=retained_ids, columns=wide.columns)
    clr_table.index.name = "id"
    clr_table.to_csv(diagnostics / "dnam_cell_clr.tsv", sep="\t")
    loadings = pd.DataFrame(pca.components_.T, index=wide.columns, columns=columns)
    loadings.index.name = "cell_type"
    loadings.to_csv(diagnostics / "dnam_cell_pca_loadings.tsv", sep="\t")
    pd.DataFrame({
        "component": columns,
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "zero_replacement": zero_replacement,
    }).to_csv(diagnostics / "dnam_cell_pca_summary.tsv", sep="\t", index=False)
    return pcs, f"passed; excluded {len(absent)} failed samples without imputing fractions"


def add_rank_safe_pcs(base: pd.DataFrame, pcs: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Drop highest PCs only when required to retain a full-rank design."""
    numeric_base = base.apply(pd.to_numeric, errors="coerce")
    if numeric_base.isna().any().any():
        raise SystemExit("Cannot assess M6d rank because M3a contains non-numeric or missing values")
    for n_keep in range(pcs.shape[1], 0, -1):
        candidate = base.join(pcs.iloc[:, :n_keep], how="left")
        design = np.column_stack([
            np.ones(candidate.shape[0]),
            candidate.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float),
        ])
        if np.isfinite(design).all() and np.linalg.matrix_rank(design) == design.shape[1]:
            return candidate, n_keep
    raise SystemExit("M6d is rank deficient even after reducing DNAm composition to one PC")


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
    if args.population == "AA":
        prep = base / "prepared"
        outdir = base / "covariate_sensitivity"
        latent_dir = outdir / "latent"
        preflight = base / "preflight" / "sample_inclusion_primary.tsv"
    else:
        prep = base / "prepared" / args.population
        outdir = prep  # write EA model matrices next to EA phenotypes
        latent_dir = prep / "latent"
        preflight = base / "preflight" / args.population / "sample_inclusion_primary.tsv"
    outdir.mkdir(parents=True, exist_ok=True)

    inclusion = pd.read_csv(preflight, sep="\t")
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
    # Prefer existing prepared covariates for exact M0 match (strip methPCs if present)
    existing = prep / "covariates.txt"
    if existing.exists():
        m0 = pd.read_csv(existing, sep="\t", index_col=0)
        m0.index = m0.index.astype(str)
        m0 = m0.loc[keep]
        meth_cols = [c for c in m0.columns if str(c).startswith("methPC")]
        if meth_cols:
            m0 = m0.drop(columns=meth_cols)

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
    m3a = None
    for k, mid in [(5, "M3a"), (10, "M3b"), (15, "M3c")]:
        lat_path = latent_dir / f"latent_factors_k{k}.tsv"
        if not lat_path.exists():
            raise SystemExit(f"Missing {lat_path}; run 09_estimate_latent_factors.py first")
        lat = pd.read_csv(lat_path, sep="\t", index_col=0)
        lat.index = lat.index.astype(str)
        mx = m0.join(lat.loc[keep], how="left")
        write_cov(mx, outdir / f"covariates_{mid}.txt")
        manifest.append({"model_id": mid, "n_covariates": mx.shape[1], "description": f"M0 + methPC1–{k}"})
        if mid == "M3a":
            m3a = mx

    # M5: cell composition PCs
    cell = load_cell_pcs(region, keep, n_pcs=3)
    if cell is not None:
        m5 = m0.join(cell, how="left")
        write_cov(m5, outdir / "covariates_M5.txt")
        manifest.append({"model_id": "M5", "n_covariates": m5.shape[1], "description": "M0 + legacy RNA MuSiC cellPC1–3"})
    else:
        print("WARNING: cell proportions unavailable; skipping M5")

    # M6d: locked M3a + DNAm scMD CLR-PCA, only for regions passing integration QC.
    if args.population == "AA":
        dnam_pcs, dnam_status = load_dnam_cell_pcs(region, keep, outdir, n_pcs=3)
        if dnam_pcs is not None:
            if m3a is None:
                raise SystemExit("Internal error: M3a was not built before M6d")
            matched_m3a = m3a.loc[dnam_pcs.index]
            m6d, n_dnam_pcs = add_rank_safe_pcs(matched_m3a, dnam_pcs)
            write_cov(matched_m3a, outdir / "covariates_M3a_dnam_matched.txt")
            write_cov(m6d, outdir / "covariates_M6d.txt")
            sample_file = outdir / "dnam_cell_pca" / "m6d_sample_ids.txt"
            sample_file.write_text("\n".join(m6d.index.astype(str)) + "\n")
            manifest.append({
                "model_id": "M3a_dnam_matched",
                "n_covariates": matched_m3a.shape[1],
                "description": "M3a restricted to samples with QC-passing DNAm fractions",
                "n_samples_override": matched_m3a.shape[0],
            })
            manifest.append({
                "model_id": "M6d",
                "n_covariates": m6d.shape[1],
                "description": (
                    f"M3a + dnamCellPC1–{n_dnam_pcs} "
                    f"(scMD fractions; zero-replaced, closed, CLR-PCA; QC gated; {dnam_status})"
                ),
                "n_samples_override": m6d.shape[0],
            })
        else:
            print(f"WARNING: skipping M6d for {region}: {dnam_status}")

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
        row["n_samples"] = row.pop("n_samples_override", len(keep))
        row["path"] = str(outdir / f"covariates_{row['model_id']}.txt")
    write_tsv(outdir / "covariate_model_manifest.tsv", manifest)
    print(f"Wrote covariate models for {region} under {outdir}")


if __name__ == "__main__":
    main()
