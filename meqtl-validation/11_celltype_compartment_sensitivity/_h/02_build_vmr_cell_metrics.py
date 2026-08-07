#!/usr/bin/env python3
"""Build VMR mean-methylation matrices and cell-composition metrics.

For each VMR × AA meQTL sample:
  - mean CpG methylation within the VMR
Per VMR summaries:
  - mean / variance of methylation
  - R² of meth ~ MuSiC cellPC1–3
  - R² of meth ~ Oligo fraction
  - corr(meth, Oligo)
  - residual variance after cellPC adjustment
  - (caudate) R² of meth ~ dnamCellPC1–3
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
MODULE = PROJECT / "meqtl-validation/11_celltype_compartment_sensitivity/_m"
PHASE6 = PROJECT / "meqtl-validation/07_repeat_mappability_sensitivity/_m"
SEED = 20260807


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions", nargs="+", default=["caudate", "dlpfc", "hippocampus"])
    return p.parse_args()


def coord_id(chrom, start, end) -> str:
    c = str(chrom).replace("chr", "")
    return f"{c}:{int(start)}-{int(end)}"


def load_sample_ids(region: str) -> list[str]:
    cov = pd.read_csv(
        PROJECT / f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/prepared/covariates.txt",
        sep="\t",
    )
    return cov["id"].astype(str).tolist()


def load_music_wide(region: str, sample_ids: list[str]) -> pd.DataFrame:
    path = PROJECT / f"inputs/cell_proportions/_m/music-proportions-{region}.tsv"
    long = pd.read_csv(path, sep="\t")
    wide = long.pivot_table(index="sample_id", columns="cell_type", values="proportion", aggfunc="first")
    wide.index = wide.index.astype(str)
    wide = wide.reindex(sample_ids)
    return wide.apply(pd.to_numeric, errors="coerce")


def load_cell_pcs_from_m5(region: str, sample_ids: list[str]) -> pd.DataFrame:
    path = (
        PROJECT
        / f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/covariate_sensitivity/covariates_M5.txt"
    )
    m5 = pd.read_csv(path, sep="\t")
    m5["id"] = m5["id"].astype(str)
    m5 = m5.set_index("id")
    cols = [c for c in ["cellPC1", "cellPC2", "cellPC3"] if c in m5.columns]
    return m5.reindex(sample_ids)[cols].apply(pd.to_numeric, errors="coerce")


def load_dnam_cell_pcs(region: str, sample_ids: list[str], outdir: Path) -> tuple[pd.DataFrame | None, str]:
    m6d = (
        PROJECT
        / f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/covariate_sensitivity/covariates_M6d.txt"
    )
    if m6d.exists():
        d = pd.read_csv(m6d, sep="\t")
        d["id"] = d["id"].astype(str)
        d = d.set_index("id")
        cols = [c for c in d.columns if c.startswith("dnamCellPC")]
        pcs = d.reindex([s for s in sample_ids if s in d.index])[cols].apply(pd.to_numeric, errors="coerce")
        return pcs, "from_M6d"

    # Exploratory CLR-PCA even if integration gate failed (flagged in status)
    prop = PROJECT / f"inputs/cell_proportions/_m/dnam-scmd-proportions-{region}.tsv"
    val = PROJECT / f"inputs/cell_proportions/_m/dnam-scmd-validation-{region}.tsv"
    gate = False
    if val.exists():
        v = pd.read_csv(val, sep="\t")
        gate = len(v) == 1 and str(v.iloc[0].get("integration_pass", "")).lower() in {
            "true", "1", "yes", "t"
        }
    if not prop.exists():
        return None, "missing_dnam_proportions"
    long = pd.read_csv(prop, sep="\t")
    wide = long.pivot_table(index="sample_id", columns="cell_type", values="proportion", aggfunc="first")
    wide.index = wide.index.astype(str)
    keep = [s for s in sample_ids if s in wide.index]
    if len(keep) < 20:
        return None, "too_few_dnam_samples"
    wide = wide.reindex(keep).apply(pd.to_numeric, errors="coerce")
    if wide.isna().any().any():
        return None, "dnam_missing_values"
    positive = wide.to_numpy(dtype=float)
    pos = positive[positive > 0]
    zero_rep = float(pos.min() / 2.0)
    positive[positive <= 0] = zero_rep
    closed = positive / positive.sum(axis=1, keepdims=True)
    logged = np.log(closed)
    clr = logged - logged.mean(axis=1, keepdims=True)
    n_comp = min(3, clr.shape[1] - 1, clr.shape[0] - 1)
    scores = PCA(n_components=n_comp, random_state=SEED).fit_transform(clr)
    cols = [f"dnamCellPC{i}" for i in range(1, n_comp + 1)]
    pcs = pd.DataFrame(scores, index=keep, columns=cols)
    tag = "exploratory_clr_pca" if not gate else "clr_pca_gate_pass"
    pcs.to_csv(outdir / f"dnam_cell_pcs_{tag}.tsv", sep="\t")
    return pcs, tag


def build_cpg_to_vmr(region: str, keep_keys: set[str]) -> dict[str, str]:
    """Map phenotype_id -> vmr key (task_id preferred, else coord)."""
    prep = PROJECT / f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/prepared"
    out: dict[str, str] = {}
    for mpath in sorted(prep.glob("cpg_vmr_map.chr*.tsv")):
        m = pd.read_csv(mpath, sep="\t", dtype={"vmr_id": str, "phenotype_id": str})
        for pid, vid in zip(m["phenotype_id"], m["vmr_id"]):
            if vid in keep_keys:
                out[pid] = vid
    return out


def load_vmr_keys(region: str) -> tuple[pd.DataFrame, set[str]]:
    tech = pd.read_csv(PHASE6 / region / "vmr_technical_annotations.tsv", sep="\t")
    tech["coord_id"] = [
        coord_id(c, s, e) for c, s, e in zip(tech["chrom"], tech["start"], tech["end"])
    ]
    keys = set(tech["coord_id"].astype(str))
    if "task_id" in tech.columns:
        tids = tech["task_id"].dropna()
        keys |= set(tids.astype(float).astype(int).astype(str))
    if "interval_id" in tech.columns:
        keys |= set(tech["interval_id"].dropna().astype(str))
    return tech, keys


def accumulate_vmr_means(
    region: str, sample_ids: list[str], cpg_to_vmr: dict[str, str]
) -> pd.DataFrame:
    bed = (
        PROJECT
        / f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/prepared/"
        "cpg_phenotypes.all_autosomes.bed.gz"
    )
    sample_index = {s: i for i, s in enumerate(sample_ids)}
    n = len(sample_ids)
    sums: dict[str, np.ndarray] = {}
    counts: dict[str, np.ndarray] = {}

    with gzip.open(bed, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        bed_samples = header[4:]
        bed_pos_to_si = [
            (i, sample_index[s]) for i, s in enumerate(bed_samples) if s in sample_index
        ]
        if not bed_pos_to_si:
            raise SystemExit(f"No overlapping samples between covariates and {bed}")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            pid = parts[3]
            vid = cpg_to_vmr.get(pid)
            if vid is None:
                continue
            if vid not in sums:
                sums[vid] = np.zeros(n, dtype=np.float64)
                counts[vid] = np.zeros(n, dtype=np.float64)
            for bi, si in bed_pos_to_si:
                v = parts[4 + bi]
                if v in ("", "NA", "NaN"):
                    continue
                sums[vid][si] += float(v)
                counts[vid][si] += 1.0

    if not sums:
        raise SystemExit(f"No VMR methylation accumulated for {region}")
    mat = {}
    for vid, s in sums.items():
        with np.errstate(invalid="ignore", divide="ignore"):
            mat[vid] = np.where(counts[vid] > 0, s / counts[vid], np.nan)
    return pd.DataFrame(mat, index=sample_ids)


def r2_ols(y: np.ndarray, X: np.ndarray) -> float:
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if ok.sum() < 20:
        return np.nan
    yy = y[ok]
    XX = np.column_stack([np.ones(ok.sum()), X[ok]])
    try:
        beta, _, _, _ = np.linalg.lstsq(XX, yy, rcond=None)
    except np.linalg.LinAlgError:
        return np.nan
    resid = yy - XX @ beta
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    if ss_tot <= 0:
        return np.nan
    return 1.0 - ss_res / ss_tot


def residual_var(y: np.ndarray, X: np.ndarray) -> float:
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if ok.sum() < 20:
        return np.nan
    yy = y[ok]
    XX = np.column_stack([np.ones(ok.sum()), X[ok]])
    beta, _, _, _ = np.linalg.lstsq(XX, yy, rcond=None)
    resid = yy - XX @ beta
    return float(np.var(resid, ddof=0))


def summarize_vmrs(
    meth: pd.DataFrame,
    cell_pcs: pd.DataFrame,
    oligo: pd.Series,
    dnam_pcs: pd.DataFrame | None,
    tech: pd.DataFrame,
) -> pd.DataFrame:
    # Map keys back to coord_id
    key_to_coord: dict[str, str] = {}
    for _, r in tech.iterrows():
        cid = str(r["coord_id"])
        key_to_coord[cid] = cid
        if pd.notna(r.get("task_id")):
            key_to_coord[str(int(float(r["task_id"])))] = cid
        if pd.notna(r.get("interval_id")):
            key_to_coord[str(r["interval_id"])] = cid

    pc_mat = cell_pcs.to_numpy(dtype=float)
    oligo_v = oligo.to_numpy(dtype=float)
    dnam_mat = dnam_pcs.to_numpy(dtype=float) if dnam_pcs is not None and len(dnam_pcs) else None
    dnam_index = set(dnam_pcs.index) if dnam_pcs is not None else set()

    rows = []
    for vid in meth.columns:
        y = meth[vid].to_numpy(dtype=float)
        mean_m = float(np.nanmean(y))
        var_m = float(np.nanvar(y))
        cell_r2 = r2_ols(y, pc_mat)
        oligo_r2 = r2_ols(y, oligo_v.reshape(-1, 1))
        ok = np.isfinite(y) & np.isfinite(oligo_v)
        oligo_corr = float(np.corrcoef(y[ok], oligo_v[ok])[0, 1]) if ok.sum() >= 20 else np.nan
        resid_v = residual_var(y, pc_mat)
        dnam_r2 = np.nan
        if dnam_mat is not None:
            # align to dnam sample subset
            idx = [i for i, s in enumerate(meth.index) if s in dnam_index]
            if len(idx) >= 20:
                y_d = y[idx]
                # dnam_pcs row order by meth.index filtered
                Xd = dnam_pcs.reindex(meth.index[idx]).to_numpy(dtype=float)
                dnam_r2 = r2_ols(y_d, Xd)
        rows.append(
            {
                "vmr_key": vid,
                "coord_id": key_to_coord.get(str(vid), ""),
                "n_samples_meth": int(np.isfinite(y).sum()),
                "mean_methylation": mean_m,
                "var_methylation": var_m,
                "cellPC_r2": cell_r2,
                "oligo_r2": oligo_r2,
                "oligo_corr": oligo_corr,
                "abs_oligo_corr": abs(oligo_corr) if np.isfinite(oligo_corr) else np.nan,
                "residual_var_cellPC": resid_v,
                "dnamCellPC_r2": dnam_r2,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    MODULE.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    for region in args.regions:
        print(f"==== {region} ====")
        outdir = MODULE / region
        outdir.mkdir(parents=True, exist_ok=True)
        sample_ids = load_sample_ids(region)
        tech, keep_keys = load_vmr_keys(region)
        print(f"  samples={len(sample_ids)} keep_keys={len(keep_keys)}")

        cpg_to_vmr = build_cpg_to_vmr(region, keep_keys)
        print(f"  cpg mapped={len(cpg_to_vmr)}")
        meth = accumulate_vmr_means(region, sample_ids, cpg_to_vmr)
        meth_path = outdir / "vmr_mean_methylation.tsv.gz"
        meth.to_csv(meth_path, sep="\t", compression="gzip")
        print(f"  meth matrix {meth.shape} -> {meth_path}")

        cell_pcs = load_cell_pcs_from_m5(region, sample_ids)
        music = load_music_wide(region, sample_ids)
        oligo_col = "Oligo" if "Oligo" in music.columns else None
        if oligo_col is None:
            raise SystemExit(f"No Oligo column in MuSiC for {region}")
        oligo = music[oligo_col]

        # neuron aggregate if available
        neuron_cols = [c for c in music.columns if c in {"D1-SPN", "D2-SPN", "Excit", "Inhib"}]
        neuron = music[neuron_cols].sum(axis=1) if neuron_cols else pd.Series(np.nan, index=sample_ids)

        dnam_pcs, dnam_status = load_dnam_cell_pcs(region, sample_ids, outdir)
        print(f"  dnam pcs: {dnam_status}")

        metrics = summarize_vmrs(meth, cell_pcs, oligo, dnam_pcs, tech)
        # attach tech flags
        tech_small = tech[
            [c for c in ["coord_id", "task_id", "line_l1_frac", "umap_k24_mean", "high_mappability", "overlaps_segdup", "overlaps_snp_prox", "length"] if c in tech.columns]
        ].drop_duplicates("coord_id")
        metrics = metrics.merge(tech_small, on="coord_id", how="left")
        metrics_path = outdir / "vmr_cell_metrics.tsv.gz"
        metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

        # save sample covariates used
        scov = pd.DataFrame(
            {
                "id": sample_ids,
                "Oligo": oligo.reindex(sample_ids).to_numpy(),
                "neuron_agg": neuron.reindex(sample_ids).to_numpy(),
            }
        )
        for c in cell_pcs.columns:
            scov[c] = cell_pcs[c].to_numpy()
        if dnam_pcs is not None:
            for c in dnam_pcs.columns:
                scov[c] = dnam_pcs.reindex(sample_ids)[c].to_numpy()
        scov.to_csv(outdir / "sample_cell_covariates.tsv", sep="\t", index=False)

        summary_rows.append(
            {
                "region": region,
                "n_samples": len(sample_ids),
                "n_vmrs_with_meth": int(meth.shape[1]),
                "n_metrics_with_coord": int(metrics["coord_id"].ne("").sum()),
                "median_cellPC_r2": float(metrics["cellPC_r2"].median(skipna=True)),
                "median_oligo_r2": float(metrics["oligo_r2"].median(skipna=True)),
                "median_abs_oligo_corr": float(metrics["abs_oligo_corr"].median(skipna=True)),
                "median_dnamCellPC_r2": float(metrics["dnamCellPC_r2"].median(skipna=True)),
                "dnam_status": dnam_status,
            }
        )
        print(f"  wrote {metrics_path}")

    write_tsv(MODULE / "vmr_cell_metrics_summary.tsv", summary_rows)
    print(pd.DataFrame(summary_rows).to_string(index=False))


if __name__ == "__main__":
    main()
