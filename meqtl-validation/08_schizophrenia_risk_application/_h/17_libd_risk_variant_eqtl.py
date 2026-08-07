#!/usr/bin/env python3
"""Level 3: SCZ index SNPs × linked genes in regenerated LIBD AA eQTL cohort.

Separate FDR family within Level-3 SNP–gene tests (not genome-wide cis FDR).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
LEVEL3 = ROOT / "meqtl-validation/08_schizophrenia_risk_application/_m/level3"
LIBD_EQTL = ROOT / "meqtl-validation/09_libd_eqtl_mapping/_m/caudate/genes"
# Fallback during transition if outputs still live under Phase 7
LIBD_EQTL_LEGACY = (
    ROOT
    / "meqtl-validation/08_schizophrenia_risk_application/_m/level3/libd_eqtl/caudate/genes"
)
PRIOR = (
    ROOT
    / "meqtl-validation/08_schizophrenia_risk_application/_m/prioritized/prioritized_loci.tsv"
)
INDEX = (
    ROOT
    / "meqtl-validation/08_schizophrenia_risk_application/_m/caudate/scz_index_snps_hg38.tsv"
)
MEQTL = (
    ROOT
    / "meqtl-validation/08_schizophrenia_risk_application/_m/caudate/risk_variant_cpg_meqtl.tsv.gz"
)


def resolve_libd_genes_dir() -> Path:
    if (LIBD_EQTL / "prepared/normalized_expression.tsv.gz").exists():
        return LIBD_EQTL
    if (LIBD_EQTL_LEGACY / "prepared/normalized_expression.tsv.gz").exists():
        return LIBD_EQTL_LEGACY
    return LIBD_EQTL


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


def load_covariates(path: Path, sample_order: list[str]) -> pd.DataFrame:
    cov = pd.read_csv(path, sep="\t", index_col=0)
    sample_set = set(sample_order)
    if set(cov.index.astype(str)).intersection(sample_set):
        cov.index = cov.index.astype(str)
    elif set(cov.columns.astype(str)).intersection(sample_set):
        cov = cov.T
        cov.index = cov.index.astype(str)
    else:
        raise SystemExit("No covariate/phenotype sample overlap")
    cov = cov.loc[[s for s in sample_order if s in cov.index]]
    return cov.apply(pd.to_numeric, errors="coerce")


def fit_snp_gene(y: np.ndarray, g: np.ndarray, cov: np.ndarray) -> dict:
    n = len(y)
    X = np.column_stack([np.ones(n), g, cov])
    ok = np.isfinite(X).all(axis=1) & np.isfinite(y) & np.isfinite(g)
    X, y = X[ok], y[ok]
    if X.shape[0] < X.shape[1] + 5:
        return {
            "n": int(X.shape[0]),
            "beta": np.nan,
            "se": np.nan,
            "tstat": np.nan,
            "pval": np.nan,
            "r2": np.nan,
        }
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    df = X.shape[0] - rank
    if df <= 0:
        return {
            "n": int(X.shape[0]),
            "beta": np.nan,
            "se": np.nan,
            "tstat": np.nan,
            "pval": np.nan,
            "r2": np.nan,
        }
    sigma2 = float(np.sum(resid**2) / df)
    try:
        xtx_inv = np.linalg.inv(X.T @ X)
        se = float(np.sqrt(max(sigma2 * xtx_inv[1, 1], 0.0)))
    except np.linalg.LinAlgError:
        return {
            "n": int(X.shape[0]),
            "beta": float(beta[1]),
            "se": np.nan,
            "tstat": np.nan,
            "pval": np.nan,
            "r2": np.nan,
        }
    tstat = float(beta[1] / se) if se > 0 else np.nan
    pval = float(2 * stats.t.sf(abs(tstat), df)) if np.isfinite(tstat) else np.nan
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid**2)) / ss_tot if ss_tot > 0 else np.nan
    return {
        "n": int(X.shape[0]),
        "beta": float(beta[1]),
        "se": se,
        "tstat": tstat,
        "pval": pval,
        "r2": r2,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    libd_default = resolve_libd_genes_dir()
    ap.add_argument("--prepared-dir", default=str(libd_default / "prepared"))
    ap.add_argument("--covariates", default=str(libd_default / "standard/covariates.txt"))
    ap.add_argument(
        "--genotype-prefix",
        default=str(ROOT / "inputs/genotypes/TOPMed_LIBD.AA"),
    )
    ap.add_argument(
        "--cis-qtl",
        default=str(libd_default / "tensorqtl/libd_aa_caudate_genes_standard.cis_qtl.txt.gz"),
    )
    ap.add_argument("--outdir", default=str(LEVEL3))
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--maf", type=float, default=0.01)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prepared = Path(args.prepared_dir)

    from tensorqtl import pgen

    targets = pd.read_csv(LEVEL3 / "level3_gene_targets.tsv", sep="\t")
    prior = pd.read_csv(PRIOR, sep="\t")
    idx = pd.read_csv(INDEX, sep="\t")
    idx = idx[idx["index_snp"].isin(set(prior["index_snp"]))].copy()

    expr = pd.read_csv(prepared / "normalized_expression.tsv.gz", sep="\t").set_index(
        "feature_id"
    )
    nov_to_feat: dict[str, str] = {}
    for fid in expr.index.astype(str):
        nov_to_feat.setdefault(fid.split(".")[0], fid)

    pair_rows = []
    for _, t in targets.iterrows():
        gid = str(t.get("gene_id", "") or "")
        if not gid or gid == "nan":
            continue
        fid = nov_to_feat.get(gid.split(".")[0])
        if fid is None:
            continue
        pair_rows.append(
            {
                "locus_id": t["locus_id"],
                "index_snp": t["index_snp"],
                "gene_id": gid.split(".")[0],
                "gene_symbol": t.get("gene_symbol", ""),
                "feature_id": fid,
                "modality_source": t.get("modality", ""),
            }
        )
    pairs = pd.DataFrame(pair_rows).drop_duplicates(["index_snp", "feature_id"])
    if pairs.empty:
        raise SystemExit("No gene targets overlap expression matrix")

    probe = pgen.PgenReader(args.genotype_prefix)
    geno_ids = set(map(str, probe.sample_ids))
    pheno_ids = [c for c in expr.columns.astype(str) if c in geno_ids]
    cov = load_covariates(Path(args.covariates), pheno_ids)
    pheno_ids = [s for s in pheno_ids if s in cov.index]
    pgr = pgen.PgenReader(args.genotype_prefix, select_samples=pheno_ids)
    sample_ids = list(map(str, pgr.sample_ids))
    expr = expr[sample_ids]
    cov = cov.loc[sample_ids]
    if cov.isnull().any().any():
        cov = cov.dropna(axis=1, how="any")
    cov_mat = cov.to_numpy(dtype=float)
    panel_ids = set(map(str, pgr.variant_ids))

    # Resolve genotype variant IDs
    if "genotype_variant_id" not in idx.columns or idx["genotype_variant_id"].isna().all():
        pvar = pgr.pvar_df.copy()
        pvar["id"] = pvar["id"].astype(str)
        pvar["rsid"] = pvar["id"].str.extract(r"(rs\d+)$", expand=False)
        id_by_rs = (
            pvar.dropna(subset=["rsid"])
            .drop_duplicates("rsid")
            .set_index("rsid")["id"]
            .to_dict()
        )
        idx["genotype_variant_id"] = idx["index_snp"].map(id_by_rs)

    meq = pd.read_csv(MEQTL, sep="\t")
    meq_sig = meq[meq["significant_fdr"] == True].copy()  # noqa: E712
    meq_best = (
        meq_sig.sort_values("pval_nominal")
        .groupby("index_snp", as_index=False)
        .first()[
            [
                "index_snp",
                "beta",
                "effect_allele",
                "risk_allele",
                "qval",
                "vmr_id",
            ]
        ]
        .rename(
            columns={
                "beta": "meqtl_beta",
                "effect_allele": "meqtl_effect_allele",
                "risk_allele": "meqtl_risk_allele",
                "qval": "meqtl_qval",
                "vmr_id": "meqtl_vmr_id",
            }
        )
    )

    results = []
    for _, pair in pairs.iterrows():
        snp_row = idx[idx["index_snp"] == pair["index_snp"]]
        if snp_row.empty:
            results.append({**pair.to_dict(), "status": "index_missing", "pval": np.nan})
            continue
        snp_row = snp_row.iloc[0]
        vid = snp_row.get("genotype_variant_id")
        if pd.isna(vid) or str(vid) not in panel_ids:
            results.append(
                {
                    **pair.to_dict(),
                    "variant_id": vid,
                    "status": "variant_missing",
                    "pval": np.nan,
                }
            )
            continue
        try:
            dos = np.asarray(pgr.read_dosages(str(vid)), dtype=float).ravel()
        except Exception as exc:  # noqa: BLE001
            results.append(
                {
                    **pair.to_dict(),
                    "variant_id": vid,
                    "status": f"dosage_failed:{exc}",
                    "pval": np.nan,
                }
            )
            continue
        if dos.shape[0] != len(sample_ids):
            results.append(
                {
                    **pair.to_dict(),
                    "variant_id": vid,
                    "status": "dosage_length_mismatch",
                    "pval": np.nan,
                }
            )
            continue
        if np.isnan(dos).any():
            dos = np.where(np.isnan(dos), np.nanmean(dos), dos)
        maf = float(np.nanmean(dos) / 2.0)
        maf = min(maf, 1.0 - maf)
        if not np.isfinite(maf) or maf < args.maf:
            results.append(
                {
                    **pair.to_dict(),
                    "variant_id": vid,
                    "status": "low_maf",
                    "maf": maf,
                    "pval": np.nan,
                }
            )
            continue
        y = expr.loc[pair["feature_id"]].to_numpy(dtype=float)
        fit = fit_snp_gene(y, dos, cov_mat)
        results.append(
            {
                **pair.to_dict(),
                "variant_id": vid,
                "status": "tested",
                "chrom": snp_row.get("chrom"),
                "pos": snp_row.get("pos_hg38"),
                "ref": snp_row.get("REF"),
                "alt": snp_row.get("ALT"),
                "risk_allele": snp_row.get("risk_allele"),
                "gwas_or": snp_row.get("gwas_or"),
                "maf": maf,
                "effect_allele": snp_row.get("ALT"),
                "n": fit["n"],
                "eqtl_beta": fit["beta"],
                "eqtl_se": fit["se"],
                "eqtl_tstat": fit["tstat"],
                "eqtl_pval": fit["pval"],
                "eqtl_r2": fit["r2"],
            }
        )

    res = pd.DataFrame(results)
    res["qval"] = bh_fdr(res["eqtl_pval"].to_numpy()) if "eqtl_pval" in res.columns else np.nan
    res["significant_fdr"] = res["qval"] <= args.fdr
    res = res.merge(meq_best, on="index_snp", how="left")

    if "eqtl_beta" in res.columns:
        flip_e = res["effect_allele"].astype(str) != res["risk_allele"].astype(str)
        res["eqtl_beta_risk"] = res["eqtl_beta"].where(~flip_e, -res["eqtl_beta"])
        flip_m = res["meqtl_effect_allele"].astype(str) != res["risk_allele"].astype(str)
        res["meqtl_beta_risk"] = res["meqtl_beta"].where(~flip_m, -res["meqtl_beta"])
        res["same_sign_meqtl_eqtl_risk"] = np.sign(res["eqtl_beta_risk"]) == np.sign(
            res["meqtl_beta_risk"]
        )

    if args.cis_qtl and Path(args.cis_qtl).exists():
        cis = pd.read_csv(args.cis_qtl, sep="\t")
        if cis.columns[0] != "phenotype_id" and "phenotype_id" not in cis.columns:
            cis = cis.rename(columns={cis.columns[0]: "phenotype_id"})
        cis = cis.set_index("phenotype_id")
        res["cis_qval"] = res["feature_id"].map(cis["qval"]) if "qval" in cis.columns else np.nan
        res["cis_egene_fdr"] = res["cis_qval"] <= args.fdr

    res.to_csv(outdir / "libd_risk_variant_eqtl.tsv.gz", sep="\t", index=False)
    snap = pd.DataFrame(
        [
            {
                "n_pairs_tested": int((res["status"] == "tested").sum()),
                "n_sig_fdr": int(res["significant_fdr"].fillna(False).sum()),
                "n_loci_with_sig_eqtl": int(
                    res.loc[res["significant_fdr"].fillna(False), "locus_id"].nunique()
                ),
                "n_pairs_same_sign_meqtl": int(
                    res.loc[
                        res["significant_fdr"].fillna(False)
                        & res.get("same_sign_meqtl_eqtl_risk", False).fillna(False),
                    ].shape[0]
                )
                if "same_sign_meqtl_eqtl_risk" in res.columns
                else 0,
                "n_samples": len(sample_ids),
                "fdr": args.fdr,
                "maf": args.maf,
            }
        ]
    )
    snap.to_csv(outdir / "libd_risk_variant_eqtl_summary.tsv", sep="\t", index=False)
    print(snap.to_string(index=False))
    print(f"Wrote {outdir}")


if __name__ == "__main__":
    main()
