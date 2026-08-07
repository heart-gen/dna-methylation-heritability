#!/usr/bin/env python3
"""Analysis 3: targeted SCZ risk-variant × CpG cis-meQTL tests.

Uses a separate FDR family from genome-wide meQTL mapping.
Primary: caudate AA with locked M3a covariates.
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
WINDOW = 500_000
SEED = 20260801


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", default="caudate")
    p.add_argument("--population", default="AA", choices=["AA", "EA"])
    p.add_argument("--window", type=int, default=WINDOW)
    p.add_argument("--maf", type=float, default=0.05)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--covariates", default="")
    p.add_argument("--phenotype-bed", default="")
    p.add_argument("--genotype-prefix", default="")
    p.add_argument("--outdir", default="")
    p.add_argument("--seed", type=int, default=SEED)
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
        raise SystemExit("No overlapping samples between covariates and phenotypes")
    cov = cov.loc[[s for s in sample_order if s in cov.index]]
    return cov.apply(pd.to_numeric, errors="coerce")


def residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Return residuals of y on columns of X (with intercept). y is n×p or n."""
    n = X.shape[0]
    Xd = np.column_stack([np.ones(n), X])
    # Drop constant/collinear columns lightly via least squares
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


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir) if args.outdir else (
        PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m" / args.region
    )
    outdir.mkdir(parents=True, exist_ok=True)

    index = pd.read_csv(outdir / "scz_index_snps_hg38.tsv", sep="\t")
    index = index[index["in_genotype_panel"].fillna(False).astype(bool)].dropna(
        subset=["genotype_variant_id", "chrom", "pos_hg38"]
    ).copy()
    if index.empty:
        raise SystemExit("No panel-matched index SNPs")

    phase1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping" / args.region / "_m"
    if args.population == "AA":
        pheno_bed = Path(args.phenotype_bed) if args.phenotype_bed else (
            phase1 / "prepared/cpg_phenotypes.all_autosomes.bed.gz"
        )
        cov_path = Path(args.covariates) if args.covariates else (
            phase1 / "prepared/covariates.txt"
        )
        geno_prefix = Path(args.genotype_prefix) if args.genotype_prefix else (
            phase1 / "genotypes/meqtl_AA"
        )
        prepared = phase1 / "prepared"
    else:
        pheno_bed = Path(args.phenotype_bed) if args.phenotype_bed else (
            phase1 / "prepared/EA/cpg_phenotypes.all_autosomes.bed.gz"
        )
        cov_path = Path(args.covariates) if args.covariates else (
            phase1 / "prepared/EA/covariates_M3a.txt"
        )
        geno_prefix = Path(args.genotype_prefix) if args.genotype_prefix else (
            phase1 / "genotypes/meqtl_EA"
        )
        prepared = phase1 / "prepared/EA"

    from tensorqtl import pgen  # noqa: WPS433
    import tensorqtl

    phenotype_df, _phenotype_pos_df = tensorqtl.read_phenotype_bed(str(pheno_bed))
    phenotype_df.columns = phenotype_df.columns.astype(str)
    probe = pgen.PgenReader(str(geno_prefix))
    geno_ids = set(map(str, probe.sample_ids))
    samples_pheno = [s for s in phenotype_df.columns if s in geno_ids]
    cov = load_covariates(cov_path, samples_pheno)
    samples_pheno = [s for s in samples_pheno if s in cov.index]

    # PgenReader(select_samples=...) keeps native genotype order, not the
    # request-list order. Align phenotypes/covariates to pgr.sample_ids.
    pgr = pgen.PgenReader(str(geno_prefix), select_samples=samples_pheno)
    samples = list(map(str, pgr.sample_ids))
    phenotype_df = phenotype_df[samples]
    cov = cov.loc[samples]
    if cov.isnull().any().any():
        cov = cov.dropna(axis=1, how="any")
    X = cov.to_numpy(dtype=float)
    panel_ids = set(map(str, pgr.variant_ids))

    maps = sorted(prepared.glob("cpg_vmr_map.chr*.tsv"))
    if not maps:
        maps = sorted((phase1 / "prepared").glob("cpg_vmr_map.chr*.tsv"))
    cpg_map = pd.concat([pd.read_csv(p, sep="\t") for p in maps], ignore_index=True)
    cpg_map["vmr_id"] = cpg_map["vmr_id"].astype(str)
    cpg_map["chrom"] = cpg_map["chrom"].astype(str)
    if not str(cpg_map["chrom"].iloc[0]).startswith("chr"):
        cpg_map["chrom"] = "chr" + cpg_map["chrom"].str.replace("^chr", "", regex=True)
    cpg_map["pos"] = pd.to_numeric(cpg_map["pos_1based"], errors="coerce")

    burden = pd.read_csv(
        PROJECT / "meqtl-validation/02_vmr_meqtl_burden/_m" / args.region / "vmr_meqtl_burden.tsv.gz",
        sep="\t",
    )
    burden["vmr_id"] = burden["vmr_id"].astype(str)
    pred_map = burden.set_index("vmr_id")["local_predictability"].to_dict() if "local_predictability" in burden.columns else {}

    # Residualize all phenotypes once (n_samples × n_cpgs is large — do per chrom)
    results = []
    failed = []

    for chrom, idx_chr in index.groupby("chrom"):
        cpg_chr = cpg_map[cpg_map["chrom"] == chrom]
        if cpg_chr.empty:
            continue
        pheno_ids_chr = set(cpg_chr["phenotype_id"]).intersection(phenotype_df.index)
        if not pheno_ids_chr:
            continue
        Y = phenotype_df.loc[sorted(pheno_ids_chr)].T.to_numpy(dtype=float)  # n × p
        # Impute CpG missing with column mean
        col_mean = np.nanmean(Y, axis=0)
        inds = np.where(np.isnan(Y))
        if inds[0].size:
            Y[inds] = np.take(col_mean, inds[1])
        Y_res = residualize(Y, X)
        pheno_ids_list = sorted(pheno_ids_chr)
        pheno_pos = cpg_chr.set_index("phenotype_id").loc[pheno_ids_list]
        pheno_vmr = pheno_pos["vmr_id"].to_numpy()
        pheno_coord = pheno_pos["pos"].to_numpy(dtype=float)

        for _, r in idx_chr.iterrows():
            vid = str(r["genotype_variant_id"])
            if vid not in panel_ids:
                failed.append({"index_snp": r["index_snp"], "reason": "variant_not_in_pgen_index", "variant_id": vid})
                continue
            try:
                dos = np.asarray(pgr.read_dosages(vid), dtype=float).ravel()
            except Exception as exc:  # noqa: BLE001
                failed.append({"index_snp": r["index_snp"], "reason": f"dosage_read_failed:{exc}", "variant_id": vid})
                continue
            if dos.shape[0] != len(samples):
                failed.append({"index_snp": r["index_snp"], "reason": "dosage_length_mismatch", "variant_id": vid})
                continue
            # MAF on dosages
            maf = float(np.nanmean(dos) / 2.0)
            maf = min(maf, 1 - maf)
            if not np.isfinite(maf) or maf < args.maf:
                failed.append({"index_snp": r["index_snp"], "reason": f"maf_below_threshold:{maf:.4f}", "variant_id": vid, "maf": maf})
                continue
            if np.isnan(dos).any():
                dos = np.where(np.isnan(dos), np.nanmean(dos), dos)
            g_res = residualize(dos, X).ravel()
            if np.std(g_res) < 1e-8:
                failed.append({"index_snp": r["index_snp"], "reason": "zero_genotype_variance", "variant_id": vid})
                continue

            risk_pos = int(r["pos_hg38"])
            in_win = (pheno_coord >= risk_pos - args.window) & (pheno_coord <= risk_pos + args.window)
            if not in_win.any():
                continue
            Yw = Y_res[:, in_win]
            ids_w = np.asarray(pheno_ids_list)[in_win]
            vmr_w = pheno_vmr[in_win]
            pos_w = pheno_coord[in_win]

            # Vectorized correlations
            g = g_res - g_res.mean()
            gss = np.dot(g, g)
            Yc = Yw - Yw.mean(axis=0, keepdims=True)
            cov_gy = g @ Yc
            yss = np.sum(Yc * Yc, axis=0)
            denom = np.sqrt(gss * yss)
            with np.errstate(invalid="ignore", divide="ignore"):
                r_xy = np.where(denom > 0, cov_gy / denom, np.nan)
            # slope = cov(g,y)/var(g)
            slope = np.where(gss > 0, cov_gy / gss, np.nan)
            n = len(samples)
            df = n - 2 - X.shape[1]  # residualize approx df
            if df < 5:
                df = n - 2
            with np.errstate(invalid="ignore", divide="ignore"):
                tstat = r_xy * np.sqrt(df / np.clip(1 - r_xy ** 2, 1e-12, None))
            pval = 2 * stats.t.sf(np.abs(tstat), df)
            # SE of slope
            resid_ss = yss - (cov_gy ** 2) / gss
            slope_se = np.sqrt(np.clip(resid_ss / (df * gss), 0, None))
            r2 = r_xy ** 2

            for j in range(len(ids_w)):
                if not np.isfinite(pval[j]):
                    continue
                results.append({
                    "region": args.region,
                    "population": args.population,
                    "locus_id": r.get("locus_id", pd.NA),
                    "index_snp": r["index_snp"],
                    "variant_id": vid,
                    "risk_allele": r.get("risk_allele", ""),
                    "other_allele": r.get("other_allele", ""),
                    "effect_allele": r.get("ALT", ""),  # dosage counts ALT copies in plink/pgen typically
                    "chrom": chrom,
                    "variant_pos": risk_pos,
                    "phenotype_id": ids_w[j],
                    "cpg_pos": int(pos_w[j]),
                    "cpg_to_variant_distance": int(abs(pos_w[j] - risk_pos)),
                    "vmr_id": str(vmr_w[j]),
                    "local_predictability": pred_map.get(str(vmr_w[j]), np.nan),
                    "n": n,
                    "maf": maf,
                    "beta": float(slope[j]),
                    "se": float(slope_se[j]),
                    "tstat": float(tstat[j]),
                    "pval_nominal": float(pval[j]),
                    "r2": float(r2[j]),
                    "gwas_p": r.get("gwas_p", np.nan),
                    "gwas_or": r.get("gwas_or", np.nan),
                })

    res = pd.DataFrame(results)
    if res.empty:
        write_tsv(outdir / "risk_variant_cpg_tests_summary.tsv", [{
            "region": args.region, "population": args.population, "n_tests": 0, "n_sig_fdr": 0,
            "error": "no_tests_completed",
        }])
        pd.DataFrame(failed).to_csv(outdir / "risk_variant_cpg_failures.tsv", sep="\t", index=False)
        raise SystemExit("No risk-variant–CpG tests completed")

    res["qval"] = bh_fdr(res["pval_nominal"].to_numpy())
    res["significant_fdr"] = res["qval"] <= args.fdr
    # Attach predictability class
    if res["local_predictability"].notna().any():
        q80 = res.drop_duplicates("vmr_id")["local_predictability"].quantile(0.8)
        q20 = res.drop_duplicates("vmr_id")["local_predictability"].quantile(0.2)
        res["predictability_class"] = np.where(
            res["local_predictability"] >= q80, "high",
            np.where(res["local_predictability"] <= q20, "low", "mid"),
        )
    else:
        res["predictability_class"] = "unknown"

    out = outdir / "risk_variant_cpg_meqtl.tsv.gz"
    res.sort_values("pval_nominal").to_csv(out, sep="\t", index=False, compression="gzip")
    pd.DataFrame(failed).to_csv(outdir / "risk_variant_cpg_failures.tsv", sep="\t", index=False)

    sig = res[res["significant_fdr"]]
    write_tsv(outdir / "risk_variant_cpg_tests_summary.tsv", [{
        "region": args.region,
        "population": args.population,
        "window_bp": args.window,
        "maf_threshold": args.maf,
        "fdr_threshold": args.fdr,
        "n_index_tested": int(res["index_snp"].nunique()),
        "n_tests": int(len(res)),
        "n_sig_fdr": int(sig["significant_fdr"].sum()),
        "n_sig_cpgs": int(sig["phenotype_id"].nunique()),
        "n_sig_vmrs": int(sig["vmr_id"].nunique()),
        "n_sig_loci": int(sig["locus_id"].nunique()) if sig["locus_id"].notna().any() else 0,
        "covariates": str(cov_path),
        "genotype_prefix": str(geno_prefix),
        "phenotype_bed": str(pheno_bed),
        "output": str(out),
    }])
    print(
        f"Wrote {out}; tests={len(res)}; FDR<{args.fdr}={int(sig['significant_fdr'].sum())}; "
        f"loci={int(sig['locus_id'].nunique()) if len(sig) else 0}"
    )


if __name__ == "__main__":
    main()
