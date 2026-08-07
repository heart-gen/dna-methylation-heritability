#!/usr/bin/env python3
"""Phase 7 Analysis 5: shared-donor genotype × region interactions.

For donors with caudate + DLPFC + hippocampus samples, residualize CpG
methylation within each region on locked M3a covariates, then fit:

  resid ~ dosage + C(region) + dosage:C(region)

Primary screen: union of FDR-significant risk–CpG pairs across regions.
Uses donor-clustered (HC1) SEs. Does not claim caudate selectivity from
significance patterns alone.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
PHASE1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping"
P7 = PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m"
REGIONS = ["caudate", "dlpfc", "hippocampus"]
# Contrasts relative to caudate (reference)
REGION_NONREF = ["dlpfc", "hippocampus"]
SEED = 20260801


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=str(P7 / "gxregion"))
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument(
        "--pair-set",
        choices=["union_sig", "prioritized_all_pairs", "both"],
        default="both",
        help="Which risk–CpG pairs to test",
    )
    p.add_argument("--min-regions", type=int, default=3, help="Require this many regions with data")
    return p.parse_args()


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


def shared_donors() -> list[str]:
    ids = {}
    for r in REGIONS:
        p = PHASE1 / r / "_m/preflight/sample_inclusion_primary.tsv"
        df = pd.read_csv(p, sep="\t")
        ids[r] = set(df["brnum"].astype(str))
    return sorted(ids["caudate"] & ids["dlpfc"] & ids["hippocampus"])


def load_covariates(path: Path) -> pd.DataFrame:
    cov = pd.read_csv(path, sep="\t", index_col=0)
    if not cov.index.astype(str).str.startswith("Br").any() and cov.columns.astype(str).str.startswith("Br").any():
        cov = cov.T
    cov.index = cov.index.astype(str)
    return cov.apply(pd.to_numeric, errors="coerce")


def residualize_matrix(Y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Residualize columns of Y (n×p) on covariates X (n×k)."""
    n = X.shape[0]
    Xd = np.column_stack([np.ones(n), X])
    beta, _, _, _ = np.linalg.lstsq(Xd, Y, rcond=None)
    return Y - Xd @ beta


def collect_pairs(pair_set: str) -> pd.DataFrame:
    frames = []
    for r in REGIONS:
        df = pd.read_csv(P7 / r / "risk_variant_cpg_meqtl.tsv.gz", sep="\t")
        df["phenotype_id"] = df["phenotype_id"].astype(str)
        df["variant_id"] = df["variant_id"].astype(str)
        df["index_snp"] = df["index_snp"].astype(str)
        df["sig"] = df["significant_fdr"].astype(bool)
        df["source_region"] = r
        frames.append(df)
    allp = pd.concat(frames, ignore_index=True)

    if pair_set in {"union_sig", "both"}:
        sig = allp[allp["sig"]][
            ["index_snp", "phenotype_id", "variant_id", "locus_id", "vmr_id"]
        ].drop_duplicates()
        sig["pair_set"] = "union_sig"
    else:
        sig = pd.DataFrame()

    if pair_set in {"prioritized_all_pairs", "both"}:
        prior = pd.read_csv(P7 / "prioritized/prioritized_loci.tsv", sep="\t")
        snps = set(prior["index_snp"].astype(str))
        pr = allp[allp["index_snp"].isin(snps)][
            ["index_snp", "phenotype_id", "variant_id", "locus_id", "vmr_id"]
        ].drop_duplicates()
        pr["pair_set"] = "prioritized_all_pairs"
    else:
        pr = pd.DataFrame()

    if pair_set == "union_sig":
        out = sig
    elif pair_set == "prioritized_all_pairs":
        out = pr
    else:
        out = pd.concat([sig, pr], ignore_index=True).drop_duplicates(
            subset=["index_snp", "phenotype_id", "variant_id"]
        )
        # mark if in either set
        keys = set(map(tuple, sig[["index_snp", "phenotype_id", "variant_id"]].to_numpy()))
        out["in_union_sig"] = [
            (a, b, c) in keys
            for a, b, c in out[["index_snp", "phenotype_id", "variant_id"]].to_numpy()
        ]
    return out.reset_index(drop=True)


def load_region_matrices(donors: list[str], phenotypes: list[str]):
    """Return per-region residual DataFrames (donors × phenotypes) and common donors."""
    from tensorqtl import pgen
    import tensorqtl

    region_resid = {}
    region_present = {}
    donor_set = set(donors)

    for r in REGIONS:
        bed = PHASE1 / r / "_m/prepared/cpg_phenotypes.all_autosomes.bed.gz"
        cov_path = PHASE1 / r / "_m/prepared/covariates.txt"
        geno = str(PHASE1 / r / "_m/genotypes/meqtl_AA")

        pheno_df, _ = tensorqtl.read_phenotype_bed(str(bed))
        pheno_df.columns = pheno_df.columns.astype(str)
        keep_ph = [p for p in phenotypes if p in pheno_df.index]
        if not keep_ph:
            region_resid[r] = pd.DataFrame(index=[])
            region_present[r] = set()
            print(f"{r}: no overlapping phenotypes")
            continue
        pheno_df = pheno_df.loc[keep_ph]

        cov = load_covariates(cov_path)
        probe = pgen.PgenReader(geno)
        geno_ids = set(map(str, probe.sample_ids))
        present = [d for d in donors if d in pheno_df.columns and d in cov.index and d in geno_ids]
        if not present:
            region_resid[r] = pd.DataFrame(index=[])
            region_present[r] = set()
            print(f"{r}: no shared donors after filters")
            continue

        pgr = pgen.PgenReader(geno, select_samples=present)
        samples = list(map(str, pgr.sample_ids))  # native genotype order
        assert set(samples) <= donor_set

        Y = pheno_df[samples].to_numpy(dtype=float).T  # n × p
        col_mean = np.nanmean(Y, axis=0)
        nan_idx = np.where(np.isnan(Y))
        if nan_idx[0].size:
            Y[nan_idx] = np.take(col_mean, nan_idx[1])
        X = cov.loc[samples]
        if X.isnull().any().any():
            X = X.dropna(axis=1, how="any")
        Y_res = residualize_matrix(Y, X.to_numpy(dtype=float))
        region_resid[r] = pd.DataFrame(Y_res, index=samples, columns=keep_ph)
        region_present[r] = set(samples)
        print(f"{r}: shared donors with complete pheno/cov/geno = {len(samples)}; phenotypes={len(keep_ph)}")

    common = set(donors)
    for r in REGIONS:
        common &= region_present[r]
    common = sorted(common)
    print(f"Donors with all {len(REGIONS)} regions after data filters: {len(common)}")
    return region_resid, common


def fit_gxregion(long: pd.DataFrame) -> dict:
    """Fit resid ~ g + C(region) + g:C(region) with caudate reference; cluster by donor."""
    d = long.dropna(subset=["resid", "g", "region", "donor"]).copy()
    if d["region"].nunique() < 2 or d["donor"].nunique() < 20:
        return {"error": "insufficient_data", "n_obs": len(d), "n_donors": int(d["donor"].nunique()) if len(d) else 0}
    d["region"] = pd.Categorical(d["region"], categories=REGIONS, ordered=True)
    # Require genotype variance
    if d["g"].std(ddof=0) < 1e-8:
        return {"error": "zero_genotype_variance", "n_obs": len(d), "n_donors": int(d["donor"].nunique())}

    try:
        model = smf.ols("resid ~ g + C(region) + g:C(region)", data=d)
        # cluster-robust by donor
        res = model.fit(cov_type="cluster", cov_kwds={"groups": d["donor"]})
    except Exception as exc:  # noqa: BLE001
        return {"error": f"fit_failed:{exc}", "n_obs": len(d), "n_donors": int(d["donor"].nunique())}

    params = res.params
    bse = res.bse
    pvals = res.pvalues

    # Main effect of g is caudate (reference)
    out = {
        "n_obs": int(res.nobs),
        "n_donors": int(d["donor"].nunique()),
        "n_regions": int(d["region"].nunique()),
        "beta_g_caudate": float(params.get("g", np.nan)),
        "se_g_caudate": float(bse.get("g", np.nan)),
        "pval_g_caudate": float(pvals.get("g", np.nan)),
        "converged": True,
    }
    # Interaction terms: g:C(region)[T.dlpfc], etc.
    inter_names = []
    for r in REGION_NONREF:
        key = f"g:C(region)[T.{r}]"
        # statsmodels may name as g:C(region)[T.dlpfc]
        if key not in params:
            # try alternate
            alts = [k for k in params.index if k.startswith("g:C(region)") and r in k]
            key = alts[0] if alts else key
        out[f"beta_g_x_{r}"] = float(params.get(key, np.nan))
        out[f"se_g_x_{r}"] = float(bse.get(key, np.nan))
        out[f"pval_g_x_{r}"] = float(pvals.get(key, np.nan))
        # region-specific total effect = g + interaction
        out[f"beta_g_{r}"] = out["beta_g_caudate"] + out[f"beta_g_x_{r}"]
        inter_names.append(key)

    # Joint Wald test for all interaction terms
    inter_names = [k for k in inter_names if k in params.index]
    if inter_names:
        try:
            w = res.wald_test(",".join([f"({k} = 0)" for k in inter_names]), scalar=True)
            out["wald_g_x_region_stat"] = float(np.asarray(w.statistic).squeeze())
            out["wald_g_x_region_df"] = float(np.asarray(w.df_num).squeeze()) if hasattr(w, "df_num") else float(len(inter_names))
            out["wald_g_x_region_pval"] = float(np.asarray(w.pvalue).squeeze())
        except Exception:
            # fallback manual
            try:
                R = np.zeros((len(inter_names), len(params)))
                for i, k in enumerate(inter_names):
                    R[i, list(params.index).index(k)] = 1.0
                w = res.wald_test(R, scalar=True)
                out["wald_g_x_region_stat"] = float(np.asarray(w.statistic).squeeze())
                out["wald_g_x_region_pval"] = float(np.asarray(w.pvalue).squeeze())
                out["wald_g_x_region_df"] = float(len(inter_names))
            except Exception as exc:  # noqa: BLE001
                out["wald_g_x_region_pval"] = np.nan
                out["wald_error"] = str(exc)
    else:
        out["wald_g_x_region_pval"] = np.nan

    # Pairwise simple effects already in beta_g_*
    return out


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    donors = shared_donors()
    write_tsv(outdir / "shared_donors_input.tsv", [{"brnum": d} for d in donors])
    print(f"Shared donors (inclusion lists): {len(donors)}")

    pairs = collect_pairs(args.pair_set)
    pairs.to_csv(outdir / "pairs_tested_definition.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"Candidate pairs: {len(pairs)}")

    phenotypes = sorted(pairs["phenotype_id"].astype(str).unique())
    variants = sorted(pairs["variant_id"].astype(str).unique())

    region_resid, common = load_region_matrices(donors, phenotypes)
    if len(common) < 50:
        raise SystemExit(f"Too few complete shared donors after filters: {len(common)}")
    write_tsv(outdir / "shared_donors_complete.tsv", [{"brnum": d} for d in common])

    # Load dosages once from caudate panel (same genotypes across regions for a donor)
    from tensorqtl import pgen
    geno = str(PHASE1 / "caudate/_m/genotypes/meqtl_AA")
    pgr = pgen.PgenReader(geno, select_samples=common)
    samples = list(map(str, pgr.sample_ids))
    # Align common to genotype order
    common = samples
    panel = set(map(str, pgr.variant_ids))
    var_ok = [v for v in variants if v in panel]
    print(f"Variants in panel: {len(var_ok)}/{len(variants)}")
    G = {
        v: np.asarray(pgr.read_dosages(v), dtype=float).ravel()
        for v in var_ok
    }
    for v in list(G):
        g = G[v]
        if np.isnan(g).any():
            g = np.where(np.isnan(g), np.nanmean(g), g)
            G[v] = g

    # Restrict residual matrices to common donor order
    for r in REGIONS:
        region_resid[r] = region_resid[r].reindex(index=common)

    results = []
    n_skip = 0
    for i, row in pairs.iterrows():
        pid = str(row["phenotype_id"])
        vid = str(row["variant_id"])
        if vid not in G:
            n_skip += 1
            continue
        # Build long table
        chunks = []
        for r in REGIONS:
            if pid not in region_resid[r].columns:
                continue
            y = region_resid[r][pid]
            if y.isna().all():
                continue
            chunks.append(pd.DataFrame({
                "donor": common,
                "region": r,
                "resid": y.to_numpy(dtype=float),
                "g": G[vid],
            }))
        if not chunks:
            n_skip += 1
            continue
        long = pd.concat(chunks, ignore_index=True)
        n_reg = long.groupby("region").ngroups
        if n_reg < args.min_regions:
            n_skip += 1
            continue
        fit = fit_gxregion(long)
        fit.update({
            "index_snp": row["index_snp"],
            "phenotype_id": pid,
            "variant_id": vid,
            "locus_id": row.get("locus_id", pd.NA),
            "vmr_id": row.get("vmr_id", pd.NA),
            "pair_set": row.get("pair_set", args.pair_set),
            "in_union_sig": row.get("in_union_sig", True),
        })
        results.append(fit)
        if len(results) % 50 == 0:
            print(f"  fitted {len(results)} pairs...")

    res = pd.DataFrame(results)
    if res.empty:
        raise SystemExit("No G×region models fitted")

    # FDR on joint interaction p-values among successfully fitted models
    ok = res["wald_g_x_region_pval"].notna() if "wald_g_x_region_pval" in res.columns else pd.Series(False, index=res.index)
    res["qval_g_x_region"] = np.nan
    if ok.any():
        res.loc[ok, "qval_g_x_region"] = bh_fdr(res.loc[ok, "wald_g_x_region_pval"].to_numpy())
    res["sig_g_x_region_fdr"] = res["qval_g_x_region"] <= args.fdr

    # Also FDR on pairwise interaction terms
    for r in REGION_NONREF:
        col = f"pval_g_x_{r}"
        qcol = f"qval_g_x_{r}"
        if col in res.columns:
            mask = res[col].notna()
            res[qcol] = np.nan
            res.loc[mask, qcol] = bh_fdr(res.loc[mask, col].to_numpy())
            res[f"sig_g_x_{r}_fdr"] = res[qcol] <= args.fdr

    out_path = outdir / "gxregion_pair_results.tsv.gz"
    res.sort_values("wald_g_x_region_pval", na_position="last").to_csv(
        out_path, sep="\t", index=False, compression="gzip"
    )

    n_sig = int(res["sig_g_x_region_fdr"].fillna(False).sum())
    write_tsv(outdir / "gxregion_summary.tsv", [{
        "n_shared_donors_input": len(donors),
        "n_shared_donors_complete": len(common),
        "n_candidate_pairs": int(len(pairs)),
        "n_pairs_fitted": int(len(res)),
        "n_pairs_skipped": int(n_skip),
        "n_sig_joint_gxregion_fdr": n_sig,
        "n_sig_g_x_dlpfc_fdr": int(res.get("sig_g_x_dlpfc_fdr", pd.Series(dtype=bool)).fillna(False).sum()) if "sig_g_x_dlpfc_fdr" in res.columns else 0,
        "n_sig_g_x_hippocampus_fdr": int(res.get("sig_g_x_hippocampus_fdr", pd.Series(dtype=bool)).fillna(False).sum()) if "sig_g_x_hippocampus_fdr" in res.columns else 0,
        "fdr": args.fdr,
        "model": "resid ~ g + C(region) + g:C(region); cluster(donor); within-region M3a residualization",
        "reference_region": "caudate",
        "output": str(out_path),
    }])

    # Prioritized locus rollup
    prior_path = P7 / "prioritized/prioritized_loci.tsv"
    if prior_path.exists():
        prior = pd.read_csv(prior_path, sep="\t")
        snps = set(prior["index_snp"].astype(str))
        sub = res[res["index_snp"].isin(snps)].copy()
        if len(sub):
            roll = sub.groupby("index_snp", as_index=False).agg(
                n_pairs=("phenotype_id", "size"),
                n_sig_joint=("sig_g_x_region_fdr", "sum"),
                min_wald_p=("wald_g_x_region_pval", "min"),
                min_wald_q=("qval_g_x_region", "min"),
                median_beta_caudate=("beta_g_caudate", "median"),
                median_beta_dlpfc=("beta_g_dlpfc", "median"),
                median_beta_hippocampus=("beta_g_hippocampus", "median"),
            )
            roll = roll.merge(prior[["index_snp", "locus_id", "rank", "priority_score"]], on="index_snp", how="left")
            roll.to_csv(outdir / "gxregion_prioritized_locus_summary.tsv", sep="\t", index=False)

    # Claim snapshot
    # Among union_sig pairs with 3-region data, what fraction show joint interaction
    write_tsv(outdir / "gxregion_claim_snapshot.tsv", [{
        "n_pairs_fitted": int(len(res)),
        "n_sig_joint_interaction_fdr": n_sig,
        "frac_sig_joint_interaction": float(n_sig / len(res)) if len(res) else np.nan,
        "n_sig_interaction_dlpfc_fdr": int(res["sig_g_x_dlpfc_fdr"].fillna(False).sum()) if "sig_g_x_dlpfc_fdr" in res.columns else 0,
        "n_sig_interaction_hippocampus_fdr": int(res["sig_g_x_hippocampus_fdr"].fillna(False).sum()) if "sig_g_x_hippocampus_fdr" in res.columns else 0,
        "interpretation": (
            "Significant genotype×region interactions support regional effect heterogeneity "
            "among shared donors; absence of interaction means region differences in discovery "
            "are not detectable as effect heterogeneity in the shared-donor subset."
            if n_sig > 0 else
            "No FDR-significant genotype×region interactions in the shared-donor subset; "
            "regional discovery differences should not be described as caudate-selective "
            "effect heterogeneity without further evidence."
        ),
    }])
    print(f"Wrote {out_path}; fitted={len(res)}; joint FDR sig={n_sig}; skipped={n_skip}")


if __name__ == "__main__":
    main()
