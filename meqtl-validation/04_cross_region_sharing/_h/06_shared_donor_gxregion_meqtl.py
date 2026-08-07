#!/usr/bin/env python3
"""Phase 4 Experiment 3: shared-donor genotype × region interactions (architecture).

For donors with caudate + DLPFC + hippocampus samples, residualize CpG
methylation within each region on locked M3a covariates, then fit:

  resid ~ dosage + C(region) + dosage:C(region)

Primary screen: stratified sample of FDR-significant lead SNP–CpG pairs from
Phase 1 (architecture-scale; not SCZ-targeted). Uses donor-clustered (HC1) SEs.
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
PHASE2 = PROJECT / "meqtl-validation/02_vmr_meqtl_burden/_m"
CROSS = PROJECT / "meqtl-validation/04_cross_region_sharing/_m"
OUTDIR = CROSS / "gxregion"
REGIONS = ["caudate", "dlpfc", "hippocampus"]
REGION_NONREF = ["dlpfc", "hippocampus"]
SEED = 20260805
FDR = 0.05
MAX_PAIRS = 3000


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=str(OUTDIR))
    p.add_argument("--fdr", type=float, default=FDR)
    p.add_argument("--max-pairs", type=int, default=MAX_PAIRS)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--min-regions", type=int, default=3)
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
    if (
        not cov.index.astype(str).str.startswith("Br").any()
        and cov.columns.astype(str).str.startswith("Br").any()
    ):
        cov = cov.T
    cov.index = cov.index.astype(str)
    return cov.apply(pd.to_numeric, errors="coerce")


def residualize_matrix(Y: np.ndarray, X: np.ndarray) -> np.ndarray:
    n = X.shape[0]
    Xd = np.column_stack([np.ones(n), X])
    beta, _, _, _ = np.linalg.lstsq(Xd, Y, rcond=None)
    return Y - Xd @ beta


def load_leads(region: str, fdr: float) -> pd.DataFrame:
    path = PHASE1 / region / "_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz"
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["phenotype_id", "variant_id", "slope", "slope_se", "qval", "pval_nominal"],
    )
    df["phenotype_id"] = df["phenotype_id"].astype(str)
    df["variant_id"] = df["variant_id"].astype(str)
    df["qval"] = pd.to_numeric(df["qval"], errors="coerce")
    df["sig"] = df["qval"].le(fdr).fillna(False)
    df["region"] = region
    return df


def load_cpg_vmr_predictability() -> pd.DataFrame:
    """Map phenotype_id → vmr_id → local_predictability (caudate burden)."""
    maps = sorted((PHASE1 / "caudate/_m/prepared").glob("cpg_vmr_map.chr*.tsv"))
    if not maps:
        return pd.DataFrame(columns=["phenotype_id", "vmr_id", "local_predictability"])
    cpg = pd.concat([pd.read_csv(p, sep="\t") for p in maps], ignore_index=True)
    cpg["phenotype_id"] = cpg["phenotype_id"].astype(str)
    cpg["vmr_id"] = cpg["vmr_id"].astype(str)
    bur = pd.read_csv(PHASE2 / "caudate/vmr_meqtl_burden.tsv.gz", sep="\t")
    bur["vmr_id"] = bur["vmr_id"].astype(str)
    bur["local_predictability"] = pd.to_numeric(bur["local_predictability"], errors="coerce")
    out = cpg.merge(bur[["vmr_id", "local_predictability"]], on="vmr_id", how="left")
    return out[["phenotype_id", "vmr_id", "local_predictability"]].drop_duplicates("phenotype_id")


def select_pairs(fdr: float, max_pairs: int, seed: int) -> pd.DataFrame:
    """Stratified sample of FDR-sig leads; prefer caudate lead SNP when available."""
    leads = {r: load_leads(r, fdr) for r in REGIONS}
    # Phenotypes tested in all three regions
    common_ph = (
        set(leads["caudate"]["phenotype_id"])
        & set(leads["dlpfc"]["phenotype_id"])
        & set(leads["hippocampus"]["phenotype_id"])
    )
    # Union of significant phenotypes
    sig_ph = set()
    for r, df in leads.items():
        sig_ph |= set(df.loc[df["sig"], "phenotype_id"])
    candidates = sorted(common_ph & sig_ph)
    print(f"Shared-tested phenotypes with ≥1-region FDR sig: {len(candidates)}")

    caud = leads["caudate"].drop_duplicates("phenotype_id").set_index("phenotype_id")
    dlpfc_idx = leads["dlpfc"].drop_duplicates("phenotype_id").set_index("phenotype_id")
    hip_idx = leads["hippocampus"].drop_duplicates("phenotype_id").set_index("phenotype_id")
    pred = load_cpg_vmr_predictability().set_index("phenotype_id")

    rows = []
    for pid in candidates:
        crow = caud.loc[pid]
        n_sig_regions = int(bool(crow["sig"]))
        if pid in dlpfc_idx.index:
            n_sig_regions += int(bool(dlpfc_idx.loc[pid, "sig"]))
        if pid in hip_idx.index:
            n_sig_regions += int(bool(hip_idx.loc[pid, "sig"]))
        lp = float(pred.loc[pid, "local_predictability"]) if pid in pred.index else np.nan
        vmr = str(pred.loc[pid, "vmr_id"]) if pid in pred.index else ""
        rows.append({
            "phenotype_id": pid,
            "variant_id": str(crow["variant_id"]),
            "vmr_id": vmr,
            "local_predictability": lp,
            "caudate_sig": bool(crow["sig"]),
            "n_sig_regions": n_sig_regions,
            "caudate_qval": float(crow["qval"]),
        })
    pool = pd.DataFrame(rows)
    if pool.empty:
        raise SystemExit("No candidate pairs for G×region")

    # Terciles of predictability among finite values
    finite = pool["local_predictability"].notna()
    pool["pred_tercile"] = "unknown"
    if finite.sum() >= 30:
        try:
            pool.loc[finite, "pred_tercile"] = pd.qcut(
                pool.loc[finite, "local_predictability"],
                q=3,
                labels=["low", "mid", "high"],
                duplicates="drop",
            ).astype(str)
        except ValueError:
            pool.loc[finite, "pred_tercile"] = "unbinned"

    rng = np.random.default_rng(seed)
    per_bin = max(max_pairs // max(pool["pred_tercile"].nunique(), 1), 1)
    sampled = []
    for terc, grp in pool.groupby("pred_tercile"):
        # Prefer multi-region significant and caudate-significant
        g = grp.sort_values(["n_sig_regions", "caudate_sig", "caudate_qval"], ascending=[False, False, True])
        take_n = min(per_bin, len(g))
        # Keep top half by multi-region support, random fill from remainder if needed
        top = g.head(max(take_n // 2, 1))
        rest = g.iloc[len(top):]
        n_fill = take_n - len(top)
        if n_fill > 0 and len(rest):
            idx = rng.choice(rest.index.to_numpy(), size=min(n_fill, len(rest)), replace=False)
            fill = rest.loc[idx]
            chosen = pd.concat([top, fill], ignore_index=True)
        else:
            chosen = top.head(take_n)
        chosen = chosen.copy()
        chosen["pred_tercile"] = terc
        sampled.append(chosen)
    out = pd.concat(sampled, ignore_index=True).drop_duplicates("phenotype_id")
    if len(out) > max_pairs:
        out = out.sample(n=max_pairs, random_state=seed)
    out["pair_set"] = "architecture_stratified_sig"
    print(
        f"Selected {len(out)} pairs "
        f"(tercile counts: {out['pred_tercile'].value_counts().to_dict()})"
    )
    return out.reset_index(drop=True)


def load_region_matrices(donors: list[str], phenotypes: list[str]):
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
        keep_ph = [p for p in phenotypes if p in pheno_df.index.astype(str)]
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
            continue

        pgr = pgen.PgenReader(geno, select_samples=present)
        samples = list(map(str, pgr.sample_ids))
        assert set(samples) <= donor_set

        Y = pheno_df[samples].to_numpy(dtype=float).T
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
        print(f"{r}: shared donors with data = {len(samples)}; phenotypes={len(keep_ph)}")

    common = set(donors)
    for r in REGIONS:
        common &= region_present[r]
    common = sorted(common)
    print(f"Donors with all {len(REGIONS)} regions after filters: {len(common)}")
    return region_resid, common


def fit_gxregion(long: pd.DataFrame) -> dict:
    d = long.dropna(subset=["resid", "g", "region", "donor"]).copy()
    if d["region"].nunique() < 2 or d["donor"].nunique() < 20:
        return {"error": "insufficient_data", "n_obs": len(d), "n_donors": int(d["donor"].nunique()) if len(d) else 0}
    d["region"] = pd.Categorical(d["region"], categories=REGIONS, ordered=True)
    if d["g"].std(ddof=0) < 1e-8:
        return {"error": "zero_genotype_variance", "n_obs": len(d), "n_donors": int(d["donor"].nunique())}

    try:
        model = smf.ols("resid ~ g + C(region) + g:C(region)", data=d)
        res = model.fit(cov_type="cluster", cov_kwds={"groups": d["donor"]})
    except Exception as exc:  # noqa: BLE001
        return {"error": f"fit_failed:{exc}", "n_obs": len(d), "n_donors": int(d["donor"].nunique())}

    params = res.params
    bse = res.bse
    pvals = res.pvalues
    out = {
        "n_obs": int(res.nobs),
        "n_donors": int(d["donor"].nunique()),
        "n_regions": int(d["region"].nunique()),
        "beta_g_caudate": float(params.get("g", np.nan)),
        "se_g_caudate": float(bse.get("g", np.nan)),
        "pval_g_caudate": float(pvals.get("g", np.nan)),
        "converged": True,
    }
    inter_names = []
    for r in REGION_NONREF:
        key = f"g:C(region)[T.{r}]"
        if key not in params:
            alts = [k for k in params.index if k.startswith("g:C(region)") and r in k]
            key = alts[0] if alts else key
        out[f"beta_g_x_{r}"] = float(params.get(key, np.nan))
        out[f"se_g_x_{r}"] = float(bse.get(key, np.nan))
        out[f"pval_g_x_{r}"] = float(pvals.get(key, np.nan))
        out[f"beta_g_{r}"] = out["beta_g_caudate"] + out[f"beta_g_x_{r}"]
        inter_names.append(key)

    inter_names = [k for k in inter_names if k in params.index]
    if inter_names:
        try:
            w = res.wald_test(",".join([f"({k} = 0)" for k in inter_names]), scalar=True)
            out["wald_g_x_region_stat"] = float(np.asarray(w.statistic).squeeze())
            out["wald_g_x_region_df"] = float(len(inter_names))
            out["wald_g_x_region_pval"] = float(np.asarray(w.pvalue).squeeze())
        except Exception:
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
    return out


def update_pending(status: str, reason: str) -> None:
    path = CROSS / "pending_analyses.tsv"
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t")
    mask = df["analysis"] == "shared_donor_genotype_by_region"
    if mask.any():
        df.loc[mask, "status"] = status
        df.loc[mask, "reason"] = reason
        df.to_csv(path, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    donors = shared_donors()
    write_tsv(outdir / "shared_donors_input.tsv", [{"brnum": d} for d in donors])
    print(f"Shared donors (inclusion lists): {len(donors)}")

    pairs = select_pairs(args.fdr, args.max_pairs, args.seed)
    pairs.to_csv(outdir / "pairs_tested_definition.tsv.gz", sep="\t", index=False, compression="gzip")

    phenotypes = sorted(pairs["phenotype_id"].astype(str).unique())
    variants = sorted(pairs["variant_id"].astype(str).unique())

    region_resid, common = load_region_matrices(donors, phenotypes)
    if len(common) < 50:
        raise SystemExit(f"Too few complete shared donors after filters: {len(common)}")
    write_tsv(outdir / "shared_donors_complete.tsv", [{"brnum": d} for d in common])

    from tensorqtl import pgen
    geno = str(PHASE1 / "caudate/_m/genotypes/meqtl_AA")
    pgr = pgen.PgenReader(geno, select_samples=common)
    samples = list(map(str, pgr.sample_ids))
    common = samples
    panel = set(map(str, pgr.variant_ids))
    var_ok = [v for v in variants if v in panel]
    print(f"Variants in panel: {len(var_ok)}/{len(variants)}")
    G = {v: np.asarray(pgr.read_dosages(v), dtype=float).ravel() for v in var_ok}
    for v in list(G):
        g = G[v]
        if np.isnan(g).any():
            G[v] = np.where(np.isnan(g), np.nanmean(g), g)

    for r in REGIONS:
        region_resid[r] = region_resid[r].reindex(index=common)

    results = []
    n_skip = 0
    for _, row in pairs.iterrows():
        pid = str(row["phenotype_id"])
        vid = str(row["variant_id"])
        if vid not in G:
            n_skip += 1
            continue
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
        if long.groupby("region").ngroups < args.min_regions:
            n_skip += 1
            continue
        fit = fit_gxregion(long)
        fit.update({
            "phenotype_id": pid,
            "variant_id": vid,
            "vmr_id": row.get("vmr_id", pd.NA),
            "local_predictability": row.get("local_predictability", np.nan),
            "pred_tercile": row.get("pred_tercile", ""),
            "n_sig_regions": row.get("n_sig_regions", np.nan),
            "caudate_sig": row.get("caudate_sig", False),
            "pair_set": row.get("pair_set", "architecture_stratified_sig"),
        })
        results.append(fit)
        if len(results) % 100 == 0:
            print(f"  fitted {len(results)} pairs...")

    res = pd.DataFrame(results)
    if res.empty:
        raise SystemExit("No G×region models fitted")

    ok = res["wald_g_x_region_pval"].notna() if "wald_g_x_region_pval" in res.columns else pd.Series(False, index=res.index)
    res["qval_g_x_region"] = np.nan
    if ok.any():
        res.loc[ok, "qval_g_x_region"] = bh_fdr(res.loc[ok, "wald_g_x_region_pval"].to_numpy())
    res["sig_g_x_region_fdr"] = res["qval_g_x_region"] <= args.fdr

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
        "frac_sig_joint": float(n_sig / len(res)) if len(res) else np.nan,
        "n_sig_g_x_dlpfc_fdr": int(res["sig_g_x_dlpfc_fdr"].fillna(False).sum()) if "sig_g_x_dlpfc_fdr" in res.columns else 0,
        "n_sig_g_x_hippocampus_fdr": int(res["sig_g_x_hippocampus_fdr"].fillna(False).sum()) if "sig_g_x_hippocampus_fdr" in res.columns else 0,
        "fdr": args.fdr,
        "max_pairs": args.max_pairs,
        "seed": args.seed,
        "model": "resid ~ g + C(region) + g:C(region); cluster(donor); within-region M3a residualization",
        "reference_region": "caudate",
        "output": str(out_path),
    }])

    # Stratify by predictability tercile
    if "pred_tercile" in res.columns:
        roll = (
            res.groupby("pred_tercile", as_index=False)
            .agg(
                n_pairs=("phenotype_id", "size"),
                n_sig_joint=("sig_g_x_region_fdr", "sum"),
                median_wald_p=("wald_g_x_region_pval", "median"),
                median_beta_caudate=("beta_g_caudate", "median"),
            )
        )
        roll["frac_sig_joint"] = roll["n_sig_joint"] / roll["n_pairs"]
        roll.to_csv(outdir / "gxregion_by_predictability_tercile.tsv", sep="\t", index=False)

    write_tsv(outdir / "gxregion_claim_snapshot.tsv", [{
        "n_pairs_fitted": int(len(res)),
        "n_sig_joint_interaction_fdr": n_sig,
        "frac_sig_joint_interaction": float(n_sig / len(res)) if len(res) else np.nan,
        "n_sig_interaction_dlpfc_fdr": int(res["sig_g_x_dlpfc_fdr"].fillna(False).sum()) if "sig_g_x_dlpfc_fdr" in res.columns else 0,
        "n_sig_interaction_hippocampus_fdr": int(res["sig_g_x_hippocampus_fdr"].fillna(False).sum()) if "sig_g_x_hippocampus_fdr" in res.columns else 0,
        "interpretation": (
            "Significant genotype×region interactions support regional effect heterogeneity "
            "among shared donors for a minority of architecture-screened pairs; absence of "
            "interaction for most pairs means regional discovery differences should not be "
            "described as widespread caudate-selective effect heterogeneity."
            if n_sig > 0 else
            "No FDR-significant genotype×region interactions in the shared-donor architecture "
            "screen; regional discovery differences should not be described as caudate-selective "
            "effect heterogeneity without further evidence."
        ),
    }])

    update_pending(
        status="done",
        reason=(
            f"Shared-donor G×region complete (fitted={len(res)}, joint FDR sig={n_sig}, "
            f"shared_donors={len(common)})"
        ),
    )
    print(f"Wrote {out_path}; fitted={len(res)}; joint FDR sig={n_sig}; skipped={n_skip}")


if __name__ == "__main__":
    main()
