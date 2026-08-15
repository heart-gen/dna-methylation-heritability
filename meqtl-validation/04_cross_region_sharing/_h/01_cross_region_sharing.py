#!/usr/bin/env python3
"""Phase 4 cross-region meQTL architecture.

Uses Phase 1 lead cis-meQTL tables and Phase 2 VMR burden tables (AA discovery).
Does not re-run tensorQTL. Formal caudate donor-downsampling remapping is
reported as pending; shared-CpG and N-normalized discovery rates are provided
as power-aware sensitivities.
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import ANALYSIS_SCHEMA_VERSION, canonical_vmr_id, parse_vmr_coordinate, write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
PHASE1 = PROJECT / "meqtl-validation" / "01_cpg_meqtl_mapping"
PHASE2 = PROJECT / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m"
OUTDIR = PROJECT / "meqtl-validation" / "04_cross_region_sharing" / "_m"
REGIONS = ["caudate", "dlpfc", "hippocampus"]
FDR = 0.05
SEED = 20260730
SAMPLE_N = {"caudate": 153, "dlpfc": 111, "hippocampus": 116}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fdr", type=float, default=FDR)
    p.add_argument("--seed", type=int, default=SEED)
    return p.parse_args()


def _z(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x * 0.0
    return (x - x.mean()) / sd


def load_burden(region: str) -> pd.DataFrame:
    path = PHASE2 / region / "vmr_meqtl_burden.tsv.gz"
    df = pd.read_csv(path, sep="\t")
    if "analysis_schema_version" not in df or not df["analysis_schema_version"].eq(ANALYSIS_SCHEMA_VERSION).all():
        raise SystemExit(f"{region}: stale Phase 2 burden; regenerate repair schema v2 before Phase 4")
    df["vmr_id"] = df["vmr_id"].astype(str)
    df["meqtl_supported"] = pd.to_numeric(df["n_cpgs_with_sig_meqtl"], errors="coerce").fillna(0).gt(0)
    df["local_predictability"] = pd.to_numeric(df["local_predictability"], errors="coerce")
    df["proportion_cpgs_with_sig_meqtl"] = pd.to_numeric(
        df["proportion_cpgs_with_sig_meqtl"], errors="coerce"
    )
    # Region-local task IDs are not cross-region locus identifiers. Resolve every
    # row to its hg38 interval before any sharing analysis.
    pred_path = (
        PROJECT / "heritability/elastic_net_model/all_individuals" / region / "_m"
        / f"{region}_summary_elastic-net_AA.tsv"
    )
    pred = pd.read_csv(pred_path, sep="\t")
    pred["task_id"] = pred["task_id"].astype(str)
    pred["vmr_coord_id"] = [
        canonical_vmr_id(c, s, e) for c, s, e in zip(pred["chrom"], pred["start"], pred["end"])
    ]
    task_to_coord = pred.drop_duplicates("task_id").set_index("task_id")["vmr_coord_id"]
    df["vmr_coord_id"] = df["vmr_id"].map(task_to_coord)
    unresolved = df["vmr_coord_id"].isna()
    if unresolved.any():
        df.loc[unresolved, "vmr_coord_id"] = df.loc[unresolved, "vmr_id"].map(
            lambda value: (
                canonical_vmr_id(*parsed) if (parsed := parse_vmr_coordinate(value)) else np.nan
            )
        )
    if df["vmr_coord_id"].isna().any():
        raise SystemExit(
            f"{region}: {int(df['vmr_coord_id'].isna().sum())} burden rows lack canonical coordinates"
        )
    return df


def load_leads(region: str, fdr: float) -> pd.DataFrame:
    path = PHASE1 / region / "_m" / "tensorqtl" / "qc" / "lead_snp_per_cpg.tsv.gz"
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["phenotype_id", "variant_id", "slope", "slope_se", "qval", "pval_nominal"],
    )
    df["phenotype_id"] = df["phenotype_id"].astype(str)
    df["variant_id"] = df["variant_id"].astype(str)
    df["slope"] = pd.to_numeric(df["slope"], errors="coerce")
    df["slope_se"] = pd.to_numeric(df["slope_se"], errors="coerce")
    df["qval"] = pd.to_numeric(df["qval"], errors="coerce")
    df["sig"] = df["qval"].le(fdr)
    df["z"] = df["slope"] / df["slope_se"].replace(0, np.nan)
    return df


def load_shared_donors() -> dict:
    ids = {}
    for region in REGIONS:
        p = PHASE1 / region / "_m" / "preflight" / "sample_inclusion_primary.tsv"
        df = pd.read_csv(p, sep="\t")
        ids[region] = set(df["brnum"].astype(str))
    out = {"n_by_region": {r: len(ids[r]) for r in REGIONS}}
    for a, b in combinations(REGIONS, 2):
        out[f"n_shared_{a}_{b}"] = len(ids[a] & ids[b])
    out["n_shared_all3"] = len(ids["caudate"] & ids["dlpfc"] & ids["hippocampus"])
    return out


def pairwise_vmr_sharing(burdens: dict[str, pd.DataFrame]) -> tuple[list[dict], pd.DataFrame]:
    rows = []
    # wide support table for triple / predictability models
    base = None
    for region, df in burdens.items():
        sub = df[["vmr_coord_id", "vmr_id", "meqtl_supported", "local_predictability",
                  "proportion_cpgs_with_sig_meqtl", "n_tested_cpgs"]].copy()
        sub = sub.rename(columns={
            "vmr_id": f"vmr_id_{region}",
            "meqtl_supported": f"supported_{region}",
            "local_predictability": f"pred_{region}",
            "proportion_cpgs_with_sig_meqtl": f"prop_{region}",
            "n_tested_cpgs": f"n_cpg_{region}",
        })
        base = sub if base is None else base.merge(sub, on="vmr_coord_id", how="outer")

    for a, b in combinations(REGIONS, 2):
        sa, sb = f"supported_{a}", f"supported_{b}"
        use = base.dropna(subset=[sa, sb]).copy()
        sa_b = use[sa].astype(bool)
        sb_b = use[sb].astype(bool)
        both = sa_b & sb_b
        either = sa_b | sb_b
        rows.append({
            "contrast": f"{a}_vs_{b}",
            "level": "vmr",
            "n_shared_tested": int(len(use)),
            "n_supported_a": int(sa_b.sum()),
            "n_supported_b": int(sb_b.sum()),
            "n_supported_both": int(both.sum()),
            "jaccard_support": float(both.sum() / either.sum()) if either.sum() else np.nan,
            "frac_a_replicated_in_b": float((both.sum() / sa_b.sum()) if sa_b.sum() else np.nan),
            "frac_b_replicated_in_a": float((both.sum() / sb_b.sum()) if sb_b.sum() else np.nan),
        })
        # predictability (mean of available) vs shared support
        pred = use[[f"pred_{a}", f"pred_{b}"]].mean(axis=1, skipna=True)
        y = both.astype(int)
        ok = pred.notna()
        if ok.sum() > 50 and y[ok].nunique() > 1:
            X = sm.add_constant(_z(pred[ok]).to_frame("predictability"))
            try:
                res = sm.GLM(y[ok], X, family=sm.families.Binomial()).fit()
                rows[-1]["logistic_or_pred_shared"] = float(np.exp(res.params["predictability"]))
                rows[-1]["logistic_p_pred_shared"] = float(res.pvalues["predictability"])
                rows[-1]["spearman_pred_shared"] = float(spearmanr(pred[ok], y[ok]).statistic)
            except Exception as exc:  # noqa: BLE001
                rows[-1]["logistic_error"] = str(exc)

    # triple sharing
    cols = [f"supported_{r}" for r in REGIONS]
    trip = base.dropna(subset=cols).copy()
    all_sup = trip[cols[0]].astype(bool) & trip[cols[1]].astype(bool) & trip[cols[2]].astype(bool)
    any_sup = trip[cols[0]].astype(bool) | trip[cols[1]].astype(bool) | trip[cols[2]].astype(bool)
    rows.append({
        "contrast": "all3",
        "level": "vmr",
        "n_shared_tested": int(len(trip)),
        "n_supported_both": int(all_sup.sum()),
        "jaccard_support": float(all_sup.sum() / any_sup.sum()) if any_sup.sum() else np.nan,
        "n_supported_any": int(any_sup.sum()),
    })
    return rows, base


def pairwise_cpg_concordance(leads: dict[str, pd.DataFrame], fdr: float) -> list[dict]:
    rows = []
    for a, b in combinations(REGIONS, 2):
        da = leads[a].set_index("phenotype_id")
        db = leads[b].set_index("phenotype_id")
        shared = da.index.intersection(db.index)
        ma = da.loc[shared]
        mb = db.loc[shared]
        both_sig = ma["sig"] & mb["sig"]
        either_sig = ma["sig"] | mb["sig"]
        same_lead = ma["variant_id"] == mb["variant_id"]
        # Effect sizes and signs are commensurate only for the identical tested
        # variant with a common REF/ALT encoding. Independent lead variants are
        # retained for discovery-overlap counts but never used for concordance.
        both_same = both_sig & same_lead
        dir_same = np.sign(ma.loc[both_same, "slope"]) == np.sign(mb.loc[both_same, "slope"])

        def _corr(x, y):
            m = x.notna() & y.notna()
            if m.sum() < 20:
                return np.nan, np.nan, int(m.sum())
            if x[m].nunique() < 2 or y[m].nunique() < 2:
                return np.nan, np.nan, int(m.sum())
            r, p = pearsonr(x[m], y[m])
            return float(r), float(p), int(m.sum())

        r_slope, p_slope, n_slope = _corr(
            ma.loc[same_lead, "slope"], mb.loc[same_lead, "slope"]
        )
        r_z, p_z, n_z = _corr(ma.loc[same_lead, "z"], mb.loc[same_lead, "z"])
        r_slope_sig, p_slope_sig, n_slope_sig = _corr(
            ma.loc[either_sig & same_lead, "slope"], mb.loc[either_sig & same_lead, "slope"]
        )
        r_z_sig, p_z_sig, n_z_sig = _corr(
            ma.loc[either_sig & same_lead, "z"], mb.loc[either_sig & same_lead, "z"]
        )

        rows.append({
            "contrast": f"{a}_vs_{b}",
            "level": "cpg",
            "n_shared_tested": int(len(shared)),
            "n_sig_a": int(ma["sig"].sum()),
            "n_sig_b": int(mb["sig"].sum()),
            "n_sig_both": int(both_sig.sum()),
            "jaccard_sig": float(both_sig.sum() / either_sig.sum()) if either_sig.sum() else np.nan,
            "frac_a_replicated_in_b": float(both_sig.sum() / ma["sig"].sum()) if ma["sig"].sum() else np.nan,
            "frac_b_replicated_in_a": float(both_sig.sum() / mb["sig"].sum()) if mb["sig"].sum() else np.nan,
            "n_same_lead_variant": int(same_lead.sum()),
            "n_both_sig_same_lead": int(both_same.sum()),
            "n_both_sig_different_lead_not_compared": int((both_sig & ~same_lead).sum()),
            "direction_concordance_both_sig": float(dir_same.mean()) if both_same.sum() else np.nan,
            "direction_concordance_both_sig_same_lead": float(dir_same.mean()) if both_same.sum() else np.nan,
            "pearson_slope_all": r_slope,
            "pearson_slope_all_p": p_slope,
            "n_pearson_slope_all": n_slope,
            "pearson_z_all": r_z,
            "pearson_z_all_p": p_z,
            "n_pearson_z_all": n_z,
            "pearson_slope_either_sig": r_slope_sig,
            "pearson_slope_either_sig_p": p_slope_sig,
            "n_pearson_slope_either_sig": n_slope_sig,
            "pearson_z_either_sig": r_z_sig,
            "pearson_z_either_sig_p": p_z_sig,
            "n_pearson_z_either_sig": n_z_sig,
            "fdr": fdr,
            "effect_comparison_scope": "identical_lead_variant_only",
        })
    return rows


def discovery_power_summary(leads: dict[str, pd.DataFrame], burdens: dict[str, pd.DataFrame]) -> list[dict]:
    rows = []
    for region in REGIONS:
        L = leads[region]
        B = burdens[region]
        n = SAMPLE_N[region]
        n_tested = len(L)
        n_sig = int(L["sig"].sum())
        n_vmr = len(B)
        n_vmr_sup = int(B["meqtl_supported"].sum())
        rows.append({
            "region": region,
            "n_samples": n,
            "n_cpgs_tested": n_tested,
            "n_cpgs_sig": n_sig,
            "frac_cpgs_sig": n_sig / n_tested if n_tested else np.nan,
            "n_cpgs_sig_per_sqrtN": n_sig / np.sqrt(n),
            "n_vmrs": n_vmr,
            "n_vmrs_meqtl_supported": n_vmr_sup,
            "frac_vmrs_supported": n_vmr_sup / n_vmr if n_vmr else np.nan,
            "n_vmrs_supported_per_sqrtN": n_vmr_sup / np.sqrt(n),
            "note": "N-normalized rates are descriptive; formal donor-downsample remapping not run",
        })
    # Shared-CpG discovery contrast (equal phenotype universe)
    for a, b in combinations(REGIONS, 2):
        ia = set(leads[a]["phenotype_id"])
        ib = set(leads[b]["phenotype_id"])
        shared = ia & ib
        sa = leads[a].set_index("phenotype_id").loc[list(shared), "sig"]
        sb = leads[b].set_index("phenotype_id").loc[list(shared), "sig"]
        rows.append({
            "region": f"shared_cpgs_{a}_vs_{b}",
            "n_samples": f"{SAMPLE_N[a]};{SAMPLE_N[b]}",
            "n_cpgs_tested": len(shared),
            "n_cpgs_sig": f"{int(sa.sum())};{int(sb.sum())}",
            "frac_cpgs_sig": f"{sa.mean():.6f};{sb.mean():.6f}",
            "note": "Equal CpG set; residual sample-size power differences remain",
        })
    return rows


def burden_gradient_across_regions() -> list[dict]:
    rows = []
    for region in REGIONS:
        p = PHASE2 / region / "burden_model_results.tsv"
        if not p.exists():
            continue
        df = pd.read_csv(p, sep="\t")
        for _, r in df.iterrows():
            rows.append({
                "region": region,
                "model": r.get("model", ""),
                "coef_predictability": r.get("coef_predictability", ""),
                "se_predictability": r.get("se_predictability", ""),
                "pval_predictability": r.get("pval_predictability", ""),
                "n_vmrs": r.get("n_vmrs", ""),
            })
    return rows


def write_shared_vmr_table(base: pd.DataFrame) -> None:
    cols = [f"supported_{r}" for r in REGIONS]
    out = base.dropna(subset=cols).copy()
    out["n_regions_supported"] = out[cols].astype(bool).sum(axis=1)
    out["supported_all3"] = out["n_regions_supported"].eq(3)
    keep = ["vmr_coord_id", "n_regions_supported", "supported_all3"] + cols
    keep.extend([f"vmr_id_{r}" for r in REGIONS])
    for r in REGIONS:
        keep.append(f"pred_{r}")
        keep.append(f"prop_{r}")
    out[keep].to_csv(OUTDIR / "vmr_cross_region_support.tsv.gz", sep="\t", index=False, compression="gzip")


def main() -> None:
    args = parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    np.random.seed(args.seed)

    print("Loading burden tables...")
    burdens = {r: load_burden(r) for r in REGIONS}
    print("Loading lead cis-meQTL tables...")
    leads = {r: load_leads(r, args.fdr) for r in REGIONS}

    vmr_rows, base = pairwise_vmr_sharing(burdens)
    write_tsv(OUTDIR / "cross_region_vmr_sharing.tsv", vmr_rows)
    write_shared_vmr_table(base)

    cpg_rows = pairwise_cpg_concordance(leads, args.fdr)
    write_tsv(OUTDIR / "cross_region_cpg_concordance.tsv", cpg_rows)

    disc = discovery_power_summary(leads, burdens)
    pd.DataFrame(disc).to_csv(OUTDIR / "discovery_power_summary.tsv", sep="\t", index=False)

    grad = burden_gradient_across_regions()
    write_tsv(OUTDIR / "burden_gradient_by_region.tsv", grad)

    donors = load_shared_donors()
    donor_row = {
        "n_caudate": donors["n_by_region"]["caudate"],
        "n_dlpfc": donors["n_by_region"]["dlpfc"],
        "n_hippocampus": donors["n_by_region"]["hippocampus"],
        "n_shared_caudate_dlpfc": donors["n_shared_caudate_dlpfc"],
        "n_shared_caudate_hippocampus": donors["n_shared_caudate_hippocampus"],
        "n_shared_dlpfc_hippocampus": donors["n_shared_dlpfc_hippocampus"],
        "n_shared_all3": donors["n_shared_all3"],
    }
    write_tsv(OUTDIR / "shared_donor_counts.tsv", [donor_row])

    # Pending / done markers for Experiment 3 (preserve done status if already run)
    pending_path = OUTDIR / "pending_analyses.tsv"
    prior = {}
    if pending_path.exists():
        try:
            prev = pd.read_csv(pending_path, sep="\t")
            for _, r in prev.iterrows():
                prior[str(r["analysis"])] = {
                    "status": str(r.get("status", "pending")),
                    "reason": str(r.get("reason", "")),
                }
        except Exception:
            prior = {}

    def _row(name: str, default_reason: str) -> dict:
        if name in prior and prior[name]["status"] == "done":
            return {
                "analysis": name,
                "status": "done",
                "reason": prior[name]["reason"],
                "n_shared_all3": donors["n_shared_all3"],
            }
        return {
            "analysis": name,
            "status": "pending",
            "reason": default_reason,
            "n_shared_all3": donors["n_shared_all3"],
        }

    write_tsv(OUTDIR / "pending_analyses.tsv", [
        _row(
            "caudate_donor_downsample_remap",
            "Requires N-matched caudate lead-SNP retention (step_3_downsample.sh)",
        ),
        _row(
            "shared_donor_genotype_by_region",
            "Requires region×genotype models on repeated donors (step_4_gxregion.sh)",
        ),
    ])

    print(f"Wrote Phase 4 cross-region outputs under {OUTDIR}")
    print(pd.DataFrame(vmr_rows)[["contrast", "n_shared_tested", "n_supported_both", "jaccard_support"]].to_string(index=False))
    print(pd.DataFrame(cpg_rows)[[
        "contrast", "n_shared_tested", "n_sig_both", "direction_concordance_both_sig", "pearson_z_either_sig"
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
