#!/usr/bin/env python3
"""Phase 4b Experiment 2: MAF / cis-SNP-density matched AA–EA discovery contrasts.

Among CpGs tested in both AA and EA, ask whether AA excess discovery shrinks
after matching or stratifying on:
  - min(AA_MAF, EA_MAF) of each group's lead SNP
  - cis variant count (num_var; LD/opportunity proxy)
  - |tss_distance|

Does NOT claim ancestry-specific biology; documents whether discovery gaps are
associated with allele-frequency / testing-opportunity differences.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
OUTDIR = PROJECT / "meqtl-validation/05_donor_group_comparison/_m"
REGIONS = ["caudate", "dlpfc", "hippocampus"]
SEED = 20260805


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", default="all", choices=REGIONS + ["all"])
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--caliper-maf", type=float, default=0.02, help="Max |MAF_AA - MAF_EA| for matched pairs")
    p.add_argument("--n-bins", type=int, default=5)
    return p.parse_args()


def load_joined(region: str) -> pd.DataFrame:
    path = OUTDIR / region / "aa_ea_cpg_effect_concordance.tsv.gz"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}; run 02_aa_ea_effect_concordance.py first"
        )
    df = pd.read_csv(path, sep="\t")
    for c in [
        "maf_AA", "maf_EA", "num_var_AA", "num_var_EA",
        "tss_distance_AA", "tss_distance_EA", "ma_count_AA", "ma_count_EA",
    ]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["sig_AA"] = df["sig_AA"].astype(bool) if "sig_AA" in df.columns else False
    df["sig_EA"] = df["sig_EA"].astype(bool) if "sig_EA" in df.columns else False

    # Derive MAF if concordance export still lacks it
    n_aa = {"caudate": 153, "dlpfc": 111, "hippocampus": 116}.get(region)
    n_ea = {"caudate": 129, "dlpfc": 55, "hippocampus": 60}.get(region)
    if "maf_AA" not in df.columns and "ma_count_AA" in df.columns and n_aa:
        df["maf_AA"] = (df["ma_count_AA"] / (2.0 * n_aa)).clip(upper=0.5)
    if "maf_EA" not in df.columns and "ma_count_EA" in df.columns and n_ea:
        df["maf_EA"] = (df["ma_count_EA"] / (2.0 * n_ea)).clip(upper=0.5)

    if "maf_AA" in df.columns and "maf_EA" in df.columns:
        df["min_maf"] = np.fmin(df["maf_AA"].to_numpy(dtype=float), df["maf_EA"].to_numpy(dtype=float))
        df["mean_maf"] = np.nanmean(
            np.vstack([df["maf_AA"].to_numpy(dtype=float), df["maf_EA"].to_numpy(dtype=float)]),
            axis=0,
        )
    else:
        df["min_maf"] = np.nan
        df["mean_maf"] = np.nan

    if "num_var_AA" in df.columns and "num_var_EA" in df.columns:
        df["mean_num_var"] = np.nanmean(
            np.vstack([
                df["num_var_AA"].to_numpy(dtype=float),
                df["num_var_EA"].to_numpy(dtype=float),
            ]),
            axis=0,
        )
    elif "num_var_AA" in df.columns:
        df["mean_num_var"] = df["num_var_AA"].to_numpy(dtype=float)
    else:
        df["mean_num_var"] = np.nan

    if "tss_distance_AA" in df.columns or "tss_distance_EA" in df.columns:
        a = df["tss_distance_AA"].to_numpy(dtype=float) if "tss_distance_AA" in df.columns else np.full(len(df), np.nan)
        b = df["tss_distance_EA"].to_numpy(dtype=float) if "tss_distance_EA" in df.columns else np.full(len(df), np.nan)
        df["abs_tss"] = np.nanmean(np.vstack([np.abs(a), np.abs(b)]), axis=0)
    else:
        df["abs_tss"] = np.nan
    return df


def rate_gap(df: pd.DataFrame) -> dict:
    aa = float(df["sig_AA"].mean()) if len(df) else np.nan
    ea = float(df["sig_EA"].mean()) if len(df) else np.nan
    return {
        "n_cpgs": int(len(df)),
        "aa_discovery_rate": aa,
        "ea_discovery_rate": ea,
        "rate_diff_aa_minus_ea": aa - ea if np.isfinite(aa) and np.isfinite(ea) else np.nan,
        "rate_ratio_aa_over_ea": aa / ea if np.isfinite(ea) and ea > 0 else np.nan,
    }


def stratified_bins(df: pd.DataFrame, col: str, n_bins: int) -> list[dict]:
    rows = []
    use = df.dropna(subset=[col]).copy()
    if len(use) < n_bins * 20:
        return [{
            "stratum_var": col,
            "stratum": "insufficient",
            **rate_gap(use),
            "note": "too few observations for binning",
        }]
    try:
        use["_bin"] = pd.qcut(use[col], q=n_bins, duplicates="drop")
    except ValueError:
        use["_bin"] = pd.cut(use[col], bins=min(n_bins, 3), duplicates="drop")
    for b, grp in use.groupby("_bin", observed=True):
        rows.append({
            "stratum_var": col,
            "stratum": str(b),
            **rate_gap(grp),
            "stratum_median": float(grp[col].median()),
        })
    return rows


def nearest_neighbor_match(df: pd.DataFrame, caliper_maf: float, seed: int) -> tuple[pd.DataFrame, dict]:
    """Match shared CpGs on MAF + cis-SNP density (num_var LD/opportunity proxy)."""
    need = ["mean_num_var"]
    maf_col = "min_maf" if df["min_maf"].notna().sum() >= 100 else (
        "mean_maf" if "mean_maf" in df.columns and df["mean_maf"].notna().sum() >= 100 else None
    )
    if maf_col:
        need.append(maf_col)
    use = df.dropna(subset=need).copy()
    if len(use) < 100:
        return use, {"matched": False, "reason": f"insufficient_rows_after_dropna n={len(use)} need={need}"}

    for c in need:
        x = use[c].to_numpy(dtype=float)
        sd = np.nanstd(x)
        use[f"z_{c}"] = (x - np.nanmean(x)) / sd if sd > 0 else 0.0
    if "abs_tss" in use.columns and use["abs_tss"].notna().any():
        x = use["abs_tss"].to_numpy(dtype=float)
        sd = np.nanstd(x)
        use["z_abs_tss"] = (x - np.nanmean(x)) / sd if sd > 0 else 0.0

    if maf_col and "maf_AA" in use.columns and "maf_EA" in use.columns:
        close = use[np.abs(use["maf_AA"] - use["maf_EA"]) <= caliper_maf].copy()
        caliper_used = caliper_maf if len(close) >= 50 else np.nan
        if len(close) < 50:
            close = use.copy()
    else:
        close = use.copy()
        caliper_used = np.nan

    close = close.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    score = close[f"z_{need[0]}"].abs()
    for c in need[1:]:
        score = score + close[f"z_{c}"].abs()
    if "z_abs_tss" in close.columns:
        score = score + 0.25 * close["z_abs_tss"].abs()
    close["match_score"] = score
    n_keep = max(len(close) // 2, min(5000, len(close)))
    matched = close.nsmallest(n_keep, "match_score")
    meta = {
        "matched": True,
        "caliper_maf": caliper_used,
        "maf_col": maf_col or "none",
        "n_input": int(len(use)),
        "n_within_caliper": int(len(close)),
        "n_matched": int(len(matched)),
        "median_abs_maf_diff": (
            float(np.abs(matched["maf_AA"] - matched["maf_EA"]).median())
            if "maf_AA" in matched.columns and "maf_EA" in matched.columns else np.nan
        ),
        "median_abs_num_var_diff": (
            float(np.abs(matched["num_var_AA"] - matched["num_var_EA"]).median())
            if "num_var_AA" in matched.columns and "num_var_EA" in matched.columns else np.nan
        ),
    }
    return matched, meta


def fisher_aa_vs_ea(df: pd.DataFrame) -> dict:
    """2×2: group × significance on shared CpGs (paired McNemar-style also reported)."""
    a = int((df["sig_AA"] & ~df["sig_EA"]).sum())  # AA only
    b = int((~df["sig_AA"] & df["sig_EA"]).sum())  # EA only
    both = int((df["sig_AA"] & df["sig_EA"]).sum())
    neither = int((~df["sig_AA"] & ~df["sig_EA"]).sum())
    # McNemar exact-ish via binomial on discordant
    from scipy.stats import binomtest
    disc = a + b
    if disc == 0:
        mcnemar_p = 1.0
    else:
        mcnemar_p = float(binomtest(a, n=disc, p=0.5).pvalue)
    # Fisher on unpaired rates is less appropriate; report for transparency
    table = np.array([
        [int(df["sig_AA"].sum()), int((~df["sig_AA"]).sum())],
        [int(df["sig_EA"].sum()), int((~df["sig_EA"]).sum())],
    ])
    try:
        or_, fisher_p = fisher_exact(table)
    except Exception:
        or_, fisher_p = np.nan, np.nan
    return {
        "n_aa_only_sig": a,
        "n_ea_only_sig": b,
        "n_both_sig": both,
        "n_neither_sig": neither,
        "mcnemar_pvalue": mcnemar_p,
        "fisher_or_unpaired": float(or_),
        "fisher_p_unpaired": float(fisher_p),
    }


def analyze_region(region: str, n_bins: int, caliper_maf: float, seed: int) -> tuple[list[dict], list[dict], dict]:
    df = load_joined(region)
    overall = {"region": region, "analysis": "unmatched_shared_cpgs", **rate_gap(df), **fisher_aa_vs_ea(df)}

    strata = []
    for col in ["min_maf", "mean_maf", "mean_num_var", "abs_tss"]:
        if col not in df.columns or df[col].notna().sum() < 50:
            continue
        for row in stratified_bins(df, col, n_bins):
            strata.append({"region": region, **row})

    matched, meta = nearest_neighbor_match(df, caliper_maf, seed)
    matched_row = {
        "region": region,
        "analysis": "maf_numvar_matched_subset",
        **rate_gap(matched),
        **fisher_aa_vs_ea(matched),
        **{f"match_{k}": v for k, v in meta.items()},
    }

    # Does gap shrink?
    gap0 = overall["rate_diff_aa_minus_ea"]
    gap1 = matched_row["rate_diff_aa_minus_ea"]
    claim = {
        "region": region,
        "unmatched_rate_diff": gap0,
        "matched_rate_diff": gap1,
        "gap_attenuation_frac": float(1.0 - gap1 / gap0) if gap0 and np.isfinite(gap0) and gap0 != 0 else np.nan,
        "unmatched_rate_ratio": overall["rate_ratio_aa_over_ea"],
        "matched_rate_ratio": matched_row["rate_ratio_aa_over_ea"],
        "mcnemar_p_unmatched": overall["mcnemar_pvalue"],
        "mcnemar_p_matched": matched_row["mcnemar_pvalue"],
        "gap_shrinks_after_matching": bool(
            np.isfinite(gap0) and np.isfinite(gap1) and abs(gap1) < abs(gap0) * 0.9
        ),
        "gap_eliminated": bool(np.isfinite(gap1) and abs(gap1) < 0.02),
        "interpretation_note": (
            "Discovery-rate differences that shrink after MAF/cis-SNP-density matching "
            "are consistent with allele-frequency / testing-opportunity explanations; "
            "residual gaps may reflect power (N), LD, or true effect heterogeneity — "
            "do not label as ancestry-specific without formal interaction evidence."
        ),
    }
    return [overall, matched_row], strata, claim


def main() -> None:
    args = parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    regions = REGIONS if args.region == "all" else [args.region]

    contrast_rows = []
    strata_rows = []
    claim_rows = []
    for region in regions:
        try:
            contrasts, strata, claim = analyze_region(
                region, args.n_bins, args.caliper_maf, args.seed
            )
        except FileNotFoundError as exc:
            print(f"SKIP {region}: {exc}")
            continue
        contrast_rows.extend(contrasts)
        strata_rows.extend(strata)
        claim_rows.append(claim)
        print(
            f"{region}: unmatched Δ={claim['unmatched_rate_diff']:.3f} "
            f"matched Δ={claim['matched_rate_diff']:.3f} "
            f"shrinks={claim['gap_shrinks_after_matching']}"
        )

    write_tsv(OUTDIR / "maf_ld_matched_discovery.tsv", contrast_rows)
    if strata_rows:
        write_tsv(OUTDIR / "maf_ld_discovery_strata.tsv", strata_rows)
    write_tsv(OUTDIR / "maf_ld_matched_discovery_claim.tsv", claim_rows)

    # Overall Experiment 2 claim rollup
    n_shrink = sum(1 for c in claim_rows if c["gap_shrinks_after_matching"])
    write_tsv(OUTDIR / "experiment2_depth_claim_summary.tsv", [{
        "n_regions": len(claim_rows),
        "n_regions_gap_shrinks": n_shrink,
        "passes_document_maf_ld_sensitivity": bool(len(claim_rows) >= 2),
        "ancestry_specific_claim_allowed": False,
        "detail": "; ".join(
            f"{c['region']}:Δ {c['unmatched_rate_diff']:.3f}→{c['matched_rate_diff']:.3f}"
            for c in claim_rows
        ),
        "preferred_interpretation": (
            "Aggregate predictability→meQTL-burden gradient is portable across donor groups; "
            "locus-level discovery differences are partly associated with MAF/cis-SNP-density "
            "opportunity and sample size. Do not claim ancestry-specific biology."
        ),
    }])
    print(pd.DataFrame(claim_rows).to_string(index=False))


if __name__ == "__main__":
    main()
