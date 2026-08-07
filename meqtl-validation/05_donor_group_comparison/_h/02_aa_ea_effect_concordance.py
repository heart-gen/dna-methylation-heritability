#!/usr/bin/env python3
"""Phase 4b Experiment 2: AA–EA CpG lead-effect concordance.

Join AA (M3a) and EA (M0) Phase 1 lead tables on shared phenotype_id.
Report direction concordance, slope/z correlations, and identical-lead-SNP rates.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
PHASE1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping"
OUTDIR = PROJECT / "meqtl-validation/05_donor_group_comparison/_m"
REGIONS = ["caudate", "dlpfc", "hippocampus"]
FDR = 0.05


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", default="all", choices=REGIONS + ["all"])
    p.add_argument("--fdr", type=float, default=FDR)
    return p.parse_args()


def lead_path(region: str, population: str) -> Path:
    if population == "AA":
        return PHASE1 / region / "_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz"
    return PHASE1 / region / "_m/tensorqtl" / population / "qc/lead_snp_per_cpg.tsv.gz"


def load_leads(region: str, population: str, fdr: float) -> pd.DataFrame:
    path = lead_path(region, population)
    if not path.exists():
        raise FileNotFoundError(path)
    cols = [
        "phenotype_id", "variant_id", "slope", "slope_se", "qval",
        "pval_nominal", "maf", "tss_distance", "num_var", "ma_count", "ma_samples",
    ]
    peek = pd.read_csv(path, sep="\t", nrows=0)
    usecols = [c for c in cols if c in peek.columns]
    df = pd.read_csv(path, sep="\t", usecols=usecols)
    df["phenotype_id"] = df["phenotype_id"].astype(str)
    df["variant_id"] = df["variant_id"].astype(str)
    for c in ["slope", "slope_se", "qval", "pval_nominal", "maf", "tss_distance", "num_var", "ma_count", "ma_samples"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    # Derive MAF when tensorQTL QC export omitted the maf column
    n_by_pop_region = {
        ("AA", "caudate"): 153, ("AA", "dlpfc"): 111, ("AA", "hippocampus"): 116,
        ("EA", "caudate"): 129, ("EA", "dlpfc"): 55, ("EA", "hippocampus"): 60,
    }
    if "maf" not in df.columns and "ma_count" in df.columns:
        n = n_by_pop_region.get((population, region))
        if n:
            df["maf"] = (df["ma_count"] / (2.0 * n)).clip(upper=0.5)
        elif "ma_samples" in df.columns:
            df["maf"] = (df["ma_count"] / (2.0 * df["ma_samples"].clip(lower=1))).clip(upper=0.5)
    df["sig"] = df["qval"].le(fdr).fillna(False) if "qval" in df.columns else False
    df["z"] = df["slope"] / df["slope_se"].replace(0, np.nan)
    return df


def summarize_region(region: str, fdr: float) -> tuple[dict, pd.DataFrame]:
    aa = load_leads(region, "AA", fdr).add_suffix("_AA")
    aa = aa.rename(columns={"phenotype_id_AA": "phenotype_id"})
    ea = load_leads(region, "EA", fdr).add_suffix("_EA")
    ea = ea.rename(columns={"phenotype_id_EA": "phenotype_id"})
    m = aa.merge(ea, on="phenotype_id", how="inner")
    if m.empty:
        return {"region": region, "n_shared_cpgs": 0, "error": "no_overlap"}, m

    m["same_lead_snp"] = m["variant_id_AA"] == m["variant_id_EA"]
    m["dir_concordant"] = np.sign(m["slope_AA"]) == np.sign(m["slope_EA"])
    # Avoid zero-slope ties
    m.loc[(m["slope_AA"] == 0) | (m["slope_EA"] == 0), "dir_concordant"] = False

    both_sig = m["sig_AA"] & m["sig_EA"]
    either_sig = m["sig_AA"] | m["sig_EA"]
    aa_only = m["sig_AA"] & ~m["sig_EA"]
    ea_only = m["sig_EA"] & ~m["sig_AA"]

    def corr(x, y):
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 20:
            return np.nan, np.nan, np.nan, np.nan
        rp, pp = pearsonr(x[ok], y[ok])
        rs, ps = spearmanr(x[ok], y[ok])
        return float(rp), float(pp), float(rs), float(ps)

    r_slope, p_slope, rs_slope, ps_slope = corr(m["slope_AA"].to_numpy(), m["slope_EA"].to_numpy())
    r_z, p_z, rs_z, ps_z = corr(m["z_AA"].to_numpy(), m["z_EA"].to_numpy())

    # Among identical lead SNP pairs
    same = m[m["same_lead_snp"]]
    if len(same) >= 20:
        r_same, p_same, _, _ = corr(same["slope_AA"].to_numpy(), same["slope_EA"].to_numpy())
        dir_same = float(same["dir_concordant"].mean())
    else:
        r_same, p_same, dir_same = np.nan, np.nan, np.nan

    summary = {
        "region": region,
        "n_shared_cpgs": int(len(m)),
        "n_sig_AA": int(m["sig_AA"].sum()),
        "n_sig_EA": int(m["sig_EA"].sum()),
        "n_sig_both": int(both_sig.sum()),
        "n_sig_AA_only": int(aa_only.sum()),
        "n_sig_EA_only": int(ea_only.sum()),
        "frac_same_lead_snp": float(m["same_lead_snp"].mean()),
        "direction_concordance_all": float(m["dir_concordant"].mean()),
        "direction_concordance_both_sig": float(m.loc[both_sig, "dir_concordant"].mean()) if both_sig.any() else np.nan,
        "direction_concordance_either_sig": float(m.loc[either_sig, "dir_concordant"].mean()) if either_sig.any() else np.nan,
        "direction_concordance_same_lead": dir_same,
        "pearson_slope_all": r_slope,
        "pearson_slope_all_p": p_slope,
        "spearman_slope_all": rs_slope,
        "spearman_slope_all_p": ps_slope,
        "pearson_z_all": r_z,
        "pearson_z_all_p": p_z,
        "pearson_slope_same_lead": r_same,
        "pearson_slope_same_lead_p": p_same,
        "aa_discovery_rate": float(m["sig_AA"].mean()),
        "ea_discovery_rate": float(m["sig_EA"].mean()),
        "aa_ea_discovery_rate_ratio": float(m["sig_AA"].mean() / max(m["sig_EA"].mean(), 1e-12)),
    }
    return summary, m


def main() -> None:
    args = parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    regions = REGIONS if args.region == "all" else [args.region]
    summaries = []
    for region in regions:
        summary, m = summarize_region(region, args.fdr)
        summaries.append(summary)
        if not m.empty:
            reg_out = OUTDIR / region
            reg_out.mkdir(parents=True, exist_ok=True)
            # Compact joined table for downstream matching
            keep = [
                c for c in m.columns
                if c.startswith((
                    "phenotype", "variant", "slope", "slope_se", "z_", "qval", "sig",
                    "maf", "tss", "num_var", "ma_", "dir", "same",
                ))
                or c == "phenotype_id"
            ]
            m[keep].to_csv(
                reg_out / "aa_ea_cpg_effect_concordance.tsv.gz",
                sep="\t", index=False, compression="gzip",
            )
        print(f"{region}: shared={summary.get('n_shared_cpgs')} both_sig={summary.get('n_sig_both')} "
              f"dir_both={summary.get('direction_concordance_both_sig')}")

    write_tsv(OUTDIR / "aa_ea_effect_concordance_summary.tsv", summaries)
    print(pd.DataFrame(summaries)[[
        "region", "n_shared_cpgs", "n_sig_both", "direction_concordance_both_sig",
        "pearson_z_all", "frac_same_lead_snp", "aa_ea_discovery_rate_ratio",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
