#!/usr/bin/env python3
"""Phase 7 Analysis 5 (partial): cross-region SCZ risk–CpG concordance.

Compares targeted risk-variant–CpG tests across caudate, DLPFC, and hippocampus.
Does not re-run TensorQTL. Formal caudate donor-downsampling remapping and
shared-donor genotype×region models are flagged as pending.
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
P7 = PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m"
REGIONS = ["caudate", "dlpfc", "hippocampus"]
SAMPLE_N = {"caudate": 153, "dlpfc": 111, "hippocampus": 116}
SEED = 20260801


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=str(P7 / "cross_region"))
    p.add_argument("--fdr", type=float, default=0.05)
    return p.parse_args()


def load_region(region: str) -> pd.DataFrame:
    path = P7 / region / "risk_variant_cpg_meqtl.tsv.gz"
    if not path.exists():
        raise SystemExit(f"Missing {path}")
    df = pd.read_csv(path, sep="\t")
    df["region"] = region
    df["index_snp"] = df["index_snp"].astype(str)
    df["phenotype_id"] = df["phenotype_id"].astype(str)
    df["variant_id"] = df["variant_id"].astype(str)
    df["vmr_id"] = df["vmr_id"].astype(str)
    for c in ["beta", "se", "tstat", "pval_nominal", "qval", "r2", "maf", "n", "local_predictability"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["sig"] = df["significant_fdr"].astype(bool) if "significant_fdr" in df.columns else df["qval"].le(0.05)
    df["z"] = df["beta"] / df["se"].replace(0, np.nan)
    df["z_sqrtn"] = df["z"] / np.sqrt(df["n"].replace(0, np.nan))
    df["pair_id"] = df["index_snp"] + "::" + df["phenotype_id"]
    return df


def pairwise_concordance(a: pd.DataFrame, b: pd.DataFrame, ra: str, rb: str) -> dict:
    m = a.merge(b, on="pair_id", suffixes=("_a", "_b"), how="inner")
    if m.empty:
        return {
            "region_a": ra, "region_b": rb, "n_shared_tested_pairs": 0,
            "error": "no_overlapping_pairs",
        }
    both_sig = m["sig_a"] & m["sig_b"]
    either_sig = m["sig_a"] | m["sig_b"]
    conc_both = np.nan
    if both_sig.any():
        conc_both = float(np.mean(np.sign(m.loc[both_sig, "beta_a"]) == np.sign(m.loc[both_sig, "beta_b"])))
    conc_either = np.nan
    if either_sig.any():
        conc_either = float(np.mean(np.sign(m.loc[either_sig, "beta_a"]) == np.sign(m.loc[either_sig, "beta_b"])))

    def _corr(x, y):
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 5:
            return np.nan, np.nan, int(ok.sum())
        r, p = pearsonr(x[ok], y[ok])
        return float(r), float(p), int(ok.sum())

    r_z, p_z, n_z = _corr(m.loc[either_sig, "z_a"].to_numpy(), m.loc[either_sig, "z_b"].to_numpy()) if either_sig.any() else (np.nan, np.nan, 0)
    r_b, p_b, n_b = _corr(m["beta_a"].to_numpy(), m["beta_b"].to_numpy())
    r_zs, p_zs, n_zs = _corr(
        m.loc[either_sig, "z_sqrtn_a"].to_numpy(),
        m.loc[either_sig, "z_sqrtn_b"].to_numpy(),
    ) if either_sig.any() else (np.nan, np.nan, 0)

    # Locus / VMR sharing among significant sets
    loci_a = set(a.loc[a["sig"], "locus_id"].dropna().astype(str))
    loci_b = set(b.loc[b["sig"], "locus_id"].dropna().astype(str))
    vmr_a = set(a.loc[a["sig"], "vmr_id"])
    vmr_b = set(b.loc[b["sig"], "vmr_id"])
    snp_a = set(a.loc[a["sig"], "index_snp"])
    snp_b = set(b.loc[b["sig"], "index_snp"])

    def jacc(x, y):
        u = x | y
        return float(len(x & y) / len(u)) if u else np.nan

    return {
        "region_a": ra,
        "region_b": rb,
        "n_a": SAMPLE_N[ra],
        "n_b": SAMPLE_N[rb],
        "n_shared_tested_pairs": int(len(m)),
        "n_both_sig": int(both_sig.sum()),
        "n_either_sig": int(either_sig.sum()),
        "n_sig_only_a": int((m["sig_a"] & ~m["sig_b"]).sum()),
        "n_sig_only_b": int((~m["sig_a"] & m["sig_b"]).sum()),
        "direction_concordance_both_sig": conc_both,
        "direction_concordance_either_sig": conc_either,
        "pearson_z_either_sig": r_z,
        "pearson_z_either_sig_p": p_z,
        "n_pairs_z_corr": n_z,
        "pearson_beta_all_shared": r_b,
        "pearson_beta_all_shared_p": p_b,
        "n_pairs_beta_corr": n_b,
        "pearson_z_sqrtn_either_sig": r_zs,
        "pearson_z_sqrtn_either_sig_p": p_zs,
        "n_pairs_z_sqrtn_corr": n_zs,
        "jaccard_sig_loci": jacc(loci_a, loci_b),
        "n_shared_sig_loci": int(len(loci_a & loci_b)),
        "n_sig_loci_a": int(len(loci_a)),
        "n_sig_loci_b": int(len(loci_b)),
        "jaccard_sig_vmrs": jacc(vmr_a, vmr_b),
        "n_shared_sig_vmrs": int(len(vmr_a & vmr_b)),
        "jaccard_sig_index_snps": jacc(snp_a, snp_b),
        "n_shared_sig_index_snps": int(len(snp_a & snp_b)),
    }


def locus_region_matrix(dfs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    all_loci = set()
    for df in dfs.values():
        all_loci |= set(df["locus_id"].dropna().astype(str))
    for locus in sorted(all_loci, key=lambda x: (len(str(x)), str(x))):
        row = {"locus_id": locus}
        n_regions_sig = 0
        n_regions_tested = 0
        betas = {}
        for r, df in dfs.items():
            sub = df[df["locus_id"].astype(str) == str(locus)]
            row[f"n_tests_{r}"] = int(len(sub))
            row[f"n_sig_{r}"] = int(sub["sig"].sum()) if len(sub) else 0
            row[f"min_pval_{r}"] = float(sub["pval_nominal"].min()) if len(sub) else np.nan
            row[f"min_qval_{r}"] = float(sub["qval"].min()) if len(sub) else np.nan
            if len(sub):
                n_regions_tested += 1
                # strongest absolute beta among sig, else among all
                use = sub[sub["sig"]] if sub["sig"].any() else sub
                i = use["pval_nominal"].idxmin()
                row[f"best_index_snp_{r}"] = use.at[i, "index_snp"]
                row[f"best_beta_{r}"] = float(use.at[i, "beta"])
                row[f"best_se_{r}"] = float(use.at[i, "se"])
                row[f"best_z_{r}"] = float(use.at[i, "z"]) if np.isfinite(use.at[i, "z"]) else np.nan
                betas[r] = row[f"best_beta_{r}"]
                if row[f"n_sig_{r}"] > 0:
                    n_regions_sig += 1
            else:
                row[f"best_index_snp_{r}"] = pd.NA
                row[f"best_beta_{r}"] = np.nan
                row[f"best_se_{r}"] = np.nan
                row[f"best_z_{r}"] = np.nan
        row["n_regions_tested"] = n_regions_tested
        row["n_regions_with_sig"] = n_regions_sig
        # direction concordance among regions with finite best beta
        signs = [np.sign(betas[r]) for r in betas if np.isfinite(betas[r]) and betas[r] != 0]
        row["n_regions_with_nonzero_beta"] = len(signs)
        row["all_regions_same_beta_sign"] = bool(len(signs) >= 2 and len(set(signs)) == 1)
        # caudate-only significance pattern (not selectivity proof)
        row["sig_caudate_only"] = bool(
            row.get("n_sig_caudate", 0) > 0
            and row.get("n_sig_dlpfc", 0) == 0
            and row.get("n_sig_hippocampus", 0) == 0
        )
        rows.append(row)
    return pd.DataFrame(rows)


def discovery_power_summary(dfs: dict[str, pd.DataFrame]) -> list[dict]:
    rows = []
    by: dict[str, dict] = {}
    for r, df in dfs.items():
        n = SAMPLE_N[r]
        n_sig = int(df["sig"].sum())
        n_loci = int(df.loc[df["sig"], "locus_id"].nunique())
        rec = {
            "row_type": "region",
            "region": r,
            "n_samples": n,
            "n_tests": int(len(df)),
            "n_sig_pairs": n_sig,
            "n_sig_loci": n_loci,
            "sig_pairs_per_sqrtN": n_sig / np.sqrt(n),
            "sig_loci_per_sqrtN": n_loci / np.sqrt(n),
            "raw_sig_pair_ratio_vs_other": np.nan,
            "sqrtN_sig_pair_ratio_vs_other": np.nan,
            "raw_sig_loci_ratio_vs_other": np.nan,
            "sqrtN_sig_loci_ratio_vs_other": np.nan,
            "note": "sqrtN-normalized discovery rate; not a formal downsample remap",
        }
        by[r] = rec
        rows.append(rec)
    for other in ["dlpfc", "hippocampus"]:
        rows.append({
            "row_type": "contrast",
            "region": f"caudate_vs_{other}",
            "n_samples": SAMPLE_N["caudate"],
            "n_tests": np.nan,
            "n_sig_pairs": np.nan,
            "n_sig_loci": np.nan,
            "sig_pairs_per_sqrtN": np.nan,
            "sig_loci_per_sqrtN": np.nan,
            "raw_sig_pair_ratio_vs_other": by["caudate"]["n_sig_pairs"] / max(by[other]["n_sig_pairs"], 1),
            "sqrtN_sig_pair_ratio_vs_other": (
                by["caudate"]["sig_pairs_per_sqrtN"] / max(by[other]["sig_pairs_per_sqrtN"], 1e-12)
            ),
            "raw_sig_loci_ratio_vs_other": by["caudate"]["n_sig_loci"] / max(by[other]["n_sig_loci"], 1),
            "sqrtN_sig_loci_ratio_vs_other": (
                by["caudate"]["sig_loci_per_sqrtN"] / max(by[other]["sig_loci_per_sqrtN"], 1e-12)
            ),
            "note": "ratio near 1 after sqrtN suggests N-driven discovery differences",
        })
    return rows


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dfs = {r: load_region(r) for r in REGIONS}
    pair_rows = []
    for ra, rb in combinations(REGIONS, 2):
        pair_rows.append(pairwise_concordance(dfs[ra], dfs[rb], ra, rb))
    write_tsv(outdir / "pairwise_concordance.tsv", pair_rows)

    locus = locus_region_matrix(dfs)
    locus.to_csv(outdir / "locus_by_region.tsv.gz", sep="\t", index=False, compression="gzip")

    # Shared significant pairs detail
    shared = None
    for r, df in dfs.items():
        sig = df[df["sig"]][[
            "pair_id", "index_snp", "phenotype_id", "locus_id", "vmr_id",
            "beta", "se", "z", "qval", "local_predictability",
        ]].copy()
        sig = sig.rename(columns={
            "index_snp": f"index_snp_{r}",
            "phenotype_id": f"phenotype_id_{r}",
            "beta": f"beta_{r}",
            "se": f"se_{r}",
            "z": f"z_{r}",
            "qval": f"qval_{r}",
            "vmr_id": f"vmr_id_{r}",
            "local_predictability": f"pred_{r}",
            "locus_id": f"locus_id_{r}",
        })
        if shared is None:
            shared = sig
        else:
            shared = shared.merge(sig, on="pair_id", how="outer")
    assert shared is not None
    # Harmonize identity columns
    shared["index_snp"] = shared[[c for c in shared.columns if c.startswith("index_snp_")]].bfill(axis=1).iloc[:, 0]
    shared["phenotype_id"] = shared[[c for c in shared.columns if c.startswith("phenotype_id_")]].bfill(axis=1).iloc[:, 0]
    shared["n_regions_sig"] = shared[[c for c in shared.columns if c.startswith("qval_")]].notna().sum(axis=1)
    shared.to_csv(outdir / "sig_pairs_by_region.tsv.gz", sep="\t", index=False, compression="gzip")

    write_tsv(outdir / "discovery_power_summary.tsv", discovery_power_summary(dfs))

    # Claim snapshot
    n_loci_multi = int((locus["n_regions_with_sig"] >= 2).sum())
    n_loci_all3 = int((locus["n_regions_with_sig"] == 3).sum())
    n_caud_only = int(locus["sig_caudate_only"].sum())
    write_tsv(outdir / "cross_region_claim_snapshot.tsv", [{
        "n_loci_sig_any_region": int((locus["n_regions_with_sig"] >= 1).sum()),
        "n_loci_sig_ge2_regions": n_loci_multi,
        "n_loci_sig_all3_regions": n_loci_all3,
        "n_loci_sig_caudate_only_pattern": n_caud_only,
        "caudate_dlpfc_dir_conc_both_sig": next(
            r["direction_concordance_both_sig"] for r in pair_rows
            if r["region_a"] == "caudate" and r["region_b"] == "dlpfc"
        ),
        "caudate_hip_dir_conc_both_sig": next(
            r["direction_concordance_both_sig"] for r in pair_rows
            if r["region_a"] == "caudate" and r["region_b"] == "hippocampus"
        ),
        "pending_caudate_downsample_remap": True,
        "pending_shared_donor_Gxregion": True,
        "interpretation": (
            "Significance-only caudate enrichment is not sufficient for caudate-selectivity; "
            "use effect-size concordance and pending downsample/GxR analyses."
        ),
    }])
    print(f"Wrote cross-region concordance under {outdir}")
    print(f"loci sig in ≥2 regions={n_loci_multi}; all3={n_loci_all3}; caudate-only pattern={n_caud_only}")


if __name__ == "__main__":
    main()
