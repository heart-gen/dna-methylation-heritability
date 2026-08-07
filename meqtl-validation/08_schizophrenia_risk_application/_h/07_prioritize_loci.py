#!/usr/bin/env python3
"""Prioritize ≤5 illustrative SCZ-risk loci for multi-omic follow-up.

Prespecified ranking uses caudate as the primary region:
fine-mapping support, internal risk–CpG meQTL, VMR predictability,
transcriptional coupling, mappability, and cross-region support.
Does not use diagnosis or mediation results.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
P7 = PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m"
RAW = P7 / "raw"
TECH = PROJECT / "meqtl-validation/07_repeat_mappability_sensitivity/_m"
MAX_LOCI = 5
SEED = 20260801


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--primary-region", default="caudate")
    p.add_argument("--max-loci", type=int, default=MAX_LOCI)
    p.add_argument("--outdir", default=str(P7 / "prioritized"))
    p.add_argument("--umap-min", type=float, default=0.8)
    return p.parse_args()


def load_finemap() -> pd.DataFrame:
    path = RAW / "pgc3_finemap_st11a_credible_sets.tsv"
    if not path.exists():
        raise SystemExit(f"Missing {path}; export FINEMAP ST11a first")
    fm = pd.read_csv(path, sep="\t")
    fm["rsid"] = fm["rsid"].astype(str)
    fm["index_snp"] = fm["index_snp"].astype(str)
    fm["finemap_posterior_probability"] = pd.to_numeric(
        fm["finemap_posterior_probability"], errors="coerce"
    )
    # Per index SNP: max PP among credible-set members; whether index itself is in CS
    by_index = fm.groupby("index_snp", as_index=False).agg(
        max_cs_pip=("finemap_posterior_probability", "max"),
        n_cs_snps=("rsid", "nunique"),
        finemap_genes=("gene_symbol", lambda s: ",".join(sorted({str(x) for x in s.dropna() if str(x) and str(x) != "nan"}))[:500]),
    )
    # Best PP for the index SNP rsid itself when present as a CS member
    self_pp = (
        fm[fm["rsid"] == fm["index_snp"]]
        .groupby("index_snp", as_index=False)["finemap_posterior_probability"]
        .max()
        .rename(columns={"finemap_posterior_probability": "index_snp_cs_pip"})
    )
    by_index = by_index.merge(self_pp, on="index_snp", how="left")
    by_index["index_in_credible_set"] = by_index["index_snp_cs_pip"].notna()
    # Small CS flag from ST11c
    small = RAW / "pgc3_finemap_st11c_small_credible_sets.tsv"
    if small.exists():
        sm = pd.read_csv(small, sep="\t")
        by_index["in_small_credible_set_table"] = by_index["index_snp"].isin(set(sm["index_snp"].astype(str)))
    else:
        by_index["in_small_credible_set_table"] = False
    return by_index


def load_tech(region: str) -> pd.DataFrame:
    path = TECH / region / "vmr_technical_annotations.tsv"
    if not path.exists():
        return pd.DataFrame(columns=["vmr_id", "umap_k24_mean", "segdup_frac", "high_mappability"])
    t = pd.read_csv(path, sep="\t")
    if "task_id" in t.columns:
        t = t.copy()
        t["vmr_id"] = pd.to_numeric(t["task_id"], errors="coerce")
        t = t.dropna(subset=["vmr_id"])
        t["vmr_id"] = t["vmr_id"].astype(int).astype(str)
    elif "vmr_id" in t.columns:
        t["vmr_id"] = t["vmr_id"].astype(str)
    else:
        return pd.DataFrame(columns=["vmr_id", "umap_k24_mean", "segdup_frac", "high_mappability"])
    keep = [c for c in ["vmr_id", "umap_k24_mean", "segdup_frac", "high_mappability", "line_l1_frac"] if c in t.columns]
    return t[keep].drop_duplicates("vmr_id", keep="first")


def load_tx(region: str) -> pd.DataFrame:
    path = P7 / region / "scz_meqtl_tx_coupled_vmrs.tsv.gz"
    if not path.exists():
        return pd.DataFrame(columns=["vmr_id", "tx_expression", "tx_psi"])
    tx = pd.read_csv(path, sep="\t")
    tx["vmr_id"] = tx["vmr_id"].astype(str)
    out = tx.groupby("vmr_id", as_index=False).agg(
        tx_modalities=("modality", lambda s: ",".join(sorted(set(map(str, s))))),
    )
    out["tx_expression"] = out["tx_modalities"].str.contains("expression", na=False)
    out["tx_psi"] = out["tx_modalities"].str.contains("psi", na=False)
    out["tx_any"] = out["tx_expression"] | out["tx_psi"]
    return out


def load_cross_region() -> pd.DataFrame:
    path = P7 / "cross_region" / "locus_by_region.tsv.gz"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", compression="infer")


def score_loci(primary: str, umap_min: float) -> pd.DataFrame:
    res = pd.read_csv(P7 / primary / "risk_variant_cpg_meqtl.tsv.gz", sep="\t")
    res["index_snp"] = res["index_snp"].astype(str)
    res["vmr_id"] = res["vmr_id"].astype(str)
    res["sig"] = res["significant_fdr"].astype(bool) if "significant_fdr" in res.columns else res["qval"].le(0.05)
    sig = res[res["sig"]].copy()
    if sig.empty:
        raise SystemExit(f"No FDR-significant risk–CpG pairs in {primary}")

    # Best pair per locus (primary ranking unit = locus_id; fall back to index_snp)
    sig["locus_key"] = sig["locus_id"].fillna(sig["index_snp"]).astype(str)
    sig = sig.sort_values(["qval", "pval_nominal"])
    best = sig.groupby("locus_key", as_index=False).first()

    # Aggregate locus stats
    agg = sig.groupby("locus_key", as_index=False).agg(
        n_sig_pairs=("phenotype_id", "size"),
        n_sig_cpgs=("phenotype_id", "nunique"),
        n_sig_vmrs=("vmr_id", "nunique"),
        n_index_snps=("index_snp", "nunique"),
        min_qval=("qval", "min"),
        min_pval=("pval_nominal", "min"),
        max_abs_beta=("beta", lambda s: float(np.nanmax(np.abs(pd.to_numeric(s, errors="coerce"))))),
        max_predictability=("local_predictability", "max"),
        mean_predictability=("local_predictability", "mean"),
        index_snp=("index_snp", "first"),
        locus_id=("locus_id", "first"),
        best_vmr_id=("vmr_id", "first"),
        best_phenotype_id=("phenotype_id", "first"),
        best_variant_id=("variant_id", "first"),
        chrom=("chrom", "first"),
        gwas_p=("gwas_p", "first"),
        gwas_or=("gwas_or", "first"),
    )

    fm = load_finemap()
    agg = agg.merge(fm, on="index_snp", how="left")

    tx = load_tx(primary)
    # TX at any significant VMR in the locus
    sig_vmr_tx = sig[["locus_key", "vmr_id"]].drop_duplicates().merge(tx, on="vmr_id", how="left")
    tx_locus = sig_vmr_tx.groupby("locus_key", as_index=False).agg(
        tx_any=("tx_any", lambda s: bool(pd.Series(s).fillna(False).any())),
        tx_expression=("tx_expression", lambda s: bool(pd.Series(s).fillna(False).any())),
        tx_psi=("tx_psi", lambda s: bool(pd.Series(s).fillna(False).any())),
        tx_vmrs=("vmr_id", lambda s: ",".join(sorted(set(map(str, s))))),
    )
    agg = agg.merge(tx_locus, on="locus_key", how="left")
    for c in ["tx_any", "tx_expression", "tx_psi"]:
        agg[c] = agg[c].fillna(False).astype(bool)

    tech = load_tech(primary)
    best_tech = best[["locus_key", "vmr_id"]].merge(tech, on="vmr_id", how="left")
    agg = agg.merge(
        best_tech.rename(columns={
            "umap_k24_mean": "best_vmr_umap",
            "segdup_frac": "best_vmr_segdup_frac",
            "high_mappability": "best_vmr_high_mappability",
            "line_l1_frac": "best_vmr_line_l1_frac",
            "vmr_id": "scored_vmr_id",
        }),
        on="locus_key",
        how="left",
    )

    xref = load_cross_region()
    if not xref.empty and "locus_id" in xref.columns:
        xref["locus_id"] = xref["locus_id"].astype(str)
        agg["locus_id"] = agg["locus_id"].astype(str)
        agg = agg.merge(
            xref[[
                c for c in [
                    "locus_id", "n_regions_with_sig", "n_regions_tested",
                    "all_regions_same_beta_sign", "sig_caudate_only",
                    "n_sig_dlpfc", "n_sig_hippocampus", "n_sig_caudate",
                ] if c in xref.columns
            ]],
            on="locus_id",
            how="left",
        )
    else:
        agg["n_regions_with_sig"] = np.nan

    # --- Scoring (transparent additive rubric) ---
    s = pd.DataFrame(index=agg.index)
    s["score_sig_meqtl"] = 2.0  # all candidates already FDR sig in primary
    s["score_finemap_index_in_cs"] = np.where(agg["index_in_credible_set"].fillna(False), 2.0, 0.0)
    s["score_finemap_max_pip"] = np.select(
        [
            agg["max_cs_pip"].fillna(0) >= 0.5,
            agg["max_cs_pip"].fillna(0) >= 0.1,
            agg["max_cs_pip"].notna(),
        ],
        [3.0, 2.0, 1.0],
        default=0.0,
    )
    s["score_small_cs"] = np.where(agg.get("in_small_credible_set_table", False), 1.0, 0.0)
    # predictability: high if max >= region 80th among candidates
    pred_q80 = agg["max_predictability"].quantile(0.8)
    pred_q50 = agg["max_predictability"].quantile(0.5)
    s["score_predictability"] = np.select(
        [agg["max_predictability"] >= pred_q80, agg["max_predictability"] >= pred_q50],
        [2.0, 1.0],
        default=0.0,
    )
    s["score_tx_expression"] = np.where(agg["tx_expression"], 2.0, 0.0)
    s["score_tx_psi"] = np.where(agg["tx_psi"], 2.0, 0.0)
    s["score_tx_both"] = np.where(agg["tx_expression"] & agg["tx_psi"], 1.0, 0.0)
    umap = pd.to_numeric(agg["best_vmr_umap"], errors="coerce")
    seg = pd.to_numeric(agg["best_vmr_segdup_frac"], errors="coerce").fillna(0)
    s["score_mappability"] = np.where(umap.fillna(0) >= umap_min, 1.0, 0.0)
    s["score_low_segdup"] = np.where(seg <= 0.05, 1.0, 0.0)
    nreg = pd.to_numeric(agg["n_regions_with_sig"], errors="coerce").fillna(1)
    s["score_multiregion"] = np.select([nreg >= 3, nreg >= 2], [2.0, 1.0], default=0.0)
    # Prefer not exclusively caudate-only pattern when other regions also tested
    caud_only = (
        agg["sig_caudate_only"].fillna(False).astype(bool)
        if "sig_caudate_only" in agg.columns
        else pd.Series(False, index=agg.index)
    )
    s["score_not_caudate_only_pattern"] = np.where(caud_only & (nreg == 1), 0.0, 0.5)
    # Stronger meQTL (lower q)
    q = pd.to_numeric(agg["min_qval"], errors="coerce")
    s["score_meqtl_strength"] = np.select(
        [q <= 1e-4, q <= 1e-3, q <= 0.01],
        [2.0, 1.0, 0.5],
        default=0.0,
    )
    # Gene annotation present from FINEMAP
    s["score_finemap_gene"] = np.where(
        agg["finemap_genes"].fillna("").astype(str).str.len().gt(0), 0.5, 0.0
    )

    agg["priority_score"] = s.sum(axis=1)
    for c in s.columns:
        agg[c] = s[c]
    agg["pred_q80_threshold"] = pred_q80
    agg["pred_q50_threshold"] = pred_q50
    return agg.sort_values(
        ["priority_score", "tx_any", "max_predictability", "min_qval"],
        ascending=[False, False, False, True],
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ranked = score_loci(args.primary_region, args.umap_min)
    ranked.to_csv(outdir / "locus_priority_ranked_all.tsv.gz", sep="\t", index=False, compression="gzip")

    top = ranked.head(args.max_loci).copy()
    top["rank"] = np.arange(1, len(top) + 1)
    for c in [
        "index_in_credible_set", "in_small_credible_set_table",
        "tx_any", "tx_expression", "tx_psi", "sig_caudate_only",
        "all_regions_same_beta_sign",
    ]:
        if c in top.columns:
            top[c] = top[c].fillna(False).astype(bool)
    # Compact presentation columns first
    front = [
        "rank", "priority_score", "locus_id", "index_snp", "chrom",
        "n_sig_pairs", "n_sig_vmrs", "min_qval", "max_abs_beta",
        "max_predictability", "best_vmr_id", "best_phenotype_id", "best_variant_id",
        "index_in_credible_set", "index_snp_cs_pip", "max_cs_pip", "in_small_credible_set_table",
        "finemap_genes", "tx_any", "tx_expression", "tx_psi",
        "best_vmr_umap", "best_vmr_segdup_frac", "n_regions_with_sig",
        "n_sig_caudate", "n_sig_dlpfc", "n_sig_hippocampus", "sig_caudate_only",
        "gwas_p", "gwas_or",
    ]
    front = [c for c in front if c in top.columns]
    rest = [c for c in top.columns if c not in front]
    top[front + rest].to_csv(outdir / "prioritized_loci.tsv", sep="\t", index=False)

    write_tsv(outdir / "prioritization_summary.tsv", [{
        "primary_region": args.primary_region,
        "n_candidate_loci": int(len(ranked)),
        "n_prioritized": int(len(top)),
        "n_prioritized_with_tx": int(top["tx_any"].sum()) if "tx_any" in top.columns else 0,
        "n_prioritized_in_credible_set": int(top["index_in_credible_set"].fillna(False).sum()),
        "n_prioritized_multiregion": int((pd.to_numeric(top.get("n_regions_with_sig"), errors="coerce").fillna(0) >= 2).sum()),
        "min_priority_score": float(top["priority_score"].min()) if len(top) else np.nan,
        "max_priority_score": float(top["priority_score"].max()) if len(top) else np.nan,
        "output": str(outdir / "prioritized_loci.tsv"),
    }])

    def _truthy(val) -> bool:
        if isinstance(val, (bool, np.bool_)):
            return bool(val)
        if val is None or (isinstance(val, float) and not np.isfinite(val)):
            return False
        if isinstance(val, str):
            return val.strip().lower() in {"1", "true", "t", "yes"}
        return bool(val)

    # Human-readable rationale per locus
    reasons = []
    for _, r in top.iterrows():
        bits = []
        if _truthy(r.get("index_in_credible_set")):
            pip = r.get("index_snp_cs_pip")
            pip_s = f"{float(pip):.3g}" if pd.notna(pip) else "NA"
            bits.append(f"index in FINEMAP CS (PIP={pip_s})")
        elif _truthy(r.get("in_small_credible_set_table")):
            bits.append("index SNP locus in small FINEMAP CS table (ST11c)")
        if pd.notna(r.get("max_cs_pip")):
            bits.append(f"locus max CS PIP={float(r.get('max_cs_pip')):.3g}")
        bits.append(f"caudate FDR meQTL q={float(r.get('min_qval')):.3g} ({int(r.get('n_sig_pairs', 0))} pairs)")
        bits.append(f"max VMR predictability={float(r.get('max_predictability')):.3g}")
        if _truthy(r.get("tx_expression")):
            bits.append("expression coupling")
        if _truthy(r.get("tx_psi")):
            bits.append("PSI coupling")
        if pd.notna(r.get("n_regions_with_sig")) and float(r.get("n_regions_with_sig")) >= 2:
            bits.append(f"sig in {int(r.get('n_regions_with_sig'))} regions")
        if pd.notna(r.get("best_vmr_umap")):
            bits.append(f"umap={float(r.get('best_vmr_umap')):.3g}")
        genes = r.get("finemap_genes")
        if isinstance(genes, str) and genes and genes not in {"-", "nan"}:
            bits.append(f"genes={genes[:120]}")
        reasons.append({
            "rank": int(r["rank"]),
            "locus_id": r.get("locus_id"),
            "index_snp": r.get("index_snp"),
            "priority_score": float(r["priority_score"]),
            "rationale": "; ".join(bits),
        })
    write_tsv(outdir / "prioritized_loci_rationale.tsv", reasons)
    print(f"Wrote top {len(top)} loci to {outdir / 'prioritized_loci.tsv'}")
    for row in reasons:
        print(f"  #{row['rank']} {row['index_snp']} locus={row['locus_id']} score={row['priority_score']:.1f}")


if __name__ == "__main__":
    main()
