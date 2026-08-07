#!/usr/bin/env python3
"""Phase 4 donor-group comparison.

Primary discovery meQTL remains AA (locked M3a). EA stratified CpG meQTL +
burden are complete under tensorqtl/EA/ and 02_vmr_meqtl_burden/_m/EA/.

This module:
  1. Records EA meQTL readiness status
  2. Tests AA vs EA local-predictability portability (elastic-net summaries)
  3. Summarizes AA/EA burden-gradient coefficients across regions
  4. Characterizes cross-group-concordant high-predictability VMRs

Experiment 2 depth (concordance + MAF/LD matching) lives in step_2.sh.

Do not label AA–EA differences as ancestry-specific without formal interaction
evidence and adequate matched meQTL power.
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
ENET = PROJECT / "heritability" / "elastic_net_model" / "all_individuals"
PHASE2 = PROJECT / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m"
OUTDIR = PROJECT / "meqtl-validation" / "05_donor_group_comparison" / "_m"
OVERLAP = PROJECT / "inputs" / "data_dictionary" / "_m" / "sample_overlap.tsv"
REGIONS = ["caudate", "dlpfc", "hippocampus"]
HIGH_PRED_Q = 0.75


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", default="all", choices=REGIONS + ["all"])
    p.add_argument("--high-pred-quantile", type=float, default=HIGH_PRED_Q)
    return p.parse_args()


def coord_id(chrom, start, end) -> str:
    c = str(chrom).replace("chr", "")
    return f"{c}:{int(start)}-{int(end)}"


def load_predictability(region: str, race: str) -> pd.DataFrame:
    path = ENET / region / "_m" / f"{region}_summary_elastic-net_{race}.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t")
    df["coord_id"] = [
        coord_id(c, s, e) for c, s, e in zip(df["chrom"], df["start"], df["end"])
    ]
    df["h2_unscaled"] = pd.to_numeric(df["h2_unscaled"], errors="coerce")
    df["race"] = race
    return df[["coord_id", "task_id", "h2_unscaled", "num_snps", "race", "chrom", "start", "end"]]


def ea_meqtl_status() -> list[dict]:
    ov = pd.read_csv(OVERLAP, sep="\t")
    rows = []
    for region in REGIONS:
        sub = ov[(ov["race"] == "EA") & (ov["region"] == region)]
        if sub.empty:
            rows.append({
                "region": region,
                "donor_group": "EA",
                "status": "blocked",
                "reason": "missing sample_overlap row",
            })
            continue
        r = sub.iloc[0]
        geno_ok = float(r.get("n_with_genotype_id_match", 0) or 0) > 0
        cpg_ok = float(r.get("n_with_cpg_matrix_id_match", 0) or 0) > 0
        status = "ready" if (geno_ok and cpg_ok) else "blocked"
        ea_cis = (
            PROJECT / "meqtl-validation" / "01_cpg_meqtl_mapping" / region / "_m"
            / "tensorqtl" / "EA" / f"cpg_meqtl_{region}_EA.cis_qtl.txt.gz"
        )
        ea_burden = (
            PROJECT / "meqtl-validation" / "02_vmr_meqtl_burden" / "_m" / "EA" / region
            / "vmr_meqtl_burden.tsv.gz"
        )
        if ea_cis.exists() and ea_burden.exists():
            status = "mapped_phase1"
            reason = "EA stratified CpG cis-meQTL + Phase 2 burden complete"
        elif status == "ready":
            reason = (
                "Paths ready: all_individuals genotypes + "
                "vmr-analysis/all_individuals CpG matrices; "
                "EA stratified meQTL not yet executed"
            )
        else:
            reason = (
                "EA genotype and/or CpG matrix BrNum ID matching is zero; "
                "check all_individuals genotype + CpG paths in paths.yml"
            )
        rows.append({
            "region": region,
            "donor_group": "EA",
            "status": status,
            "n_phenotype_samples": r.get("n_phenotype_samples", ""),
            "n_with_genotype_id_match": r.get("n_with_genotype_id_match", ""),
            "n_with_cpg_matrix_id_match": r.get("n_with_cpg_matrix_id_match", ""),
            "has_predictability_summary": r.get("has_predictability_summary", ""),
            "reason": reason,
        })    # AA ready reference
    for region in REGIONS:
        sub = ov[(ov["race"] == "AA") & (ov["region"] == region)]
        if sub.empty:
            continue
        r = sub.iloc[0]
        rows.append({
            "region": region,
            "donor_group": "AA",
            "status": "mapped_phase1",
            "n_phenotype_samples": r.get("n_phenotype_samples", ""),
            "n_with_genotype_id_match": r.get("n_with_genotype_id_match", ""),
            "n_with_cpg_matrix_id_match": r.get("n_with_cpg_matrix_id_match", ""),
            "has_predictability_summary": r.get("has_predictability_summary", ""),
            "reason": "Primary Phase 1 cis-meQTL complete",
        })
    return rows


def predictability_portability(region: str, q_high: float) -> tuple[dict, pd.DataFrame]:
    aa = load_predictability(region, "AA")
    ea = load_predictability(region, "EA")
    m = aa.merge(ea, on="coord_id", suffixes=("_AA", "_EA"), how="inner")
    m = m.dropna(subset=["h2_unscaled_AA", "h2_unscaled_EA"])
    if len(m) < 20:
        return {
            "region": region,
            "n_shared_vmrs": len(m),
            "error": "insufficient overlapping VMRs",
        }, m

    r_p, p_p = pearsonr(m["h2_unscaled_AA"], m["h2_unscaled_EA"])
    r_s, p_s = spearmanr(m["h2_unscaled_AA"], m["h2_unscaled_EA"])
    thr_aa = m["h2_unscaled_AA"].quantile(q_high)
    thr_ea = m["h2_unscaled_EA"].quantile(q_high)
    m["high_AA"] = m["h2_unscaled_AA"] >= thr_aa
    m["high_EA"] = m["h2_unscaled_EA"] >= thr_ea
    m["concordant_high"] = m["high_AA"] & m["high_EA"]
    m["concordant_low"] = (~m["high_AA"]) & (~m["high_EA"])
    both_high = int(m["concordant_high"].sum())
    either_high = int((m["high_AA"] | m["high_EA"]).sum())
    summary = {
        "region": region,
        "n_shared_vmrs": int(len(m)),
        "pearson_r": float(r_p),
        "pearson_p": float(p_p),
        "spearman_r": float(r_s),
        "spearman_p": float(p_s),
        "high_pred_quantile": q_high,
        "n_concordant_high": both_high,
        "n_discordant_high": int(((m["high_AA"] ^ m["high_EA"])).sum()),
        "jaccard_high_pred": float(both_high / either_high) if either_high else np.nan,
        "median_h2_AA": float(m["h2_unscaled_AA"].median()),
        "median_h2_EA": float(m["h2_unscaled_EA"].median()),
    }
    return summary, m


def annotate_concordant(region: str, m: pd.DataFrame) -> pd.DataFrame:
    """Join Phase 2 AA burden annotations for concordant-high VMRs."""
    burden_path = PHASE2 / region / "vmr_meqtl_burden.tsv.gz"
    if not burden_path.exists() or m.empty:
        return m
    b = pd.read_csv(burden_path, sep="\t")
    b["coord_id"] = b["vmr_id"].astype(str)
    if "meqtl_supported" not in b.columns and "n_cpgs_with_sig_meqtl" in b.columns:
        b["meqtl_supported"] = pd.to_numeric(b["n_cpgs_with_sig_meqtl"], errors="coerce").fillna(0).gt(0)
    keep = [
        c for c in [
            "coord_id", "meqtl_supported", "n_cpgs_with_sig_meqtl",
            "proportion_cpgs_with_sig_meqtl", "annot.type", "line_l1_frac",
            "umap_k24_mean", "length", "h2_category",
        ] if c in b.columns
    ]
    out = m.merge(b[list(dict.fromkeys(keep))], on="coord_id", how="left")
    return out


def genomic_enrichment_concordant(region: str, m: pd.DataFrame) -> list[dict]:
    rows = []
    if m.empty or "concordant_high" not in m.columns:
        return rows
    if "annot.type" not in m.columns:
        return [{
            "region": region,
            "comparison": "concordant_high_vs_other",
            "note": "annot.type unavailable on joined table",
        }]
    use = m.dropna(subset=["annot.type"]).copy()
    hi = use["concordant_high"].astype(bool)
    for ann, grp in use.groupby("annot.type"):
        in_hi = int(grp["concordant_high"].sum())
        in_lo = int((~grp["concordant_high"]).sum())
        rows.append({
            "region": region,
            "annot_type": ann,
            "n_concordant_high": in_hi,
            "n_other": in_lo,
            "frac_of_concordant_high": in_hi / hi.sum() if hi.sum() else np.nan,
            "frac_of_other": in_lo / (~hi).sum() if (~hi).sum() else np.nan,
        })
    # meQTL support among concordant-high vs other (AA meQTL only)
    if "meqtl_supported" in use.columns:
        from scipy.stats import fisher_exact

        x = use["concordant_high"].astype(int)
        y = use["meqtl_supported"].fillna(False).astype(bool).astype(int)
        tab = pd.crosstab(x, y).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        or_, p = fisher_exact(tab.to_numpy(), alternative="greater")
        rows.append({
            "region": region,
            "annot_type": "AA_meqtl_support_enrichment_in_concordant_high",
            "n_concordant_high": int((x.eq(1) & y.eq(1)).sum()),
            "n_other": int((x.eq(0) & y.eq(1)).sum()),
            "odds_ratio": float(or_),
            "pvalue": float(p),
            "note": "Uses AA Phase 1 meQTL only; not EA meQTL",
        })
    return rows


def aa_burden_coefficients() -> list[dict]:
    rows = []
    for region in REGIONS:
        p = PHASE2 / region / "burden_model_results.tsv"
        if not p.exists():
            continue
        df = pd.read_csv(p, sep="\t")
        for _, r in df.iterrows():
            rows.append({
                "region": region,
                "donor_group": "AA",
                "model": r.get("model", ""),
                "coef_predictability": r.get("coef_predictability", ""),
                "se_predictability": r.get("se_predictability", ""),
                "pval_predictability": r.get("pval_predictability", ""),
                "n_vmrs": r.get("n_vmrs", ""),
            })
        ea_p = PHASE2 / "EA" / region / "burden_model_results.tsv"
        if ea_p.exists():
            edf = pd.read_csv(ea_p, sep="\t")
            for _, r in edf.iterrows():
                rows.append({
                    "region": region,
                    "donor_group": "EA",
                    "model": r.get("model", ""),
                    "coef_predictability": r.get("coef_predictability", ""),
                    "se_predictability": r.get("se_predictability", ""),
                    "pval_predictability": r.get("pval_predictability", ""),
                    "n_vmrs": r.get("n_vmrs", ""),
                    "notes": "EA stratified meQTL burden (M0 covariates)",
                })
        else:
            rows.append({
                "region": region,
                "donor_group": "EA",
                "model": "not_run",
                "coef_predictability": "",
                "se_predictability": "",
                "pval_predictability": "",
                "n_vmrs": "",
                "notes": "EA stratified meQTL + burden not yet available",
            })
    return rows


def main() -> None:
    args = parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    regions = REGIONS if args.region == "all" else [args.region]

    status = ea_meqtl_status()
    write_tsv(OUTDIR / "ea_meqtl_readiness.tsv", status)

    port_rows = []
    annot_rows = []
    for region in regions:
        summary, m = predictability_portability(region, args.high_pred_quantile)
        port_rows.append(summary)
        if m.empty:
            continue
        m = annotate_concordant(region, m)
        reg_out = OUTDIR / region
        reg_out.mkdir(parents=True, exist_ok=True)
        m.to_csv(reg_out / "aa_ea_predictability_joined.tsv.gz", sep="\t", index=False, compression="gzip")
        if "concordant_high" in m.columns:
            conc = m.loc[m["concordant_high"].astype(bool)]
            conc.to_csv(
                reg_out / "concordant_high_predictability_vmrs.tsv.gz",
                sep="\t", index=False, compression="gzip",
            )
        annot_rows.extend(genomic_enrichment_concordant(region, m))
        write_tsv(reg_out / "donor_group_coefficient_comparison.tsv", [
            r for r in aa_burden_coefficients() if r["region"] == region
        ])

    write_tsv(OUTDIR / "aa_ea_predictability_portability.tsv", port_rows)
    if annot_rows:
        write_tsv(OUTDIR / "concordant_high_genomic_annotation.tsv", annot_rows)
    write_tsv(OUTDIR / "donor_group_coefficient_comparison.tsv", aa_burden_coefficients())

    print(f"Wrote donor-group outputs under {OUTDIR}")
    print(pd.DataFrame(port_rows).to_string(index=False))
    print(pd.DataFrame(status)[["region", "donor_group", "status"]].to_string(index=False))


if __name__ == "__main__":
    main()
