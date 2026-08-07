#!/usr/bin/env python3
"""Summarize Phase 4 cross-region + donor-group claim status."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

CROSS = Path(
    "/projects/b1213/users/kynon/projects/dna-methylation-heritability/"
    "meqtl-validation/04_cross_region_sharing/_m"
)
DONOR = Path(
    "/projects/b1213/users/kynon/projects/dna-methylation-heritability/"
    "meqtl-validation/05_donor_group_comparison/_m"
)


def main() -> None:
    claims = []

    cpg = pd.read_csv(CROSS / "cross_region_cpg_concordance.tsv", sep="\t")
    vmr = pd.read_csv(CROSS / "cross_region_vmr_sharing.tsv", sep="\t")
    pair_cpg = cpg[cpg["contrast"].str.contains("_vs_")]
    dir_ok = (pair_cpg["direction_concordance_both_sig"] > 0.7).sum()
    zcorr_ok = ((pair_cpg["pearson_z_either_sig"] > 0) & (pair_cpg["pearson_z_either_sig_p"] < 0.05)).sum()
    claims.append({
        "claim": "cpg_effect_direction_concordance",
        "metric": "direction_concordance_both_sig>0.7",
        "n_pairwise_pass": int(dir_ok),
        "n_pairwise": int(len(pair_cpg)),
        "passes": bool(dir_ok >= 2),
        "detail": "; ".join(
            f"{r.contrast}:{r.direction_concordance_both_sig:.3f}"
            for r in pair_cpg.itertuples()
        ),
    })
    claims.append({
        "claim": "cpg_standardized_effect_correlation",
        "metric": "pearson_z_either_sig>0 & p<0.05",
        "n_pairwise_pass": int(zcorr_ok),
        "n_pairwise": int(len(pair_cpg)),
        "passes": bool(zcorr_ok >= 2),
        "detail": "; ".join(
            f"{r.contrast}:r={r.pearson_z_either_sig:.3f}"
            for r in pair_cpg.itertuples()
        ),
    })

    pair_vmr = vmr[vmr["contrast"].str.contains("_vs_")].copy()
    pred_ok = (
        (pair_vmr.get("logistic_or_pred_shared", pd.Series(dtype=float)) > 1)
        & (pair_vmr.get("logistic_p_pred_shared", pd.Series(dtype=float)) < 0.05)
    ).sum() if "logistic_p_pred_shared" in pair_vmr.columns else 0
    claims.append({
        "claim": "predictability_associates_with_shared_meqtl_support",
        "metric": "logistic OR>1 & p<0.05",
        "n_pairwise_pass": int(pred_ok),
        "n_pairwise": int(len(pair_vmr)),
        "passes": bool(pred_ok >= 2),
        "detail": "; ".join(
            f"{r.contrast}:OR={getattr(r,'logistic_or_pred_shared', float('nan')):.3f}"
            for r in pair_vmr.itertuples()
        ) if "logistic_or_pred_shared" in pair_vmr.columns else "",
    })

    grad = pd.read_csv(CROSS / "burden_gradient_by_region.tsv", sep="\t")
    tech = grad[grad["model"] == "adjusted_technical"]
    grad_ok = ((tech["coef_predictability"].astype(float) > 0) & (tech["pval_predictability"].astype(float) < 0.05)).sum()
    claims.append({
        "claim": "aa_burden_gradient_reproducible_across_regions",
        "metric": "adjusted_technical coef>0 & p<0.05",
        "n_pairwise_pass": int(grad_ok),
        "n_pairwise": int(len(tech)),
        "passes": bool(grad_ok >= 2),
        "detail": "; ".join(f"{r.region}:{float(r.coef_predictability):.3f}" for r in tech.itertuples()),
    })

    port = pd.read_csv(DONOR / "aa_ea_predictability_portability.tsv", sep="\t")
    port_ok = ((port["spearman_r"] > 0) & (port["spearman_p"] < 0.05)).sum()
    claims.append({
        "claim": "aa_ea_predictability_portability",
        "metric": "spearman_r>0 & p<0.05",
        "n_pairwise_pass": int(port_ok),
        "n_pairwise": int(len(port)),
        "passes": bool(port_ok >= 2),
        "detail": "; ".join(f"{r.region}:rho={r.spearman_r:.3f}" for r in port.itertuples()),
    })

    ready = pd.read_csv(DONOR / "ea_meqtl_readiness.tsv", sep="\t")
    ea_rows = ready[ready["donor_group"] == "EA"]
    ea_mapped = ea_rows[ea_rows["status"] == "mapped_phase1"]
    coef = pd.read_csv(DONOR / "donor_group_coefficient_comparison.tsv", sep="\t")
    ea_tech = coef[(coef["donor_group"] == "EA") & (coef["model"] == "adjusted_technical")].copy()
    if ea_tech.empty:
        ea_tech = coef[(coef["donor_group"] == "EA") & (coef["model"] == "adjusted_minimal")].copy()
    if ea_tech.empty:
        ea_tech = coef[(coef["donor_group"] == "EA") & (coef["model"] == "unadjusted")].copy()
    ea_grad_ok = 0
    detail_bits = []
    if not ea_tech.empty:
        ea_tech["coef_predictability"] = pd.to_numeric(ea_tech["coef_predictability"], errors="coerce")
        ea_tech["pval_predictability"] = pd.to_numeric(ea_tech["pval_predictability"], errors="coerce")
        ea_grad_ok = int(
            ((ea_tech["coef_predictability"] > 0) & (ea_tech["pval_predictability"] < 0.05)).sum()
        )
        detail_bits = [
            f"{r.region}:coef={float(r.coef_predictability):.3f}"
            for r in ea_tech.itertuples()
            if pd.notna(r.coef_predictability)
        ]
    passes_ea = bool(len(ea_mapped) >= 2 and ea_grad_ok >= 2)
    claims.append({
        "claim": "ea_stratified_meqtl_burden",
        "metric": "EA meQTL mapped + positive burden~predictability in ≥2 regions",
        "n_pairwise_pass": int(ea_grad_ok),
        "n_pairwise": int(len(ea_rows)),
        "passes": passes_ea,
        "detail": (
            f"mapped={len(ea_mapped)}/{len(ea_rows)}; "
            + ("; ".join(detail_bits) if detail_bits else "EA burden not available")
        ),
    })

    # Experiment 3: caudate downsample + G×region
    pending = CROSS / "pending_analyses.tsv"
    if pending.exists():
        pend = pd.read_csv(pending, sep="\t")
        for analysis, claim_name in [
            ("caudate_donor_downsample_remap", "caudate_n_matched_downsample"),
            ("shared_donor_genotype_by_region", "shared_donor_gxregion"),
        ]:
            sub = pend[pend["analysis"] == analysis]
            done = bool(len(sub) and str(sub.iloc[0]["status"]) == "done")
            detail = str(sub.iloc[0]["reason"]) if len(sub) else "missing"
            # Enrich with claim snapshots when present
            if analysis == "caudate_donor_downsample_remap":
                snap = CROSS / "caudate_downsample/downsample_claim_snapshot.tsv"
                if snap.exists():
                    s = pd.read_csv(snap, sep="\t").iloc[0]
                    detail = (
                        f"not_solely_N={s.get('criterion_not_solely_sample_size')}; "
                        f"median_n_sig={s.get('downsample_median_n_sig')}; "
                        f"ratio_dlpfc={s.get('median_ratio_vs_dlpfc')}"
                    )
            if analysis == "shared_donor_genotype_by_region":
                snap = CROSS / "gxregion/gxregion_claim_snapshot.tsv"
                if snap.exists():
                    s = pd.read_csv(snap, sep="\t").iloc[0]
                    detail = (
                        f"fitted={s.get('n_pairs_fitted')}; "
                        f"joint_FDR_sig={s.get('n_sig_joint_interaction_fdr')}; "
                        f"frac={s.get('frac_sig_joint_interaction')}"
                    )
            claims.append({
                "claim": claim_name,
                "metric": "status=done in pending_analyses.tsv (+ claim snapshot)",
                "n_pairwise_pass": int(done),
                "n_pairwise": 1,
                "passes": done,
                "detail": detail,
            })

    # Experiment 2 depth: AA–EA concordance + MAF/LD matching
    conc = DONOR / "aa_ea_effect_concordance_summary.tsv"
    if conc.exists():
        c = pd.read_csv(conc, sep="\t")
        dir_ok = (
            (c["direction_concordance_both_sig"] > 0.7)
            & (c["n_sig_both"] > 0)
        ).sum() if "direction_concordance_both_sig" in c.columns else 0
        claims.append({
            "claim": "aa_ea_cpg_effect_concordance",
            "metric": "direction_concordance_both_sig>0.7 in ≥2 regions",
            "n_pairwise_pass": int(dir_ok),
            "n_pairwise": int(len(c)),
            "passes": bool(dir_ok >= 2),
            "detail": "; ".join(
                f"{r.region}:dir={r.direction_concordance_both_sig:.3f},n_both={int(r.n_sig_both)}"
                for r in c.itertuples()
                if hasattr(r, "direction_concordance_both_sig")
            ),
        })
    maf = DONOR / "experiment2_depth_claim_summary.tsv"
    if maf.exists():
        m = pd.read_csv(maf, sep="\t").iloc[0]
        claims.append({
            "claim": "aa_ea_maf_ld_matched_discovery",
            "metric": "documented MAF/cis-SNP-density sensitivity (no ancestry claim)",
            "n_pairwise_pass": int(m.get("n_regions_gap_shrinks", 0)),
            "n_pairwise": int(m.get("n_regions", 0)),
            "passes": bool(m.get("passes_document_maf_ld_sensitivity", False)),
            "detail": str(m.get("detail", "")),
        })

    write_tsv(CROSS / "phase4_claim_summary.tsv", claims)
    # also copy under donor folder for discoverability
    write_tsv(DONOR / "phase4_claim_summary.tsv", claims)
    print(pd.DataFrame(claims)[["claim", "passes", "n_pairwise_pass", "n_pairwise", "detail"]].to_string(index=False))


if __name__ == "__main__":
    main()
