#!/usr/bin/env python3
"""Summarize Level-3 shared genetic support across LIBD eQTL + GTEx."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
LEVEL3 = ROOT / "meqtl-validation/08_schizophrenia_risk_application/_m/level3"
PRIOR = (
    ROOT
    / "meqtl-validation/08_schizophrenia_risk_application/_m/prioritized/prioritized_loci.tsv"
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default=str(LEVEL3))
    args = ap.parse_args()
    outdir = Path(args.outdir)

    prior = pd.read_csv(PRIOR, sep="\t")
    libd_path = outdir / "libd_risk_variant_eqtl.tsv.gz"
    if not libd_path.exists():
        libd_path = outdir / "libd_eqtl/libd_risk_variant_eqtl.tsv.gz"
    gtex_path = outdir / "gtex/gtex_level3_locus_gene_summary.tsv"
    libd = pd.read_csv(libd_path, sep="\t") if libd_path.exists() else pd.DataFrame()
    gtex = pd.read_csv(gtex_path, sep="\t") if gtex_path.exists() else pd.DataFrame()

    rows = []
    for _, r in prior.iterrows():
        snp = r["index_snp"]
        loc = r["locus_id"]
        L = libd[libd["index_snp"] == snp] if len(libd) else pd.DataFrame()
        G = gtex[gtex["index_snp"] == snp] if len(gtex) else pd.DataFrame()
        libd_sig = L[L.get("significant_fdr", False) == True] if len(L) else L  # noqa: E712
        best = (
            libd_sig.sort_values("eqtl_pval").iloc[0]
            if len(libd_sig)
            else (L.sort_values("eqtl_pval").iloc[0] if len(L) and "eqtl_pval" in L else None)
        )
        rows.append(
            {
                "locus_id": loc,
                "index_snp": snp,
                "finemap_genes": r.get("finemap_genes"),
                "tx_expression": r.get("tx_expression"),
                "tx_psi": r.get("tx_psi"),
                "libd_n_tested": int((L["status"] == "tested").sum()) if len(L) else 0,
                "libd_n_sig_fdr": int(libd_sig.shape[0]) if len(L) else 0,
                "libd_best_gene": best["gene_symbol"] if best is not None else "",
                "libd_best_pval": best["eqtl_pval"] if best is not None else pd.NA,
                "libd_best_qval": best["qval"] if best is not None else pd.NA,
                "libd_best_beta_risk": best["eqtl_beta_risk"]
                if best is not None and "eqtl_beta_risk" in best
                else pd.NA,
                "libd_same_sign_meqtl": bool(best["same_sign_meqtl_eqtl_risk"])
                if best is not None and "same_sign_meqtl_eqtl_risk" in best
                else False,
                "gtex_caudate_gene_variant": bool(G["gtex_caudate_eqtl_gene_variant"].any())
                if len(G)
                else False,
                "gtex_any_eqtl_gene": bool(G["gtex_any_eqtl_gene"].any()) if len(G) else False,
                "level3_pass": bool(
                    (len(libd_sig) > 0)
                    or (len(G) and G["gtex_caudate_eqtl_gene_variant"].any())
                ),
            }
        )

    locus = pd.DataFrame(rows)
    locus.to_csv(outdir / "level3_locus_summary.tsv", sep="\t", index=False)

    claim = pd.DataFrame(
        [
            {
                "n_prioritized_loci": int(len(prior)),
                "n_loci_libd_sig": int((locus["libd_n_sig_fdr"] > 0).sum()),
                "n_loci_gtex_caudate_gene_variant": int(
                    locus["gtex_caudate_gene_variant"].sum()
                ),
                "n_loci_level3_pass": int(locus["level3_pass"].sum()),
                "success_criterion": ">=1 locus with risk->meQTL and risk->eQTL/sQTL",
                "pass": bool(locus["level3_pass"].any()),
                "libd_recipe": "AA_AgeGt13_ControlSCZD; Sex+Dx+Age+snpPC1-3+exprPCs; filterByExpr+TMM-logCPM",
            }
        ]
    )
    claim.to_csv(outdir / "level3_claim_snapshot.tsv", sep="\t", index=False)
    print(claim.to_string(index=False))
    print(locus[["index_snp", "libd_n_sig_fdr", "gtex_caudate_gene_variant", "level3_pass"]].to_string(index=False))


if __name__ == "__main__":
    main()
