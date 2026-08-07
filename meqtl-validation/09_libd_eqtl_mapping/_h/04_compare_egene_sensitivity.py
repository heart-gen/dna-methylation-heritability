#!/usr/bin/env python3
"""Compare eGene discovery: CPM primary vs RPKM sensitivity."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
M = ROOT / "meqtl-validation/09_libd_eqtl_mapping/_m/caudate"
OUT = M / "sensitivity_rpkm_vs_cpm.tsv"


def summarize(label: str, cis: Path, prep: Path, cov: Path) -> dict:
    df = pd.read_csv(cis, sep="\t")
    prep_df = pd.read_csv(prep, sep="\t") if prep.exists() else pd.DataFrame()
    cov_cols = pd.read_csv(cov, sep="\t", nrows=0).columns.tolist() if cov.exists() else []
    n_pcs = sum(c.startswith("PC") and c[2:].isdigit() for c in cov_cols)
    row = {
        "label": label,
        "cis_path": str(cis),
        "n_phenotypes_mapped": int(len(df)),
        "n_egene_fdr05": int((df["qval"] <= 0.05).sum()) if "qval" in df else 0,
        "n_egene_fdr10": int((df["qval"] <= 0.10).sum()) if "qval" in df else 0,
        "n_egene_fdr20": int((df["qval"] <= 0.20).sum()) if "qval" in df else 0,
        "min_qval": float(df["qval"].min()) if "qval" in df and len(df) else float("nan"),
        "n_expr_pcs": n_pcs,
    }
    if len(prep_df):
        for k in (
            "n_samples",
            "n_features",
            "phenotype",
            "feature_filter",
            "norm_method",
        ):
            if k in prep_df.columns:
                row[k] = prep_df.iloc[0][k]
    return row


def main() -> None:
    arms = [
        (
            "cpm_filterByExpr",
            M / "genes/tensorqtl/libd_aa_caudate_genes_standard.cis_qtl.txt.gz",
            M / "genes/prepared/prep_summary.tsv",
            M / "genes/standard/covariates.txt",
        ),
        (
            "rpkm_mean0.2",
            M / "genes_rpkm/tensorqtl/libd_aa_caudate_genes_rpkm.cis_qtl.txt.gz",
            M / "genes_rpkm/prepared/prep_summary.tsv",
            M / "genes_rpkm/standard/covariates.txt",
        ),
    ]
    rows = []
    for label, cis, prep, cov in arms:
        if not cis.exists():
            rows.append({"label": label, "cis_path": str(cis), "status": "missing"})
            continue
        rows.append({**summarize(label, cis, prep, cov), "status": "ok"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
