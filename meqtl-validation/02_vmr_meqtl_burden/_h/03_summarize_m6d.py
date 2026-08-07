#!/usr/bin/env python3
"""Apply prespecified manuscript-retention criteria to M6d versus locked M3a."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
REGIONS = ("caudate", "dlpfc", "hippocampus")


def burden_effect(path: Path) -> tuple[float, float]:
    if not path.exists():
        return np.nan, np.nan
    tab = pd.read_csv(path, sep="\t")
    row = tab.loc[tab["model"].eq("adjusted_technical")]
    if row.empty:
        row = tab.loc[tab["model"].eq("adjusted_minimal")]
    if row.empty:
        return np.nan, np.nan
    return float(row.iloc[0].get("coef_predictability", np.nan)), float(row.iloc[0].get("pval_predictability", np.nan))


rows = []
for region in REGIONS:
    comparison = (
        ROOT / "meqtl-validation/01_cpg_meqtl_mapping" / region / "_m"
        / "covariate_sensitivity/m6d_vs_m3a_matched_summary.tsv"
    )
    row = {"region": region, "status": "missing"}
    if comparison.exists():
        comp = pd.read_csv(comparison, sep="\t")
        base = comp.loc[comp["model_id"].eq("M3a_dnam_matched")]
        m6d = comp.loc[comp["model_id"].eq("M6d")]
        if len(base) == 1 and len(m6d) == 1 and m6d.iloc[0].get("status") == "ok":
            base_n = float(base.iloc[0]["n_significant_fdr"])
            m6d_n = float(m6d.iloc[0]["n_significant_fdr"])
            row.update({
                "status": "ok",
                "m3a_discoveries": int(base_n),
                "m6d_discoveries": int(m6d_n),
                "discovery_retention": m6d_n / base_n if base_n > 0 else np.nan,
                "external_enrichment_retention": float(m6d.iloc[0].get("enrichment_or_vs_baseline", np.nan)),
            })
    coef, pval = burden_effect(
        ROOT / "meqtl-validation/02_vmr_meqtl_burden/_m/M6d" / region / "burden_model_results.tsv"
    )
    row["m6d_predictability_coef"] = coef
    row["m6d_predictability_pvalue"] = pval
    row["positive_burden_effect"] = bool(np.isfinite(coef) and coef > 0)
    row["discovery_retention_pass"] = bool(np.isfinite(row.get("discovery_retention", np.nan)) and row["discovery_retention"] >= 0.50)
    row["external_enrichment_retention_pass"] = bool(np.isfinite(row.get("external_enrichment_retention", np.nan)) and row["external_enrichment_retention"] >= 0.80)
    rows.append(row)

out = pd.DataFrame(rows)
n_positive = int(out["positive_burden_effect"].sum())
out["overall_retention_pass"] = (
    (n_positive >= 2)
    & out["discovery_retention_pass"].fillna(False)
    & out["external_enrichment_retention_pass"].fillna(False)
)
dest = ROOT / "meqtl-validation/02_vmr_meqtl_burden/_m/M6d/m6d_robustness_summary.tsv"
dest.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(dest, sep="\t", index=False)
print(out.to_string(index=False))
print(f"Wrote {dest}")
