#!/usr/bin/env python3
"""Summarize Phase 5 expression/PSI enrichment across regions and modalities."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PHASE5 = Path(
    "/projects/b1213/users/kynon/projects/dna-methylation-heritability/"
    "meqtl-validation/06_transcription_splicing_integration/_m"
)
REGIONS = ["caudate", "dlpfc", "hippocampus"]
MODALITIES = ["expression", "psi", "expression_abc"]


def main() -> None:
    rows = []
    for region in REGIONS:
        for mod in MODALITIES:
            p = PHASE5 / region / f"meqtl_x_{mod}_enrichment.tsv"
            if p.exists():
                rows.append(pd.read_csv(p, sep="\t"))
    if not rows:
        raise SystemExit("No Phase 5 enrichment tables found")
    all_df = pd.concat(rows, ignore_index=True)
    all_df.to_csv(PHASE5 / "tx_enrichment_all.tsv", sep="\t", index=False)

    # Primary claim rows: fisher + adjusted logistic meQTL support
    primary = all_df[
        all_df["model"].isin([
            "fisher_meqtl_support_x_tx",
            "logistic_tx_~_meqtl_support_adjusted",
            "logistic_tx_~_predictability_adjusted",
        ])
        & (all_df["modality"].isin(["expression", "psi"]))
    ].copy()
    primary.to_csv(PHASE5 / "tx_enrichment_primary.tsv", sep="\t", index=False)

    # Both expression and PSI: VMR intersection enrichment if joined tables exist
    both_rows = []
    for region in REGIONS:
        e = PHASE5 / region / "vmr_meqtl_expression_joined.tsv.gz"
        s = PHASE5 / region / "vmr_meqtl_psi_joined.tsv.gz"
        if not (e.exists() and s.exists()):
            continue
        de = pd.read_csv(e, sep="\t")
        ds = pd.read_csv(s, sep="\t")
        m = de[["coord_id", "meqtl_supported", "tx_associated", "predictability"]].merge(
            ds[["coord_id", "tx_associated"]].rename(columns={"tx_associated": "psi_associated"}),
            on="coord_id",
            how="inner",
        )
        m["both_tx"] = m["tx_associated"].astype(bool) & m["psi_associated"].astype(bool)
        m["meqtl_supported"] = m["meqtl_supported"].astype(bool)
        from scipy.stats import fisher_exact

        x = m["meqtl_supported"].astype(int)
        y = m["both_tx"].astype(int)
        tab = pd.crosstab(x, y).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        or_, p = fisher_exact(tab.to_numpy(), alternative="greater")
        both_rows.append({
            "region": region,
            "modality": "expression_and_psi",
            "model": "fisher_meqtl_support_x_both_tx",
            "n_vmrs": len(m),
            "n_meqtl_supported": int(m["meqtl_supported"].sum()),
            "n_tx_associated": int(m["both_tx"].sum()),
            "n_both": int((m["meqtl_supported"] & m["both_tx"]).sum()),
            "estimate": float(or_),
            "pvalue": float(p),
        })
    if both_rows:
        pd.DataFrame(both_rows).to_csv(PHASE5 / "tx_enrichment_both_modalities.tsv", sep="\t", index=False)

    # Claim summary
    claim = []
    for mod in ["expression", "psi"]:
        sub = all_df[(all_df["modality"] == mod) & (all_df["model"] == "fisher_meqtl_support_x_tx")]
        n_pos = int(((sub["estimate"] > 1) & (sub["pvalue"] < 0.05)).sum())
        claim.append({
            "modality": mod,
            "n_regions_significant_positive": n_pos,
            "n_regions": len(sub),
            "passes_phase5": n_pos >= 2 or (n_pos >= 1 and (sub["estimate"] > 1).sum() >= 2),
        })
    write_tsv(PHASE5 / "phase5_claim_summary.tsv", claim)
    print(f"Wrote summaries under {PHASE5}")
    print(primary[["region", "modality", "model", "n_both", "estimate", "pvalue"]].to_string(index=False))
    print(pd.DataFrame(claim).to_string(index=False))


if __name__ == "__main__":
    main()
