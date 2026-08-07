#!/usr/bin/env python3
"""Run Phase 3 enrichment for primary external resources × regions and summarize."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
H = PROJECT / "meqtl-validation" / "03_external_meqtl_validation" / "_h"
OUT = PROJECT / "meqtl-validation" / "03_external_meqtl_validation" / "_m"

RESOURCES = [
    "jaffe_dlpfc_450k_meqtl",
    "schulz_hippocampus_array_meqtl",
    "brainseq_wgbs_meqtl_scz_subset",  # exploratory only
]
REGIONS = ["caudate", "dlpfc", "hippocampus"]


def main() -> None:
    script = H / "01_test_external_support.py"
    for resource in RESOURCES:
        for region in REGIONS:
            harm = OUT / "harmonized" / f"{resource}.{region}.vmr_support.tsv.gz"
            if not harm.exists():
                print(f"SKIP missing {harm}")
                continue
            cmd = [
                sys.executable, str(script),
                "--region", region,
                "--resource-id", resource,
            ]
            print("RUN", " ".join(cmd))
            subprocess.check_call(cmd)

    model_rows = []
    match_rows = []
    for region in REGIONS:
        for resource in RESOURCES:
            mp = OUT / region / f"external_support_model_{resource}.tsv"
            kp = OUT / region / f"external_matched_{resource}.tsv"
            if mp.exists():
                model_rows.append(pd.read_csv(mp, sep="\t"))
            if kp.exists():
                match_rows.append(pd.read_csv(kp, sep="\t"))

    if model_rows:
        models = pd.concat(model_rows, ignore_index=True)
        models.to_csv(OUT / "external_support_models_all.tsv", sep="\t", index=False)
        # Primary decision rows: preferred tissue pairings, adjusted_minimal if present
        primary = models[
            (models["analysis_role"] == "primary")
            & (models["model"].isin(["unadjusted", "adjusted_minimal", "adjusted_technical"]))
        ].copy()
        primary.to_csv(OUT / "external_support_primary_summary.tsv", sep="\t", index=False)
        print(f"Wrote {OUT / 'external_support_models_all.tsv'}")
        print(f"Wrote {OUT / 'external_support_primary_summary.tsv'}")

    if match_rows:
        matched = pd.concat(match_rows, ignore_index=True)
        matched.to_csv(OUT / "external_matched_all.tsv", sep="\t", index=False)
        print(f"Wrote {OUT / 'external_matched_all.tsv'}")

    # Compact pass/fail for §7.5 criterion 5
    verdict = []
    if model_rows:
        m = pd.concat(model_rows, ignore_index=True)
        focus = m[
            (m["analysis_role"] == "primary")
            & (m["model"] == "adjusted_minimal")
            & (~m["resource_id"].isin(["brainseq_wgbs_meqtl_scz_subset"]))
        ]
        for _, row in focus.iterrows():
            coef = row.get("coef_predictability")
            pval = row.get("pval_predictability")
            ok = pd.notna(coef) and coef > 0 and pd.notna(pval) and pval < 0.05
            verdict.append({
                "resource_id": row["resource_id"],
                "region": row["region"],
                "model": row["model"],
                "coef_predictability": coef,
                "pval_predictability": pval,
                "supports_gradient": bool(ok),
            })
    write_tsv(OUT / "phase3_criterion5_verdict.tsv", verdict)
    print(f"Wrote {OUT / 'phase3_criterion5_verdict.tsv'}")
    if verdict:
        n_pass = sum(v["supports_gradient"] for v in verdict)
        print(f"Criterion 5: {n_pass}/{len(verdict)} primary resource×region tests support gradient")


if __name__ == "__main__":
    main()
