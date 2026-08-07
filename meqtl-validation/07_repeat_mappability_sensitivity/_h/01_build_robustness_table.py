#!/usr/bin/env python3
"""Build consolidated repeat/mappability robustness table scaffold.

Produces one table with columns for original / adjusted / matched /
high-mappability / SNP-proximity-excluded / segdup-excluded estimates.
Until Phase 1–2 outputs exist, writes a template with required rows.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--annotation-enrichment-tsv",
        default="",
        help="Optional existing LINE/H3K9me3 enrichment results to seed 'original' column",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(
        "/projects/b1213/users/kynon/projects/dna-methylation-heritability/"
        "meqtl-validation/07_repeat_mappability_sensitivity/_m"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    analyses = [
        "LINE_L1_enrichment_high_vs_low_predictability",
        "H3K9me3_enrichment_high_vs_low_predictability",
        "quiescent_chromatin_enrichment_high_vs_low_predictability",
        "predictability_meqtl_burden_association",
    ]
    rows = []
    for a in analyses:
        rows.append({
            "analysis": a,
            "original_estimate": "",
            "original_pvalue": "",
            "adjusted_estimate": "",
            "adjusted_pvalue": "",
            "matched_estimate": "",
            "matched_pvalue": "",
            "high_mappability_estimate": "",
            "high_mappability_pvalue": "",
            "snp_proximity_excluded_estimate": "",
            "snp_proximity_excluded_pvalue": "",
            "segdup_excluded_estimate": "",
            "segdup_excluded_pvalue": "",
            "direction_consistent": "",
            "significance_consistent": "",
            "status": "awaiting_phase1_2_inputs",
        })

    if args.annotation_enrichment_tsv and Path(args.annotation_enrichment_tsv).exists():
        src = pd.read_csv(args.annotation_enrichment_tsv, sep="\t")
        # Best-effort: if columns look usable, note path for manual fill
        rows.append({
            "analysis": "SOURCE_FILE_NOTED",
            "original_estimate": str(args.annotation_enrichment_tsv),
            "status": "seed_path_recorded",
            "direction_consistent": "",
            "significance_consistent": "",
            "original_pvalue": "",
            "adjusted_estimate": "",
            "adjusted_pvalue": "",
            "matched_estimate": "",
            "matched_pvalue": "",
            "high_mappability_estimate": "",
            "high_mappability_pvalue": "",
            "snp_proximity_excluded_estimate": "",
            "snp_proximity_excluded_pvalue": "",
            "segdup_excluded_estimate": "",
            "segdup_excluded_pvalue": "",
        })
        _ = src  # reserved for future auto-fill

    write_tsv(outdir / "consolidated_robustness_table.tsv", rows)
    write_tsv(
        outdir / "sensitivity_plan.tsv",
        [
            {"step": "restrict_high_mappability", "threshold": ">=0.9", "status": "prespecified"},
            {"step": "exclude_segmental_duplications", "threshold": "as annotated", "status": "prespecified"},
            {"step": "exclude_snp_proximal_cpgs", "threshold": "150bp", "status": "prespecified"},
            {"step": "match_coverage_variance_gc_cpg_density_length_snp_density_annotation", "threshold": "SMD<=0.1", "status": "prespecified"},
            {"step": "repeat_LINE_L1_and_H3K9me3_after_restrictions", "threshold": "n/a", "status": "prespecified"},
        ],
    )
    print(f"Wrote robustness scaffold under {outdir}")


if __name__ == "__main__":
    main()
