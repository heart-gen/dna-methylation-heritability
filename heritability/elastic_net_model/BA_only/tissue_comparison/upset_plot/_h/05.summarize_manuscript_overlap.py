#!/usr/bin/env python
"""
Summarize F_0.25 overlap results for the manuscript-ready overlap figure.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    import session_info
except ImportError:  # pragma: no cover - optional at runtime
    session_info = None


RECIPROCAL_THRESHOLD = 25.0
OVERLAP_FLAG = "F_0.25"
REGIONS = ("caudate", "dlpfc", "hippocampus")

CATEGORY_ORDER = ["all", "heritable", "non-heritable", "low_prediction"]
CATEGORY_LABELS = {
    "all": "All VMRs",
    "heritable": "Heritable",
    "non-heritable": "Non-heritable",
    "low_prediction": "Low prediction",
}

PAIR_ORDER = [
    ("caudate", "dlpfc"),
    ("caudate", "hippocampus"),
    ("hippocampus", "dlpfc"),
]
PAIR_LABELS = {
    ("caudate", "dlpfc"): "Caudate-DLPFC",
    ("caudate", "hippocampus"): "Caudate-Hippocampus",
    ("hippocampus", "dlpfc"): "Hippocampus-DLPFC",
}
PAIR_LABELS_MULTILINE = {
    ("caudate", "dlpfc"): "Caudate\nDLPFC",
    ("caudate", "hippocampus"): "Caudate\nHippocampus",
    ("hippocampus", "dlpfc"): "Hippocampus\nDLPFC",
}

SET_ORDER = [
    "Caudate",
    "DLPFC",
    "Hippocampus",
    "Caudate&DLPFC",
    "Caudate&Hippocampus",
    "Hippocampus&DLPFC",
    "Caudate&Hippocampus&DLPFC",
]
SET_MEMBERSHIP = {
    "Caudate": {"includes_caudate": 1, "includes_dlpfc": 0, "includes_hippocampus": 0},
    "DLPFC": {"includes_caudate": 0, "includes_dlpfc": 1, "includes_hippocampus": 0},
    "Hippocampus": {"includes_caudate": 0, "includes_dlpfc": 0, "includes_hippocampus": 1},
    "Caudate&DLPFC": {"includes_caudate": 1, "includes_dlpfc": 1, "includes_hippocampus": 0},
    "Caudate&Hippocampus": {"includes_caudate": 1, "includes_dlpfc": 0, "includes_hippocampus": 1},
    "Hippocampus&DLPFC": {"includes_caudate": 0, "includes_dlpfc": 1, "includes_hippocampus": 1},
    "Caudate&Hippocampus&DLPFC": {"includes_caudate": 1, "includes_dlpfc": 1, "includes_hippocampus": 1},
}


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    """Adjust p-values using the Benjamini-Hochberg procedure."""
    pvals = np.asarray(p_values, dtype=float)
    adjusted = np.full(pvals.shape, np.nan)

    valid = np.isfinite(pvals)
    if not valid.any():
        return adjusted

    valid_pvals = pvals[valid]
    order = np.argsort(valid_pvals)
    ranked = valid_pvals[order]
    n_tests = len(ranked)

    correction = ranked * n_tests / np.arange(1, n_tests + 1)
    correction = np.minimum.accumulate(correction[::-1])[::-1]
    correction = np.clip(correction, 0, 1)

    reordered = np.empty_like(correction)
    reordered[order] = correction
    adjusted[valid] = reordered
    return adjusted


def interval_id(df: pd.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pd.Series:
    return (
        df[chrom_col].astype(str)
        + ":"
        + df[start_col].astype(int).astype(str)
        + "-"
        + df[end_col].astype(int).astype(str)
    )


def classify_vmrs(enet: pd.DataFrame) -> pd.DataFrame:
    vmr = enet.dropna().copy()
    vmr["h2_category"] = np.where(
        vmr["r_squared_cv"] <= 0.3,
        "low_prediction",
        np.where(vmr["h2_unscaled"] >= 0.1, "heritable", "non-heritable"),
    )
    vmr["interval_id"] = interval_id(vmr, "chrom", "start", "end")
    return vmr


def load_region_summaries(analysis_dir: Path) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    ba_only_dir = analysis_dir.parents[1]
    region_tables: dict[str, pd.DataFrame] = {}
    total_rows: list[dict[str, object]] = []

    for region in REGIONS:
        region_file = ba_only_dir / region / "_m" / f"{region}_summary_elastic-net.tsv"
        vmr = pd.read_csv(region_file, sep="\t", dtype={"chrom": str})
        vmr = classify_vmrs(vmr)

        if vmr["interval_id"].duplicated().any():
            duplicates = int(vmr["interval_id"].duplicated().sum())
            raise ValueError(f"{region_file} contains {duplicates} duplicated interval IDs.")

        region_tables[region] = vmr[
            ["interval_id", "chrom", "start", "end", "h2_unscaled", "r_squared_cv", "h2_category"]
        ].copy()

        for h2_category in CATEGORY_ORDER:
            total_rows.append(
                {
                    "region": region,
                    "h2_category": h2_category,
                    "vmr_total": int(len(vmr) if h2_category == "all" else (vmr["h2_category"] == h2_category).sum()),
                }
            )

    totals = pd.DataFrame(total_rows)
    totals["h2_label"] = totals["h2_category"].map(CATEGORY_LABELS)
    return region_tables, totals


def load_set_counts(output_dir: Path) -> pd.DataFrame:
    set_counts = pd.read_csv(output_dir / "overlap_summary_counts.tsv", sep="\t")
    set_counts = set_counts.loc[set_counts["flag"] == OVERLAP_FLAG].copy()
    set_counts = set_counts.loc[set_counts["h2_category"].isin(CATEGORY_ORDER)].copy()

    if len(set_counts) != len(CATEGORY_ORDER) * len(SET_ORDER):
        raise ValueError("Unexpected number of F_0.25 set counts for the manuscript summary.")

    set_counts["h2_label"] = set_counts["h2_category"].map(CATEGORY_LABELS)
    set_counts["class_order"] = set_counts["h2_category"].map({name: idx + 1 for idx, name in enumerate(CATEGORY_ORDER)})
    set_counts["set_order"] = set_counts["set"].map({name: idx + 1 for idx, name in enumerate(SET_ORDER)})
    set_counts["set_label"] = set_counts["set"]

    for field in ("includes_caudate", "includes_dlpfc", "includes_hippocampus"):
        set_counts[field] = set_counts["set"].map(lambda name: SET_MEMBERSHIP[name][field])

    set_counts["combination_size"] = (
        set_counts["includes_caudate"] + set_counts["includes_dlpfc"] + set_counts["includes_hippocampus"]
    )

    set_counts = set_counts.sort_values(["class_order", "set_order"]).reset_index(drop=True)
    return set_counts[
        [
            "flag",
            "h2_category",
            "h2_label",
            "class_order",
            "set",
            "set_label",
            "set_order",
            "count",
            "combination_size",
            "includes_caudate",
            "includes_dlpfc",
            "includes_hippocampus",
        ]
    ]


def validate_region_totals(region_totals: pd.DataFrame) -> None:
    expected_rows = len(REGIONS) * len(CATEGORY_ORDER)
    if len(region_totals) != expected_rows:
        raise ValueError("Unexpected number of region/class totals derived from the elastic-net summary tables.")

    missing = region_totals.loc[region_totals["vmr_total"].isna()]
    if not missing.empty:
        raise ValueError("Missing VMR totals detected in the main-analysis class summary.")


def load_pairwise_overlaps(percent_overlap_dir: Path, tissue1: str, tissue2: str, h2_category: str) -> pd.DataFrame:
    overlap_file = percent_overlap_dir / f"{tissue1}_{tissue2}_overlap_{h2_category}.tsv"
    overlap = pd.read_csv(overlap_file, sep="\t", dtype={"chromA": str, "chromB": str})
    overlap["interval_id_a"] = interval_id(overlap, "chromA", "startA", "endA")
    overlap["interval_id_b"] = interval_id(overlap, "chromB", "startB", "endB")
    overlap = overlap.drop_duplicates(subset=["interval_id_a", "interval_id_b"]).copy()
    overlap = overlap.loc[overlap["reciprocal"] >= RECIPROCAL_THRESHOLD].copy()
    return overlap


def summarize_reciprocal_overlap(percent_overlap_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for pair_order, (tissue1, tissue2) in enumerate(PAIR_ORDER, start=1):
        for h2_category in CATEGORY_ORDER:
            overlap = load_pairwise_overlaps(percent_overlap_dir, tissue1, tissue2, h2_category)
            rows.append(
                {
                    "flag": OVERLAP_FLAG,
                    "reciprocal_threshold": RECIPROCAL_THRESHOLD,
                    "h2_category": h2_category,
                    "h2_label": CATEGORY_LABELS[h2_category],
                    "class_order": CATEGORY_ORDER.index(h2_category) + 1,
                    "tissue1": tissue1,
                    "tissue2": tissue2,
                    "pair_order": pair_order,
                    "pair_label": PAIR_LABELS[(tissue1, tissue2)],
                    "pair_label_multiline": PAIR_LABELS_MULTILINE[(tissue1, tissue2)],
                    "n_pairs": int(len(overlap)),
                    "n_unique_tissue1": int(overlap["interval_id_a"].nunique()),
                    "n_unique_tissue2": int(overlap["interval_id_b"].nunique()),
                    "median_reciprocal": float(overlap["reciprocal"].median()),
                    "mean_reciprocal": float(overlap["reciprocal"].mean()),
                    "q25_reciprocal": float(overlap["reciprocal"].quantile(0.25)),
                    "q75_reciprocal": float(overlap["reciprocal"].quantile(0.75)),
                    "pct_pairs_ge_50": float((overlap["reciprocal"] >= 50).mean() * 100),
                    "pct_pairs_ge_75": float((overlap["reciprocal"] >= 75).mean() * 100),
                    "pct_pairs_ge_90": float((overlap["reciprocal"] >= 90).mean() * 100),
                }
            )

    summary = pd.DataFrame(rows)
    return summary.sort_values(["class_order", "pair_order"]).reset_index(drop=True)


def summarize_h2_concordance(percent_overlap_dir: Path, region_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for pair_order, (tissue1, tissue2) in enumerate(PAIR_ORDER, start=1):
        region_a = region_tables[tissue1][["interval_id", "h2_unscaled", "h2_category"]].rename(
            columns={"interval_id": "interval_id_a", "h2_unscaled": "h2_tissue1", "h2_category": "h2_category_tissue1"}
        )
        region_b = region_tables[tissue2][["interval_id", "h2_unscaled", "h2_category"]].rename(
            columns={"interval_id": "interval_id_b", "h2_unscaled": "h2_tissue2", "h2_category": "h2_category_tissue2"}
        )

        for h2_category in CATEGORY_ORDER:
            overlap = load_pairwise_overlaps(percent_overlap_dir, tissue1, tissue2, h2_category)
            merged = overlap.merge(region_a, on="interval_id_a", how="left").merge(region_b, on="interval_id_b", how="left")

            missing_pairs = int(merged[["h2_tissue1", "h2_tissue2"]].isna().any(axis=1).sum())
            if h2_category == "all" and missing_pairs > 0:
                merged = merged.dropna(subset=["h2_tissue1", "h2_tissue2"]).copy()

            if merged[["h2_tissue1", "h2_tissue2"]].isna().any().any():
                raise ValueError(f"Missing h2 values after joining {tissue1}-{tissue2} overlaps for {h2_category}.")

            if h2_category != "all":
                mismatch = merged.loc[
                    (merged["h2_category_tissue1"] != h2_category) | (merged["h2_category_tissue2"] != h2_category)
                ]
                if not mismatch.empty:
                    raise ValueError(
                        f"Category mismatch detected in {tissue1}-{tissue2} overlaps for {h2_category}."
                    )

            spearman = spearmanr(merged["h2_tissue1"], merged["h2_tissue2"], nan_policy="omit")
            rho = float(getattr(spearman, "statistic", spearman[0]))
            p_value = float(getattr(spearman, "pvalue", spearman[1]))

            rows.append(
                {
                    "flag": OVERLAP_FLAG,
                    "reciprocal_threshold": RECIPROCAL_THRESHOLD,
                    "h2_category": h2_category,
                    "h2_label": CATEGORY_LABELS[h2_category],
                    "class_order": CATEGORY_ORDER.index(h2_category) + 1,
                    "tissue1": tissue1,
                    "tissue2": tissue2,
                    "pair_order": pair_order,
                    "pair_label": PAIR_LABELS[(tissue1, tissue2)],
                    "pair_label_multiline": PAIR_LABELS_MULTILINE[(tissue1, tissue2)],
                    "n_pairs": int(len(merged)),
                    "n_pairs_missing_h2": missing_pairs,
                    "n_unique_tissue1": int(merged["interval_id_a"].nunique()),
                    "n_unique_tissue2": int(merged["interval_id_b"].nunique()),
                    "spearman_rho": rho,
                    "spearman_p_value": p_value,
                }
            )

    concordance = pd.DataFrame(rows).sort_values(["class_order", "pair_order"]).reset_index(drop=True)
    concordance["spearman_fdr"] = benjamini_hochberg(concordance["spearman_p_value"])
    return concordance


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    analysis_dir = script_dir.parent
    output_dir = analysis_dir / "_m"
    percent_overlap_dir = output_dir / OVERLAP_FLAG / "percent_overlap"
    manuscript_dir = output_dir / "manuscript_overlap_figure"
    manuscript_dir.mkdir(parents=True, exist_ok=True)

    region_tables, region_totals = load_region_summaries(analysis_dir)
    validate_region_totals(region_totals)
    set_counts = load_set_counts(output_dir)

    reciprocal_summary = summarize_reciprocal_overlap(percent_overlap_dir)
    h2_concordance = summarize_h2_concordance(percent_overlap_dir, region_tables)

    set_counts.to_csv(manuscript_dir / "F_0.25_set_counts.tsv", sep="\t", index=False)
    reciprocal_summary.to_csv(manuscript_dir / "F_0.25_reciprocal_overlap_summary.tsv", sep="\t", index=False)
    h2_concordance.to_csv(manuscript_dir / "F_0.25_h2_concordance_summary.tsv", sep="\t", index=False)

    print("Loaded main-analysis class totals from elastic-net summaries:")
    print(region_totals.to_string(index=False))
    print("Validated category assignments through the joined F_0.25 overlap records.")
    print(reciprocal_summary[["h2_category", "pair_label", "median_reciprocal", "n_pairs"]].to_string(index=False))
    print(h2_concordance[["h2_category", "pair_label", "spearman_rho", "spearman_p_value", "n_pairs"]].to_string(index=False))

    if session_info is not None:
        session_info.show()


if __name__ == "__main__":
    main()
