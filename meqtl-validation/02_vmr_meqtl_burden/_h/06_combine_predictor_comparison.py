#!/usr/bin/env python3
"""Combine region-level predictor comparisons, control FDR, and evaluate gates."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import load_yaml, write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--config", default="predictor_comparison.yml")
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def _load_config(value: str) -> dict:
    if "/" not in value:
        return load_yaml(value)
    import yaml

    with Path(value).open() as handle:
        return yaml.safe_load(handle)


def _bh(pvalues: pd.Series) -> pd.Series:
    values = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    result = np.full(len(values), np.nan)
    finite = np.flatnonzero(np.isfinite(values))
    if not len(finite):
        return pd.Series(result, index=pvalues.index)
    ordered = finite[np.argsort(values[finite])]
    ranked = values[ordered] * len(ordered) / np.arange(1, len(ordered) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result[ordered] = np.minimum(ranked, 1.0)
    return pd.Series(result, index=pvalues.index)


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def _grouped_bh(frame: pd.DataFrame, group: str, pvalue: str) -> pd.Series:
    adjusted = pd.Series(np.nan, index=frame.index, dtype=float)
    for _, indices in frame.groupby(group, sort=False).groups.items():
        adjusted.loc[indices] = _bh(frame.loc[indices, pvalue])
    return adjusted


def _read_regions(root: Path, regions: list[str], filename: str) -> pd.DataFrame:
    tables = []
    for region in regions:
        path = root / region / filename
        if not path.is_file():
            raise SystemExit(f"Missing region output: {path}")
        table = pd.read_csv(path, sep="\t")
        if "region" not in table.columns:
            table["region"] = region
        tables.append(table)
    return pd.concat(tables, ignore_index=True, sort=False)


def main() -> None:
    args = parse_args()
    config = _load_config(args.config)
    root = Path(args.input_root)
    regions = list(config["regions"])
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)

    models = _read_regions(root, regions, "burden_model_comparison.tsv")
    models["qvalue"] = _grouped_bh(models, "test_family", "pvalue")
    models.to_csv(output / "burden_model_comparison.tsv", sep="\t", index=False)

    matched = _read_regions(root, regions, "matched_comparison.tsv")
    matched["qvalue"] = _grouped_bh(matched, "analysis", "permutation_pvalue")
    matched.to_csv(output / "matched_comparison.tsv", sep="\t", index=False)

    for filename in [
        "predictor_concordance.tsv", "classification_concordance.tsv",
        "quintile_transition.tsv", "join_qc.tsv", "comparison_model_qc.tsv",
    ]:
        _read_regions(root, regions, filename).to_csv(
            output / filename, sep="\t", index=False
        )
    balance_tables = []
    for region in regions:
        path = root / region / "matched_balance.tsv"
        if path.is_file():
            balance_tables.append(pd.read_csv(path, sep="\t"))
    if balance_tables:
        pd.concat(balance_tables, ignore_index=True, sort=False).to_csv(
            output / "matched_balance.tsv", sep="\t", index=False
        )

    acceptance_cfg = config["acceptance"]
    threshold = float(acceptance_cfg["fdr"])
    calibrated = models[models["test_family"].eq("calibrated_primary")]
    calibrated_regions = int(
        (calibrated["estimate_per_sd"].gt(0) & calibrated["qvalue"].le(threshold)).sum()
    )
    calibrated_matched = matched[
        matched["analysis"].eq("calibrated_positive_vs_nonpositive")
    ]
    matched_regions = int((
        calibrated_matched["mean_difference"].gt(0)
        & calibrated_matched["qvalue"].le(threshold)
        & _as_bool(calibrated_matched["balance_pass"])
    ).sum())
    join_qc = pd.read_csv(output / "join_qc.tsv", sep="\t")
    model_qc = pd.read_csv(output / "comparison_model_qc.tsv", sep="\t")
    zero_computational = bool(join_qc["calibrated_computational_failures"].eq(0).all())
    calibrated_qc_pass = bool(_as_bool(join_qc["calibrated_overall_qc_pass"]).all())
    structural_qc = bool(
        _as_bool(join_qc["join_unique"]).all()
        and _as_bool(model_qc["all_models_converged"]).all()
    )

    criteria = [
        {
            "criterion": "calibrated_burden_positive_fdr_in_required_regions",
            "observed": calibrated_regions,
            "required": int(acceptance_cfg["min_regions_positive_fdr"]),
            "pass": calibrated_regions >= int(acceptance_cfg["min_regions_positive_fdr"]),
        },
        {
            "criterion": "calibrated_matched_contrast_positive_in_required_regions",
            "observed": matched_regions,
            "required": int(acceptance_cfg["min_regions_matched_positive"]),
            "pass": matched_regions >= int(acceptance_cfg["min_regions_matched_positive"]),
        },
        {
            "criterion": "zero_calibrated_computational_failures",
            "observed": int(zero_computational),
            "required": 1,
            "pass": zero_computational,
        },
        {
            "criterion": "calibrated_region_qc_and_comparison_structure_pass",
            "observed": int(calibrated_qc_pass and structural_qc),
            "required": 1,
            "pass": calibrated_qc_pass and structural_qc,
        },
    ]
    phase2_gate_pass = all(bool(row["pass"]) for row in criteria)
    criteria.append({
        "criterion": "annotation_sensitivity_for_full_promotion",
        "observed": "pending",
        "required": "required",
        "pass": False,
    })
    write_tsv(output / "comparison_acceptance.tsv", criteria)
    write_tsv(output / "promotion_status.tsv", [{
        "phase2_calibrated_gate_pass": phase2_gate_pass,
        "full_manuscript_promotion_ready": False,
        "full_promotion_blocker": "key LINE/L1 and repressive-chromatin calibrated sensitivities not yet evaluated",
        "calibrated_regions_positive_fdr": calibrated_regions,
        "calibrated_matched_regions_positive": matched_regions,
    }])
    print(
        f"Phase 2 calibrated gate pass={phase2_gate_pass}; "
        f"model regions={calibrated_regions}, matched regions={matched_regions}"
    )


if __name__ == "__main__":
    main()
