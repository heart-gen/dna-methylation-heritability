#!/usr/bin/env python3
"""Combine calibrated follow-up analyses and evaluate promotion gates."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--followup-config", required=True)
    parser.add_argument("--base-promotion", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def _bh(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    result = np.full(len(values), np.nan)
    finite = np.flatnonzero(np.isfinite(values))
    if not len(finite):
        return pd.Series(result, index=series.index)
    order = finite[np.argsort(values[finite])]
    adjusted = values[order] * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result[order] = np.minimum(adjusted, 1.0)
    return pd.Series(result, index=series.index)


def _grouped_bh(frame: pd.DataFrame, groups: list[str], pvalue: str) -> pd.Series:
    result = pd.Series(np.nan, index=frame.index, dtype=float)
    for _, indices in frame.groupby(groups, dropna=False, sort=False).groups.items():
        result.loc[indices] = _bh(frame.loc[indices, pvalue])
    return result


def _regions(root: Path, regions: list[str], subdir: str, filename: str) -> pd.DataFrame:
    frames = []
    for region in regions:
        path = root / subdir / region / filename
        if not path.is_file():
            raise SystemExit(f"Missing follow-up result: {path}")
        frame = pd.read_csv(path, sep="\t")
        if "region" not in frame.columns:
            frame["region"] = region
        frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False)


def main() -> None:
    args = parse_args()
    with Path(args.followup_config).open() as handle:
        config = yaml.safe_load(handle)
    regions = list(config["regions"])
    root = Path(args.input_root)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    threshold = float(config["acceptance"]["fdr"])

    aligned = _regions(root, regions, "coordinate_aligned", "burden_model_comparison.tsv")
    aligned["qvalue"] = _grouped_bh(aligned, ["test_family"], "pvalue")
    aligned.to_csv(output / "coordinate_aligned_models.tsv", sep="\t", index=False)
    _regions(root, regions, "coordinate_aligned", "coordinate_alignment_qc.tsv").to_csv(
        output / "coordinate_alignment_qc.tsv", sep="\t", index=False
    )
    _regions(root, regions, "coordinate_aligned", "join_qc.tsv").to_csv(
        output / "coordinate_aligned_join_qc.tsv", sep="\t", index=False
    )

    boundary = _regions(root, regions, "sensitivities", "boundary_excluded_models.tsv")
    boundary["qvalue"] = _grouped_bh(
        boundary, ["vmr_definition", "test_family"], "pvalue"
    )
    boundary.to_csv(output / "boundary_excluded_models.tsv", sep="\t", index=False)
    weighting = _regions(root, regions, "sensitivities", "overlap_weighted_burden.tsv")
    weighting["qvalue"] = _bh(weighting["pvalue"])
    weighting.to_csv(output / "overlap_weighted_burden.tsv", sep="\t", index=False)
    _regions(root, regions, "sensitivities", "overlap_weighted_balance.tsv").to_csv(
        output / "overlap_weighted_balance.tsv", sep="\t", index=False
    )

    annotation = _regions(root, regions, "annotation", "annotation_models.tsv")
    annotation["qvalue"] = _grouped_bh(
        annotation, ["feature", "estimator", "sensitivity"], "pvalue"
    )
    annotation.to_csv(output / "annotation_models.tsv", sep="\t", index=False)
    annotation_matched = _regions(root, regions, "annotation", "annotation_matched.tsv")
    annotation_matched["qvalue"] = _grouped_bh(
        annotation_matched, ["feature", "estimator"], "permutation_pvalue"
    )
    annotation_matched.to_csv(output / "annotation_matched.tsv", sep="\t", index=False)
    _regions(root, regions, "annotation", "annotation_join_qc.tsv").to_csv(
        output / "annotation_join_qc.tsv", sep="\t", index=False
    )

    external = _regions(root, regions, "orthogonal", "external_validation_models.tsv")
    if len(external):
        external["qvalue"] = _grouped_bh(external, ["analysis", "sensitivity"], "pvalue")
    external.to_csv(output / "external_validation_models.tsv", sep="\t", index=False)
    _regions(root, regions, "orthogonal", "external_validation_descriptive.tsv").to_csv(
        output / "external_validation_descriptive.tsv", sep="\t", index=False
    )
    transcription = _regions(root, regions, "orthogonal", "transcription_validation_models.tsv")
    transcription["qvalue"] = _grouped_bh(
        transcription, ["analysis", "estimator", "sensitivity"], "pvalue"
    )
    transcription.to_csv(output / "transcription_validation_models.tsv", sep="\t", index=False)
    _regions(root, regions, "orthogonal", "transcription_validation_descriptive.tsv").to_csv(
        output / "transcription_validation_descriptive.tsv", sep="\t", index=False
    )

    base = pd.read_csv(args.base_promotion, sep="\t")
    phase2_pass = bool(_as_bool(base["phase2_calibrated_gate_pass"]).iloc[0])
    aligned_primary = aligned[aligned["test_family"].eq("calibrated_primary")]
    aligned_regions = int((
        aligned_primary["estimate_per_sd"].gt(0) & aligned_primary["qvalue"].le(threshold)
    ).sum())
    boundary_primary = boundary[
        boundary["vmr_definition"].eq("aa_exact_intersection")
        & boundary["test_family"].eq("calibrated_primary")
    ]
    boundary_regions = int((
        boundary_primary["estimate_per_sd"].gt(0) & boundary_primary["qvalue"].le(threshold)
    ).sum())
    hip_weight = weighting[weighting["region"].eq("hippocampus")]
    hip_weight_pass = bool(len(hip_weight) == 1 and (
        hip_weight["mean_difference"].iloc[0] > 0
        and hip_weight["qvalue"].iloc[0] <= threshold
        and _as_bool(hip_weight["balance_pass"]).iloc[0]
    ))

    annotation_required = int(config["acceptance"]["min_annotation_regions_per_feature"])
    annotation_feature_pass = {}
    restriction_names = ["high_mappability", "snp_proximity_excluded", "segdup_excluded"]
    for feature in config["annotation"]["features"]:
        adjusted = annotation[
            annotation["feature"].eq(feature)
            & annotation["estimator"].eq("calibrated")
            & annotation["sensitivity"].eq("adjusted_prespecified")
        ]
        significant_regions = set(adjusted.loc[
            adjusted["estimate_log_or_per_sd"].gt(0) & adjusted["qvalue"].le(threshold), "region"
        ])
        restriction = annotation[
            annotation["feature"].eq(feature)
            & annotation["estimator"].eq("calibrated")
            & annotation["sensitivity"].isin(restriction_names)
            & annotation["estimate_log_or_per_sd"].notna()
        ]
        direction_by_region = restriction.groupby("region")["estimate_log_or_per_sd"].apply(
            lambda values: bool(values.gt(0).all())
        )
        robust_regions = {region for region in significant_regions if direction_by_region.get(region, False)}
        annotation_feature_pass[feature] = {
            "n_regions": len(robust_regions),
            "pass": len(robust_regions) >= annotation_required,
        }

    jaffe = external[
        external["analysis"].eq("external_jaffe_dlpfc_450k_meqtl")
        & external["estimator"].eq("calibrated")
        & external["sensitivity"].eq("primary")
    ]
    external_pass = bool(len(jaffe) == 1 and (
        jaffe["estimate_log_or_per_sd"].iloc[0] > 0
        and jaffe["qvalue"].iloc[0] <= threshold
    ))

    criteria = [
        {"criterion": "base_phase2_calibrated_gate", "observed": int(phase2_pass), "required": 1, "pass": phase2_pass},
        {
            "criterion": "coordinate_aligned_burden_positive_regions",
            "observed": aligned_regions,
            "required": int(config["acceptance"]["min_coordinate_aligned_regions"]),
            "pass": aligned_regions >= int(config["acceptance"]["min_coordinate_aligned_regions"]),
        },
        {
            "criterion": "boundary_excluded_burden_positive_regions",
            "observed": boundary_regions,
            "required": int(config["acceptance"]["min_boundary_excluded_regions"]),
            "pass": boundary_regions >= int(config["acceptance"]["min_boundary_excluded_regions"]),
        },
        {"criterion": "hippocampus_overlap_weighting_positive_balanced", "observed": int(hip_weight_pass), "required": 1, "pass": hip_weight_pass},
    ]
    for feature, result in annotation_feature_pass.items():
        criteria.append({
            "criterion": f"calibrated_annotation_{feature}_robust_regions",
            "observed": result["n_regions"],
            "required": annotation_required,
            "pass": result["pass"],
        })
    criteria.append({
        "criterion": "independent_jaffe_dlpfc_external_support",
        "observed": int(external_pass),
        "required": 1,
        "pass": external_pass,
    })
    full_pass = all(bool(row["pass"]) for row in criteria)
    pd.DataFrame(criteria).to_csv(output / "followup_acceptance.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "base_phase2_gate_pass": phase2_pass,
        "full_manuscript_promotion_ready": full_pass,
        "recommended_manuscript_axis": "calibrated_primary_legacy_sensitivity" if full_pass else "legacy_locked_calibrated_sensitivity",
        "coordinate_aligned_positive_regions": aligned_regions,
        "boundary_excluded_positive_regions": boundary_regions,
        "hippocampus_weighting_pass": hip_weight_pass,
        "external_jaffe_pass": external_pass,
        "failed_criteria": ",".join(row["criterion"] for row in criteria if not row["pass"]),
    }]).to_csv(output / "promotion_status.tsv", sep="\t", index=False)
    print(f"Full calibrated-estimator promotion ready={full_pass}")


if __name__ == "__main__":
    main()

