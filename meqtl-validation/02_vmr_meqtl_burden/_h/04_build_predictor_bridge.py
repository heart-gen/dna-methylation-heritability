#!/usr/bin/env python3
"""Join repaired VMR burden, legacy predictability, and calibrated estimates."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import (  # noqa: E402
    ANALYSIS_SCHEMA_VERSION,
    canonical_vmr_id,
    load_paths,
    load_yaml,
    parse_vmr_coordinate,
    write_tsv,
)

COMPARISON_SCHEMA_VERSION = 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", required=True)
    parser.add_argument("--burden", required=True)
    parser.add_argument("--legacy", default="")
    parser.add_argument("--calibrated", default="")
    parser.add_argument("--calibrated-qc", default="")
    parser.add_argument("--config", default="predictor_comparison.yml")
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def _path(project: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else project / path


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_unique(frame: pd.DataFrame, column: str, label: str) -> None:
    duplicated = frame[column].notna() & frame[column].duplicated(keep=False)
    if duplicated.any():
        examples = ", ".join(frame.loc[duplicated, column].astype(str).head(5))
        raise SystemExit(f"{label} has duplicate {column} values: {examples}")


def _canonicalize_existing(value: object) -> str | None:
    parsed = parse_vmr_coordinate(value)
    return canonical_vmr_id(*parsed) if parsed is not None else None


def build_bridge(
    burden: pd.DataFrame,
    legacy: pd.DataFrame,
    calibrated: pd.DataFrame,
    calibrated_qc: pd.DataFrame,
    *,
    region: str,
    config: dict,
) -> tuple[pd.DataFrame, dict]:
    if "analysis_schema_version" not in burden.columns or not burden[
        "analysis_schema_version"
    ].eq(ANALYSIS_SCHEMA_VERSION).all():
        raise SystemExit("Burden input is stale; schema-v2 aggregation is required")
    if "tech_join_source" not in burden.columns:
        raise SystemExit("Burden input lacks the required repaired technical join")

    qc_region = calibrated_qc[
        calibrated_qc["region"].astype(str).str.lower().eq(region)
    ].copy()
    if len(qc_region) != 1:
        raise SystemExit(f"Expected exactly one calibrated QC row for {region}")
    qc_pass = bool(_as_bool(qc_region["overall_qc_pass"]).iloc[0])
    computational_failures = int(qc_region["computational_failed_tasks"].iloc[0])
    if not qc_pass or computational_failures != 0:
        raise SystemExit(
            f"Accepted calibrated input is required for {region}: "
            f"overall_qc_pass={qc_pass}, computational_failed_tasks={computational_failures}"
        )

    legacy_cfg = config["legacy"]
    calibrated_cfg = config["calibrated"]
    legacy_estimate = legacy_cfg["estimate_column"]
    legacy_prediction = legacy_cfg["prediction_column"]
    for column in ["chrom", "start", "end", legacy_estimate, legacy_prediction]:
        if column not in legacy.columns:
            raise SystemExit(f"Legacy input is missing {column}")
    legacy = legacy.copy()
    legacy["vmr_coord_id"] = [
        canonical_vmr_id(chrom, start, end)
        for chrom, start, end in zip(legacy["chrom"], legacy["start"], legacy["end"])
    ]
    _require_unique(legacy, "vmr_coord_id", "legacy input")
    legacy_keep = legacy[["vmr_coord_id", legacy_estimate, legacy_prediction]].rename(
        columns={
            legacy_estimate: "legacy_h2_unscaled",
            legacy_prediction: "legacy_r_squared_cv",
        }
    )

    calibrated_estimate = calibrated_cfg["estimate_column"]
    signal_column = calibrated_cfg["signal_column"]
    calibrated_required = [
        "region", "chromosome", "start", "end", calibrated_estimate,
        signal_column, "calibration_status", "h2_upper_boundary_hit",
        "rho2_oof", "r2_oof",
    ]
    missing = [column for column in calibrated_required if column not in calibrated.columns]
    if missing:
        raise SystemExit(f"Calibrated input is missing: {', '.join(missing)}")
    calibrated = calibrated[
        calibrated["region"].astype(str).str.lower().eq(region)
    ].copy()
    calibrated["vmr_coord_id"] = [
        canonical_vmr_id(chrom, start, end)
        for chrom, start, end in zip(
            calibrated["chromosome"], calibrated["start"], calibrated["end"]
        )
    ]
    _require_unique(calibrated, "vmr_coord_id", "calibrated input")
    calibrated_keep = calibrated[
        ["vmr_coord_id", calibrated_estimate, signal_column, "calibration_status",
         "h2_upper_boundary_hit", "rho2_oof", "r2_oof"]
    ].rename(
        columns={
            calibrated_estimate: "calibrated_h2",
            signal_column: "calibrated_positive_signal",
            "rho2_oof": "calibrated_rho2_oof",
            "r2_oof": "calibrated_r2_oof",
        }
    )
    calibrated_keep["calibrated_positive_signal"] = _as_bool(
        calibrated_keep["calibrated_positive_signal"]
    )
    calibrated_keep["h2_upper_boundary_hit"] = _as_bool(
        calibrated_keep["h2_upper_boundary_hit"]
    )

    burden = burden.copy()
    burden["vmr_id"] = burden["vmr_id"].astype(str)
    if "vmr_coord_id" not in burden.columns:
        burden["vmr_coord_id"] = burden["vmr_id"].map(_canonicalize_existing)
    else:
        burden["vmr_coord_id"] = burden["vmr_coord_id"].map(_canonicalize_existing)
    task_to_coord = dict(zip(legacy.get("task_id", pd.Series(dtype=str)).astype(str), legacy["vmr_coord_id"]))
    missing_coord = burden["vmr_coord_id"].isna()
    burden.loc[missing_coord, "vmr_coord_id"] = burden.loc[
        missing_coord, "vmr_id"
    ].map(task_to_coord)
    _require_unique(burden, "vmr_coord_id", "burden input")

    bridge = burden.merge(legacy_keep, on="vmr_coord_id", how="left", validate="one_to_one")
    bridge = bridge.merge(
        calibrated_keep, on="vmr_coord_id", how="left", validate="one_to_one"
    )
    if "local_predictability" in bridge.columns:
        aggregate_score = pd.to_numeric(bridge["local_predictability"], errors="coerce")
        legacy_score = pd.to_numeric(bridge["legacy_h2_unscaled"], errors="coerce")
        comparable = aggregate_score.notna() & legacy_score.notna()
        if comparable.any() and not np.allclose(
            aggregate_score[comparable], legacy_score[comparable], rtol=1e-10, atol=1e-12
        ):
            raise SystemExit("Aggregated and source legacy predictability values disagree")

    bridge["legacy_high"] = (
        pd.to_numeric(bridge["legacy_r_squared_cv"], errors="coerce")
        > float(legacy_cfg["prediction_min"])
    ) & (
        pd.to_numeric(bridge["legacy_h2_unscaled"], errors="coerce")
        >= float(legacy_cfg["high_estimate_min"])
    )
    bridge["legacy_low"] = (
        pd.to_numeric(bridge["legacy_r_squared_cv"], errors="coerce")
        > float(legacy_cfg["prediction_min"])
    ) & (
        pd.to_numeric(bridge["legacy_h2_unscaled"], errors="coerce")
        < float(legacy_cfg["high_estimate_min"])
    )
    bridge["legacy_class"] = np.select(
        [bridge["legacy_high"], bridge["legacy_low"]],
        ["high", "low"],
        default="low_prediction",
    )
    required_status = str(calibrated_cfg["required_status"])
    bridge["calibrated_within_domain"] = bridge["calibration_status"].eq(required_status)
    bridge["common_complete_case"] = (
        pd.to_numeric(bridge["legacy_h2_unscaled"], errors="coerce").notna()
        & pd.to_numeric(bridge["calibrated_h2"], errors="coerce").notna()
        & bridge["calibrated_within_domain"]
        & pd.to_numeric(bridge["n_tested_cpgs"], errors="coerce").gt(0)
    )
    if bool(calibrated_cfg.get("exclude_upper_boundary_primary", False)):
        bridge["common_complete_case"] &= ~_as_bool(bridge["h2_upper_boundary_hit"])
    bridge["predictor_comparison_schema_version"] = COMPARISON_SCHEMA_VERSION

    qc = {
        "region": region,
        "analysis_schema_version": ANALYSIS_SCHEMA_VERSION,
        "predictor_comparison_schema_version": COMPARISON_SCHEMA_VERSION,
        "n_burden_vmrs": int(len(bridge)),
        "n_with_legacy": int(bridge["legacy_h2_unscaled"].notna().sum()),
        "n_with_calibrated": int(bridge["calibrated_h2"].notna().sum()),
        "n_within_calibration_domain": int(bridge["calibrated_within_domain"].sum()),
        "n_common_complete_case": int(bridge["common_complete_case"].sum()),
        "n_calibrated_boundary_hits": int(_as_bool(bridge["h2_upper_boundary_hit"]).sum()),
        "n_with_technical_join": int(bridge["tech_join_source"].notna().sum()),
        "calibrated_overall_qc_pass": qc_pass,
        "calibrated_computational_failures": computational_failures,
        "join_unique": True,
    }
    return bridge, qc


def main() -> None:
    args = parse_args()
    paths = load_paths()
    if "/" not in args.config:
        config = load_yaml(args.config)
    else:
        with Path(args.config).open() as handle:
            config = yaml.safe_load(handle)
    project = Path(paths["project_root"])
    region = args.region.lower()
    legacy_path = Path(args.legacy) if args.legacy else (
        project / paths["local_predictability_summary_template"].format(region=region)
    )
    calibrated_path = Path(args.calibrated) if args.calibrated else _path(
        project, config["calibrated"]["estimates"]
    )
    calibrated_qc_path = Path(args.calibrated_qc) if args.calibrated_qc else _path(
        project, config["calibrated"]["qc"]
    )
    burden_path = Path(args.burden)
    for path in [burden_path, legacy_path, calibrated_path, calibrated_qc_path]:
        if not path.is_file():
            raise SystemExit(f"Required input is missing: {path}")

    bridge, qc = build_bridge(
        pd.read_csv(burden_path, sep="\t"),
        pd.read_csv(legacy_path, sep="\t"),
        pd.read_csv(calibrated_path, sep="\t"),
        pd.read_csv(calibrated_qc_path, sep="\t"),
        region=region,
        config=config,
    )
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bridge_path = output_dir / "predictor_bridge.tsv.gz"
    bridge.to_csv(bridge_path, sep="\t", index=False, compression="gzip")
    qc.update({
        "burden_sha256": _sha256(burden_path),
        "legacy_sha256": _sha256(legacy_path),
        "calibrated_sha256": _sha256(calibrated_path),
        "calibrated_qc_sha256": _sha256(calibrated_qc_path),
        "output": str(bridge_path),
    })
    write_tsv(output_dir / "join_qc.tsv", [qc])
    print(f"Wrote {bridge_path}")


if __name__ == "__main__":
    main()
