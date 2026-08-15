#!/usr/bin/env python3
"""External brain meQTL and existing transcription/splicing validation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import canonical_vmr_id, write_tsv  # noqa: E402
from _lib.stats_utils import zscore  # noqa: E402

LINK_REL = {
    "expression": "{region}/AA/expression/nearest_gene_window_250kb/architecture_model_input.tsv",
    "psi": "{region}/AA/psi/window_250kb/architecture_model_input.tsv",
    "expression_abc": "{region}/AA/expression/abc/architecture_model_input.tsv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", required=True)
    parser.add_argument("--primary-bridge", required=True)
    parser.add_argument("--calibrated", required=True)
    parser.add_argument("--external-root", required=True)
    parser.add_argument("--regulatory-context-root", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def _fit_binary(
    frame: pd.DataFrame,
    *,
    outcome: str,
    predictor: str,
    estimator: str,
    analysis: str,
    covariates: list[str],
    sensitivity: str = "primary",
) -> dict:
    columns = [outcome, predictor] + covariates
    use = frame.dropna(subset=columns).copy()
    if len(use) < 50 or use[outcome].nunique() < 2:
        return {
            "analysis": analysis, "estimator": estimator, "sensitivity": sensitivity,
            "n_vmrs": int(len(use)), "error": "insufficient_rows_or_outcome_variation",
        }
    exog = pd.DataFrame({predictor: zscore(use[predictor])}, index=use.index)
    retained_covariates = []
    for column in covariates:
        values = pd.to_numeric(use[column], errors="coerce")
        if values.nunique() <= 1:
            continue
        candidate = values.astype(float) if values.nunique() == 2 else zscore(values)
        trial = pd.concat([exog, candidate.rename(column)], axis=1)
        if np.linalg.matrix_rank(trial.to_numpy(float)) > np.linalg.matrix_rank(exog.to_numpy(float)):
            exog[column] = candidate
            retained_covariates.append(column)
    exog = sm.add_constant(exog, has_constant="add")
    try:
        result = sm.GLM(use[outcome].astype(float), exog, family=sm.families.Binomial()).fit(
            cov_type="HC3", maxiter=250
        )
        estimate = float(result.params[predictor])
        se = float(result.bse[predictor])
        return {
            "analysis": analysis,
            "estimator": estimator,
            "sensitivity": sensitivity,
            "estimate_log_or_per_sd": estimate,
            "or_per_sd": float(np.exp(estimate)),
            "se": se,
            "ci_lower": float(estimate - 1.96 * se),
            "ci_upper": float(estimate + 1.96 * se),
            "pvalue": float(result.pvalues[predictor]),
            "n_vmrs": int(len(use)),
            "n_positive": int(use[outcome].sum()),
            "covariates": ",".join(retained_covariates),
            "converged": bool(result.converged),
            "error": "",
        }
    except Exception as error:  # noqa: BLE001
        return {
            "analysis": analysis, "estimator": estimator, "sensitivity": sensitivity,
            "n_vmrs": int(len(use)), "error": str(error),
        }


def canonicalize_external_vmr_keys(ext: pd.DataFrame, bridge: pd.DataFrame) -> pd.DataFrame:
    """Resolve mixed numeric task IDs and coordinate IDs to bridge coordinates."""
    result = ext.copy()
    raw_key = result["vmr_id"].astype(str).str.strip()
    task_to_coord = dict(zip(
        bridge["vmr_id"].astype(str).str.strip(),
        bridge["vmr_coord_id"].astype(str),
    ))
    coordinate_universe = set(bridge["vmr_coord_id"].astype(str))
    result["vmr_coord_id"] = raw_key.map(task_to_coord)
    coordinate_key = raw_key.str.replace("^chr", "", regex=True)
    missing = result["vmr_coord_id"].isna() & coordinate_key.isin(coordinate_universe)
    result.loc[missing, "vmr_coord_id"] = coordinate_key.loc[missing]
    return result.dropna(subset=["vmr_coord_id"]).copy()


def external_validation(region: str, bridge: pd.DataFrame, root: Path) -> tuple[list[dict], list[dict]]:
    resources = {
        "jaffe_dlpfc_450k_meqtl": ("dlpfc", True),
        "schulz_hippocampus_array_meqtl": ("hippocampus", False),
    }
    model_rows = []
    descriptive = []
    for resource, (target_region, inferential) in resources.items():
        if region != target_region:
            continue
        path = root / "harmonized" / f"{resource}.{region}.vmr_support.tsv.gz"
        if not path.is_file():
            raise SystemExit(f"Missing external support table: {path}")
        ext = pd.read_csv(path, sep="\t")
        # Harmonized resources preserve the source VMR key, which can be either
        # the numeric AA task_id or a genomic coordinate. Canonicalize both
        # representations through the accepted bridge before aggregation.
        ext = canonicalize_external_vmr_keys(ext, bridge)
        source_universe_complete = (
            _as_bool(ext["assay_universe_complete"]).all()
            if "assay_universe_complete" in ext.columns and len(ext)
            else False
        )
        ext["external_meqtl_support"] = pd.to_numeric(
            ext["external_meqtl_support"], errors="coerce"
        ).fillna(0).clip(0, 1).astype(int)
        vmr = ext.groupby("vmr_coord_id", as_index=False).agg(
            n_external_cpgs_assayed=("external_meqtl_support", "size"),
            n_external_cpgs_supported=("external_meqtl_support", "sum"),
            external_meqtl_support=("external_meqtl_support", "max"),
        )
        d = bridge.merge(vmr, on="vmr_coord_id", how="inner", validate="one_to_one")
        d = d[_as_bool(d["common_complete_case"])].copy()
        outcome_variation = bool(len(d) and d["external_meqtl_support"].nunique() >= 2)
        inferentially_usable = bool(inferential and source_universe_complete and outcome_variation)
        if not outcome_variation and len(d):
            universe_note = "no outcome variation among overlapping assayed VMRs"
        elif not source_universe_complete:
            universe_note = "source does not provide a complete assayed negative universe"
        else:
            universe_note = "inferentially usable assayed universe"
        descriptive.append({
            "region": region,
            "resource_id": resource,
            "role": "inferential_independent" if inferential else "descriptive_positive_only",
            "n_overlapping_vmrs": int(len(d)),
            "n_supported_vmrs": int(d["external_meqtl_support"].sum()),
            "support_rate": float(d["external_meqtl_support"].mean()) if len(d) else np.nan,
            "assay_universe_complete": source_universe_complete,
            "outcome_variation": outcome_variation,
            "inferentially_usable": inferentially_usable,
            "universe_note": universe_note,
        })
        if not inferential:
            continue
        candidate = [
            "n_external_cpgs_assayed", "average_cpg_coverage", "mean_cpg_variance", "length"
        ]
        covariates = [
            column for column in candidate
            if column in d.columns and d[column].notna().sum() >= max(50, int(0.8 * len(d)))
        ]
        common = d.dropna(subset=["calibrated_h2", "legacy_h2_unscaled"] + covariates)
        for estimator, predictor in [
            ("calibrated", "calibrated_h2"), ("legacy", "legacy_h2_unscaled")
        ]:
            row = _fit_binary(
                common,
                outcome="external_meqtl_support",
                predictor=predictor,
                estimator=estimator,
                analysis=f"external_{resource}",
                covariates=covariates,
            )
            row.update({"region": region, "resource_id": resource, "role": "inferential_independent"})
            model_rows.append(row)
        no_boundary = common.loc[~_as_bool(common["h2_upper_boundary_hit"])]
        row = _fit_binary(
            no_boundary,
            outcome="external_meqtl_support",
            predictor="calibrated_h2",
            estimator="calibrated",
            analysis=f"external_{resource}",
            covariates=covariates,
            sensitivity="exclude_upper_boundary",
        )
        row.update({"region": region, "resource_id": resource, "role": "inferential_independent"})
        model_rows.append(row)
    return model_rows, descriptive


def transcription_validation(
    region: str, calibrated_path: Path, root: Path
) -> tuple[list[dict], list[dict]]:
    calibrated = pd.read_csv(calibrated_path, sep="\t")
    calibrated = calibrated[
        calibrated["region"].astype(str).str.lower().eq(region)
        & calibrated["calibration_status"].eq("within_domain")
    ].copy()
    calibrated["vmr_coord_id"] = [
        canonical_vmr_id(chrom, start, end)
        for chrom, start, end in zip(
            calibrated["chromosome"], calibrated["start"], calibrated["end"]
        )
    ]
    calibrated["h2_upper_boundary_hit"] = _as_bool(calibrated["h2_upper_boundary_hit"])
    calibrated = calibrated.drop_duplicates("vmr_coord_id", keep="first")
    model_rows = []
    descriptive = []
    for modality, relative in LINK_REL.items():
        path = root / relative.format(region=region)
        if not path.is_file():
            descriptive.append({
                "region": region, "modality": modality, "status": "missing", "path": str(path)
            })
            continue
        links = pd.read_csv(path, sep="\t")
        parts = links["vmr_id"].astype(str).str.extract(r"^chr([^_]+)_([0-9]+)_([0-9]+)$")
        links["vmr_coord_id"] = [
            canonical_vmr_id(chrom, start, end)
            for chrom, start, end in zip(parts[0], parts[1], parts[2])
        ]
        links["tx_associated"] = pd.to_numeric(
            links["any_sig_fdr_05"], errors="coerce"
        ).fillna(0).eq(1).astype(int)
        links["legacy_h2"] = pd.to_numeric(links["h2_unscaled"], errors="coerce")
        links = links.drop_duplicates("vmr_coord_id", keep="first")
        d = links.merge(
            calibrated[["vmr_coord_id", "h2_en_calibrated", "h2_upper_boundary_hit"]].rename(
                columns={"h2_en_calibrated": "calibrated_h2"}
            ), on="vmr_coord_id", how="inner", validate="one_to_one",
        )
        candidate = [
            "n_features_tested", "vmr_length", "min_distance", "methylation_variance", "num_snps"
        ]
        covariates = [
            column for column in candidate
            if column in d.columns and d[column].notna().sum() >= max(50, int(0.8 * len(d)))
        ]
        common = d.dropna(subset=["calibrated_h2", "legacy_h2"] + covariates)
        descriptive.append({
            "region": region,
            "modality": modality,
            "status": "available",
            "n_common_vmrs": int(len(common)),
            "n_tx_associated": int(common["tx_associated"].sum()),
            "path": str(path),
        })
        for estimator, predictor in [("calibrated", "calibrated_h2"), ("legacy", "legacy_h2")]:
            row = _fit_binary(
                common,
                outcome="tx_associated",
                predictor=predictor,
                estimator=estimator,
                analysis=f"transcription_{modality}",
                covariates=covariates,
            )
            row.update({"region": region, "modality": modality})
            model_rows.append(row)
        row = _fit_binary(
            common.loc[~common["h2_upper_boundary_hit"]],
            outcome="tx_associated",
            predictor="calibrated_h2",
            estimator="calibrated",
            analysis=f"transcription_{modality}",
            covariates=covariates,
            sensitivity="exclude_upper_boundary",
        )
        row.update({"region": region, "modality": modality})
        model_rows.append(row)
    return model_rows, descriptive


def main() -> None:
    args = parse_args()
    region = args.region.lower()
    bridge = pd.read_csv(args.primary_bridge, sep="\t")
    external_models, external_description = external_validation(
        region, bridge, Path(args.external_root)
    )
    tx_models, tx_description = transcription_validation(
        region, Path(args.calibrated), Path(args.regulatory_context_root)
    )
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    external_model_columns = [
        "analysis", "estimator", "sensitivity", "estimate_log_or_per_sd", "or_per_sd",
        "se", "ci_lower", "ci_upper", "pvalue", "n_vmrs", "n_positive", "covariates",
        "converged", "error", "region", "resource_id", "role",
    ]
    pd.DataFrame(external_models, columns=external_model_columns).to_csv(
        output / "external_validation_models.tsv", sep="\t", index=False
    )
    write_tsv(
        output / "external_validation_descriptive.tsv",
        external_description,
        fieldnames=[
            "region", "resource_id", "role", "n_overlapping_vmrs", "n_supported_vmrs",
            "support_rate", "assay_universe_complete", "outcome_variation",
            "inferentially_usable", "universe_note",
        ],
    )
    pd.DataFrame(tx_models).to_csv(output / "transcription_validation_models.tsv", sep="\t", index=False)
    write_tsv(output / "transcription_validation_descriptive.tsv", tx_description)
    print(f"Wrote orthogonal validation under {output}")


if __name__ == "__main__":
    main()
