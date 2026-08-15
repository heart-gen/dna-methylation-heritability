#!/usr/bin/env python3
"""Tests for calibrated-estimator promotion follow-up analyses."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "meqtl-validation"))


def _load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ALIGN = _load(
    "followup_align",
    "meqtl-validation/02_vmr_meqtl_burden/_h/08_reaggregate_coordinate_aligned.py",
)
WEIGHT = _load(
    "followup_weight",
    "meqtl-validation/02_vmr_meqtl_burden/_h/09_fit_boundary_and_weighting.py",
)
ANNOTATION = _load(
    "followup_annotation",
    "meqtl-validation/02_vmr_meqtl_burden/_h/10_fit_annotation_sensitivity.py",
)
ORTHOGONAL = _load(
    "followup_orthogonal",
    "meqtl-validation/02_vmr_meqtl_burden/_h/11_fit_orthogonal_validation.py",
)
COMBINE = _load(
    "followup_combine",
    "meqtl-validation/02_vmr_meqtl_burden/_h/12_combine_followup.py",
)


def test_coordinate_assignment_is_exact_and_does_not_force_nearest_interval():
    vmrs = pd.DataFrame({
        "chrom": ["1", "1", "2"],
        "start": [100, 300, 50],
        "end": [200, 350, 75],
        "vmr_coord_id": ["1:100-200", "1:300-350", "2:50-75"],
    })
    points = pd.DataFrame({
        "chrom": ["1", "1", "1", "1", "2", "3"],
        "pos_1based": [100, 200, 201, 325, 75, 60],
    })
    assigned = ALIGN.assign_coordinates(points, vmrs)
    assert assigned.tolist()[:2] == ["1:100-200", "1:100-200"]
    assert pd.isna(assigned.iloc[2])
    assert assigned.iloc[3] == "1:300-350"
    assert assigned.iloc[4] == "2:50-75"
    assert pd.isna(assigned.iloc[5])


def test_reaggregation_denominator_contains_only_tested_rows():
    tested = pd.DataFrame({
        "vmr_coord_id": ["1:100-200", "1:100-200", "1:300-350"],
        "phenotype_id": ["a", "b", "c"],
        "qval": [0.01, 0.5, 0.04],
        "pval_beta": [0.001, 0.2, 0.003],
        "slope": [0.2, -0.1, 0.3],
        "num_var": [20, 30, 40],
        "variant_id": ["v1", "v2", "v3"],
    })
    result = ALIGN._aggregate(tested, "caudate", 0.05)
    first = result[result["vmr_coord_id"].eq("1:100-200")].iloc[0]
    assert first["n_tested_cpgs"] == 2
    assert first["n_cpgs_with_sig_meqtl"] == 1
    assert first["proportion_cpgs_with_sig_meqtl"] == 0.5


def test_overlap_weighting_recovers_balance_and_positive_contrast():
    rng = np.random.default_rng(42)
    n = 600
    x = rng.normal(size=n)
    exposure = x + rng.normal(scale=0.8, size=n) > 0
    outcome = 0.15 + 0.5 * exposure + 0.05 * x + rng.normal(scale=0.05, size=n)
    frame = pd.DataFrame({
        "calibrated_positive_signal": exposure,
        "proportion_cpgs_with_sig_meqtl": outcome,
        "coverage": x,
        "broad_genomic_annotation": np.where(np.arange(n) % 2, "genic", "other"),
    })
    result, balance = WEIGHT.overlap_weighted(frame, ["coverage"])
    assert result["mean_difference"] > 0.4
    assert result["max_abs_smd_after"] < 0.1
    assert balance["smd_overlap_weighted"].abs().max() < 0.1


def test_logistic_helpers_drop_constant_covariates_and_recover_positive_effect():
    rng = np.random.default_rng(7)
    n = 500
    predictor = rng.normal(size=n)
    probability = 1 / (1 + np.exp(-(-1 + predictor)))
    frame = pd.DataFrame({
        "outcome": rng.binomial(1, probability),
        "predictor": predictor,
        "constant": 1,
        "variable": rng.normal(size=n),
    })
    row = ORTHOGONAL._fit_binary(
        frame,
        outcome="outcome",
        predictor="predictor",
        estimator="calibrated",
        analysis="synthetic",
        covariates=["constant", "variable"],
    )
    assert row["estimate_log_or_per_sd"] > 0
    assert "constant" not in row["covariates"]
    assert row["converged"]


def test_external_mixed_task_and_coordinate_keys_are_canonicalized():
    bridge = pd.DataFrame({
        "vmr_id": ["17", "18"],
        "vmr_coord_id": ["1:100-200", "1:300-400"],
    })
    external = pd.DataFrame({
        "vmr_id": ["17", "chr1:300-400", "missing"],
        "external_meqtl_support": [1, 0, 1],
    })
    result = ORTHOGONAL.canonicalize_external_vmr_keys(external, bridge)
    assert result["vmr_coord_id"].tolist() == ["1:100-200", "1:300-400"]


def test_annotation_logistic_uses_common_rows():
    rng = np.random.default_rng(9)
    n = 300
    predictor = rng.normal(size=n)
    outcome = rng.binomial(1, 1 / (1 + np.exp(-predictor)))
    frame = pd.DataFrame({"feature": outcome, "calibrated_h2": predictor, "length": rng.normal(size=n)})
    row = ANNOTATION.fit_logistic(
        frame,
        outcome="feature",
        predictor="calibrated_h2",
        covariates=["length"],
        estimator="calibrated",
        sensitivity="adjusted_prespecified",
    )
    assert row["n_vmrs"] == n
    assert row["estimate_log_or_per_sd"] > 0
    assert row["converged"]


def test_followup_bh_preserves_missing_values():
    values = pd.Series([0.01, np.nan, 0.04], index=[4, 2, 9])
    adjusted = COMBINE._bh(values)
    assert np.allclose(adjusted.loc[[4, 9]], [0.02, 0.04])
    assert np.isnan(adjusted.loc[2])
