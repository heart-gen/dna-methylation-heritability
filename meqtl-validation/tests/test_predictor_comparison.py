#!/usr/bin/env python3
"""Unit tests for the calibrated-versus-legacy Phase 2 comparison."""

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


BRIDGE = _load(
    "predictor_bridge",
    "meqtl-validation/02_vmr_meqtl_burden/_h/04_build_predictor_bridge.py",
)
COMPARE = _load(
    "predictor_compare",
    "meqtl-validation/02_vmr_meqtl_burden/_h/05_fit_predictor_comparison.py",
)
COMBINE = _load(
    "predictor_combine",
    "meqtl-validation/02_vmr_meqtl_burden/_h/06_combine_predictor_comparison.py",
)


def _synthetic_inputs(n: int = 400):
    rng = np.random.default_rng(20260808)
    calibrated_h2 = rng.uniform(0, 0.8, n)
    legacy_h2 = np.clip(0.45 * calibrated_h2 + rng.normal(0, 0.16, n), 0, None)
    probability = 1 / (1 + np.exp(-(-2.0 + 2.0 * calibrated_h2)))
    tested = rng.integers(4, 15, n)
    significant = rng.binomial(tested, probability)
    starts = np.arange(n) * 1000 + 100
    ends = starts + 250
    coords = [f"1:{start}-{end}" for start, end in zip(starts, ends)]
    burden = pd.DataFrame({
        "vmr_id": [str(i + 1) for i in range(n)],
        "vmr_coord_id": coords,
        "analysis_schema_version": 2,
        "tech_join_source": "interval_id",
        "n_tested_cpgs": tested,
        "n_cpgs_with_sig_meqtl": significant,
        "proportion_cpgs_with_sig_meqtl": significant / tested,
        "local_predictability": legacy_h2,
        "average_cpg_coverage": rng.normal(15, 2, n),
        "mean_cpg_variance": rng.uniform(0.01, 0.08, n),
        "vmr_mean_methylation": rng.uniform(0.1, 0.9, n),
        "length": ends - starts,
        "cpg_density": tested / (ends - starts),
        "mean_num_tested_snps_per_cpg": rng.integers(20, 200, n),
        "umap_k24_mean": rng.uniform(0.85, 1.0, n),
        "line_l1_frac": rng.uniform(0, 0.3, n),
        "segdup_frac": 0.0,
        "blacklist_frac": 0.0,
        "genomic_annotation": np.resize(
            np.array(["promoter", "CpG_island", "intron", "intergenic"]), n
        ),
    })
    legacy = pd.DataFrame({
        "task_id": [str(i + 1) for i in range(n)],
        "chrom": "1",
        "start": starts,
        "end": ends,
        "h2_unscaled": legacy_h2,
        "r_squared_cv": np.where(legacy_h2 > 0.05, 0.5, 0.2),
    })
    calibrated = pd.DataFrame({
        "region": "caudate",
        "chromosome": "1",
        "start": starts,
        "end": ends,
        "h2_en_calibrated": calibrated_h2,
        "positive_signal": calibrated_h2 > 0.25,
        "calibration_status": "within_domain",
        "h2_upper_boundary_hit": False,
        "rho2_oof": calibrated_h2 * 0.7,
        "r2_oof": calibrated_h2 * 0.6,
    })
    qc = pd.DataFrame({
        "region": ["caudate"],
        "overall_qc_pass": [True],
        "computational_failed_tasks": [0],
    })
    return burden, legacy, calibrated, qc


def test_bridge_is_exact_and_common_set_is_identical():
    burden, legacy, calibrated, qc = _synthetic_inputs()
    config = {
        "legacy": {
            "estimate_column": "h2_unscaled",
            "prediction_column": "r_squared_cv",
            "prediction_min": 0.3,
            "high_estimate_min": 0.1,
        },
        "calibrated": {
            "estimate_column": "h2_en_calibrated",
            "signal_column": "positive_signal",
            "required_status": "within_domain",
            "exclude_upper_boundary_primary": False,
        },
    }
    result, summary = BRIDGE.build_bridge(
        burden, legacy, calibrated, qc, region="caudate", config=config
    )
    assert result["vmr_coord_id"].is_unique
    assert result["common_complete_case"].all()
    assert summary["n_common_complete_case"] == len(result)
    assert np.allclose(result["legacy_h2_unscaled"], burden["local_predictability"])


def test_adjusted_models_use_one_common_sample_and_recover_positive_signal():
    burden, legacy, calibrated, qc = _synthetic_inputs()
    bridge_config = {
        "legacy": {
            "estimate_column": "h2_unscaled", "prediction_column": "r_squared_cv",
            "prediction_min": 0.3, "high_estimate_min": 0.1,
        },
        "calibrated": {
            "estimate_column": "h2_en_calibrated", "signal_column": "positive_signal",
            "required_status": "within_domain", "exclude_upper_boundary_primary": False,
        },
    }
    bridge, _ = BRIDGE.build_bridge(
        burden, legacy, calibrated, qc, region="caudate", config=bridge_config
    )
    model_config = {
        "model": {
            "min_covariate_complete_fraction": 0.9,
            "candidate_covariates": ["n_tested_cpgs", "average_cpg_coverage"],
            "annotation_indicators": [
                "annotation_promoter", "annotation_cpg_island", "annotation_gene_body"
            ],
        }
    }
    common, covariates = COMPARE.prepare_common_set(bridge, model_config)
    results = pd.DataFrame(COMPARE.fit_model_comparison(common, covariates))
    assert results["n_vmrs"].nunique() == 1
    calibrated_row = results[results["test_family"].eq("calibrated_primary")].iloc[0]
    assert calibrated_row["estimate_per_sd"] > 0
    assert bool(calibrated_row["converged"])


def test_text_false_is_not_truthy_and_grouped_bh_preserves_indices():
    parsed = COMPARE._as_bool(pd.Series(["True", "False", "0", "yes", "NA"]))
    assert parsed.tolist() == [True, False, False, True, False]
    frame = pd.DataFrame({
        "family": ["a", "b", "a", "b"],
        "p": [0.01, 0.04, 0.03, 0.01],
    }, index=[4, 1, 7, 2])
    adjusted = COMBINE._grouped_bh(frame, "family", "p")
    assert adjusted.index.equals(frame.index)
    assert np.allclose(adjusted.loc[[4, 7]], [0.02, 0.03])
    assert np.allclose(adjusted.loc[[1, 2]], [0.04, 0.02])


def test_binary_matching_is_deterministic():
    burden, _, calibrated, _ = _synthetic_inputs()
    data = burden.assign(exposed=calibrated["positive_signal"])
    kwargs = dict(
        exposure="exposed",
        outcome="proportion_cpgs_with_sig_meqtl",
        numeric_covariates=["average_cpg_coverage", "mean_cpg_variance"],
        exact_covariates=["genomic_annotation"],
        caliper_sd=0.5,
        seed=99,
    )
    first, _, first_meta = COMPARE.greedy_binary_propensity_pairs(data, **kwargs)
    second, _, second_meta = COMPARE.greedy_binary_propensity_pairs(data, **kwargs)
    pd.testing.assert_frame_equal(first.reset_index(drop=True), second.reset_index(drop=True))
    assert first_meta["n_pairs"] == second_meta["n_pairs"] > 0

