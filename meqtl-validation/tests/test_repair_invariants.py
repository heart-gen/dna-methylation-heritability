#!/usr/bin/env python3
"""Regression tests for the meQTL-validation repair invariants."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
MEQTL = ROOT / "meqtl-validation"
sys.path.insert(0, str(MEQTL))

from _lib.io_utils import canonical_vmr_id, parse_vmr_coordinate  # noqa: E402
from _lib.stats_utils import paired_randomization_pvalue  # noqa: E402


def load_script(name: str, relative: str):
    path = MEQTL / relative
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


CROSS = load_script("cross_region_repair", "04_cross_region_sharing/_h/01_cross_region_sharing.py")
EXTERNAL = load_script("external_repair", "03_external_meqtl_validation/_h/01_test_external_support.py")
DOWNSAMPLE = load_script(
    "downsample_repair",
    "10_downsampling_caudate/_h/04_summarize_tensorqtl_downsample.py",
)
HARMONIZE = load_script(
    "external_harmonize_repair",
    "03_external_meqtl_validation/_h/03_harmonize_external_meqtls.py",
)


class RepairInvariantTests(unittest.TestCase):
    def test_numeric_task_id_is_not_a_coordinate(self):
        self.assertIsNone(parse_vmr_coordinate("1001"))
        self.assertEqual(canonical_vmr_id("chr1", 10, 20), "1:10-20")

    def test_vmr_sharing_uses_coordinate_not_region_local_task_id(self):
        def burden(coord, task):
            return pd.DataFrame({
                "vmr_coord_id": coord,
                "vmr_id": task,
                "meqtl_supported": [True] * len(coord),
                "local_predictability": np.arange(len(coord), dtype=float),
                "proportion_cpgs_with_sig_meqtl": [0.5] * len(coord),
                "n_tested_cpgs": [2] * len(coord),
            })

        burdens = {
            "caudate": burden(["1:10-20", "1:30-40"], ["1", "2"]),
            "dlpfc": burden(["1:10-20", "1:50-60"], ["99", "1"]),
            "hippocampus": burden(["1:10-20", "1:70-80"], ["5", "1"]),
        }
        rows, _ = CROSS.pairwise_vmr_sharing(burdens)
        pair = next(row for row in rows if row["contrast"] == "caudate_vs_dlpfc")
        self.assertEqual(pair["n_shared_tested"], 1)

    def test_effect_concordance_excludes_different_leads(self):
        n = 25
        common = [f"c{i}" for i in range(n)]
        different = [f"d{i}" for i in range(n)]

        def leads(region):
            ids = common + different
            variants = [f"v{i}" for i in range(n)]
            variants += ([f"a{i}" for i in range(n)] if region == "caudate" else [f"b{i}" for i in range(n)])
            slopes = np.r_[np.ones(n), np.ones(n) if region == "caudate" else -np.ones(n)]
            return pd.DataFrame({
                "phenotype_id": ids,
                "variant_id": variants,
                "slope": slopes,
                "slope_se": np.ones(2 * n),
                "qval": np.repeat(0.01, 2 * n),
                "sig": np.repeat(True, 2 * n),
                "z": slopes,
            })

        rows = CROSS.pairwise_cpg_concordance(
            {"caudate": leads("caudate"), "dlpfc": leads("dlpfc"), "hippocampus": leads("dlpfc")},
            0.05,
        )
        row = next(x for x in rows if x["contrast"] == "caudate_vs_dlpfc")
        self.assertEqual(row["n_same_lead_variant"], n)
        self.assertEqual(row["direction_concordance_both_sig"], 1.0)
        self.assertEqual(row["n_both_sig_different_lead_not_compared"], n)

    def test_external_aggregation_uses_only_assayed_cpgs(self):
        ext = pd.DataFrame({
            "vmr_id": ["v1", "v1", "v2"],
            "external_assayed": [1, 0, 1],
            "external_meqtl_support": [0, 1, 1],
            "assay_universe_complete": [True, True, True],
        })
        result = EXTERNAL.aggregate_vmr_external(ext).set_index("vmr_id")
        self.assertEqual(result.loc["v1", "n_cpgs_annotated"], 1)
        self.assertEqual(result.loc["v1", "n_cpgs_with_external_support"], 0)

    def test_complete_assay_universe_must_contain_negatives(self):
        """A resource labelled complete but carrying only positives is not estimable.

        Regression for the Phase 3 failure: the Jaffe GEO 'allPairs' file is a
        significant-results table, so every row was a positive while
        assay_universe_complete was hardcoded True. The outcome was constant and the
        GLM perfectly separated. Fitting must report that, not a coefficient.
        """
        positives_only = pd.DataFrame({
            "external_meqtl_support": [1] * 60,
            "local_predictability": np.linspace(0, 1, 60),
            "assay_universe_complete": [True] * 60,
            "n_cpgs_annotated": [3] * 60,
            "n_cpgs_with_external_support": [3] * 60,
        })
        rows = EXTERNAL.fit_support_models(positives_only)
        primary = next(r for r in rows if r["model"] == "unadjusted")
        self.assertIn("error", primary)
        self.assertNotIn("coef_predictability", primary)

        with_negatives = positives_only.copy()
        with_negatives.loc[:29, "external_meqtl_support"] = 0
        with_negatives.loc[:29, "n_cpgs_with_external_support"] = 0
        rows = EXTERNAL.fit_support_models(with_negatives)
        primary = next(r for r in rows if r["model"] == "unadjusted")
        self.assertNotIn("error", primary)
        self.assertTrue(np.isfinite(primary["coef_predictability"]))

    def test_harmonized_array_resource_carries_tested_negatives(self):
        """Every 450K-manifest probe becomes a row; non-significant probes are 0."""
        universe = pd.DataFrame({
            "probe_id": [f"cg{i:08d}" for i in range(5)],
            "chrom": ["chr1"] * 5,
            "pos_1based": [100, 200, 300, 400, 500],
        })
        sig = pd.DataFrame({
            "probe_id": ["cg00000001", "cg00000003"],
            "external_pvalue": [1e-9, 1e-8],
            "external_fdr": [1e-6, 1e-5],
            "external_beta": [0.4, -0.3],
            "lead_snp": ["rs1", "rs2"],
            "snp_chrom": ["chr1", "chr1"],
            "snp_pos": [150, 350],
        })
        original = HARMONIZE.ensure_450k_universe
        HARMONIZE.ensure_450k_universe = lambda: universe.copy()
        try:
            out, stats = HARMONIZE.build_array_resource("test_450k", "DLPFC", sig, 0.05)
        finally:
            HARMONIZE.ensure_450k_universe = original

        self.assertEqual(len(out), 5)
        self.assertEqual(int(out["external_meqtl_support"].sum()), 2)
        self.assertEqual(out["external_meqtl_support"].nunique(), 2)
        self.assertTrue(out["external_assayed"].eq(1).all())
        self.assertTrue(out["assay_universe_complete"].all())
        self.assertEqual(stats["n_probes_universe"], 5)
        self.assertEqual(stats["n_probes_supported"], 2)

    def test_non_450k_resource_is_refused_rather_than_given_negatives(self):
        """Guard against joining a non-array catalog to the array manifest."""
        universe = pd.DataFrame({
            "probe_id": [f"cg{i:08d}" for i in range(5)],
            "chrom": ["chr1"] * 5,
            "pos_1based": [100, 200, 300, 400, 500],
        })
        sig = pd.DataFrame({
            "probe_id": ["not_a_probe_1", "not_a_probe_2", "cg00000001"],
            "external_pvalue": [1e-9] * 3,
            "external_fdr": [1e-6] * 3,
            "external_beta": [0.1] * 3,
            "lead_snp": ["rs1"] * 3,
            "snp_chrom": ["chr1"] * 3,
            "snp_pos": [150] * 3,
        })
        original = HARMONIZE.ensure_450k_universe
        HARMONIZE.ensure_450k_universe = lambda: universe.copy()
        try:
            with self.assertRaises(SystemExit):
                HARMONIZE.build_array_resource("bogus", "DLPFC", sig, 0.05)
        finally:
            HARMONIZE.ensure_450k_universe = original

    def test_paired_randomization_is_deterministic(self):
        differences = np.array([1.0, 1.0, 1.0, -1.0])
        p1 = paired_randomization_pvalue(differences, seed=7, n_perm=1000)
        p2 = paired_randomization_pvalue(differences, seed=7, n_perm=1000)
        self.assertEqual(p1, p2)

    def test_downsampling_uses_identical_cpg_universe(self):
        downsample = pd.Series([0.01, 0.20, 0.01], index=["shared1", "shared2", "caudate_only"])
        dlpfc = pd.Series([0.01, 0.01, 0.01], index=["shared1", "shared2", "dlpfc_only"])
        hippocampus = pd.Series([0.20, 0.20, 0.01], index=["shared1", "shared2", "hip_only"])
        rates = DOWNSAMPLE.common_universe_rates(downsample, dlpfc, hippocampus, 0.05)
        self.assertEqual(rates["n_common_cpgs_all3"], 2)
        self.assertEqual(rates["caudate_downsample_common_discovery_rate"], 0.5)
        self.assertEqual(rates["dlpfc_common_discovery_rate"], 1.0)
        self.assertEqual(rates["hippocampus_common_discovery_rate"], 0.0)


if __name__ == "__main__":
    unittest.main()
