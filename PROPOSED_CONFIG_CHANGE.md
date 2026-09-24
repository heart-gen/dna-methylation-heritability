# Proposed edits to `config/covariates.yml` (PI action required)

Branch `module05/adopt-locked-covariate-model`, 2026-09-24. `config/covariates.yml`
is PI-locked, so nothing below has been applied. AGENTS.md §12: agents may
recommend, never silently decide.

Context: the PI decided (F9) that the lock is authoritative —
`primary_meqtl.locked_model: M3a` governs, and
`05_cpg_meqtl_burden/_h/01b_prepare_meqtl_inputs.py`'s `n_pc = 3` with no
methylation PC does not. The code now implements M3a. Two things the code cannot
make correct on its own are below.

---

## Edit 1 — pin the methPC recipe (required before the rerun is citable)

**Why the code alone cannot be correct.** `latent_factor_policy` names a *method*
("PCA on M0-residualized CpG phenotypes") and no specification: no file, no CpG
set, no subsample size, no seed, no standardization, and no answer to whether the
M0 residualization carries the ancestry PCs. AGENTS.md §9 requires a production
output to carry its configuration, and a method description is not a
configuration. The code therefore resolves every one of those choices explicitly,
writes them to `results/latent-factor-provenance.tsv`, and records
`recipe_source: code_default` so the manifest does not imply the config supplied
them. Applying this edit flips that field to `recipe_source: config`, which is
what makes the covariate model a PI decision rather than an agent's.

The values are the 2026-09-23 decision pilot's resolution
(`05_cpg_meqtl_burden/PILOT_COVARIATE_MODEL.md`), i.e. the v1 implementation of the
same method (`meqtl-validation/01_cpg_meqtl_mapping/_h/09_estimate_latent_factors.py`)
re-fitted on the v2 tested-CpG phenotypes. Two entries are stricter than the
pilot and are called out below the diff.

```diff
--- a/config/covariates.yml
+++ b/config/covariates.yml
@@ -36,6 +36,42 @@ primary_meqtl:
   latent_factor_policy: >
     Locked primary uses methPC1–5 (PCA on M0-residualized CpG phenotypes).
     Do not auto-include a larger surrogate-variable set without a new lock review.
+  # Pinned 2026-09-__ . latent_factor_policy names a METHOD; this is the
+  # specification. Read by 00_shared/covariate_lock.py::latent_factor_recipe()
+  # and executed by 05_cpg_meqtl_burden/_h/01a_estimate_latent_factors.py, which
+  # records every value in results/latent-factor-provenance.tsv. The block must
+  # be complete or absent: a partial block is refused, because it would look
+  # pinned without being pinned.
+  latent_factor_recipe:
+    method: pca_on_m0_residuals
+    # The v1 implementation of this method, re-fitted on the v2 tested CpGs. The
+    # v1 factor TABLES are deliberately not reused: they were estimated on the v1
+    # VMR catalog's CpGs, which is a different CpG set.
+    implementation: "v1 09_estimate_latent_factors.py, re-fitted on v2 tested CpGs"
+    seed: 20260730
+    n_cpg_subsample: 50000
+    n_factors_estimated: 15
+    # M0 = agedeath + sex + primarydx + snpPC1-5, i.e. EVERY locked ancestry PC.
+    # sensitivity_models defines M3a as "M0 + methPC1-5"; the policy line above
+    # does not say whether the residualization carries the ancestry PCs, and this
+    # is the answer.
+    residualize_on: M0
+    cpg_standardization: zscore_per_cpg
+    # The subsample is drawn by INDEX, so the row order of the pooled CpG matrix
+    # decides which CpGs are selected. AGENTS.md 10.1 requires numeric chromosome
+    # ordering; without this key the recipe is not reproducible from its own
+    # description.
+    cpg_order: chrom_numeric_then_position
+    missing_imputation: cpg_mean
+    # An explicit solver rather than a library default: sklearn's
+    # svd_solver="auto" selects the RANDOMIZED solver at these shapes
+    # (153 x 50,000), so the factors would be reproducible only against one
+    # sklearn version and a minor-release change could silently alter a
+    # production covariate.
+    pca_solver: numpy_svd_full
+    sign_convention: max_abs_loading_positive
```

Two deviations from the pilot, both deliberate, both recorded:

- **`cpg_order`.** The pilot pooled the autosomes in
  `sorted(glob("chr*.phenotype.bed.gz"))` order — chr1, chr10, chr11, …, chr2 —
  which is not reproducible from a written description and violates AGENTS.md
  §10.1. Numeric order draws a different 50,000-CpG subsample, so the factors are
  not numerically identical to the pilot's. They agree closely, which is the
  expected behaviour of a PCA over 50,000 of 191,946 CpGs: methPC1-5 explain
  **15.82%** of residual CpG variance here (6.06/4.30/2.40/1.74/1.32%) against
  the pilot's **15.8%** (5.96/4.36/2.37/1.76/1.32%).
- **`pca_solver`.** As annotated in the diff.

## Edit 2 — `cell_composition: sensitivity_only` is now false of the primary model

**Why the code cannot fix this.** It is a statement about what the primary model
is, and under M3a it is not true. The pilot measured **methPC1 as 72% explained by
this region's RNA MuSiC cell proportions** (R² = 0.721; Oligo ρ = +0.766,
p = 9.2e-31; D1-SPN ρ = −0.63; D2-SPN ρ = −0.62). Adopting M3a therefore imports
a substantial cell-composition adjustment into the primary meQTL scan, while line
39 says cell composition is a sensitivity only and lines 49-50 register the
cell-adjusted designs as the separate, gated M5 and M6d. The PI adopted M3a
knowing this; what remains is to stop the config asserting the opposite.

It also weakens the M6d contrast by moving its baseline: M6d is literally
`M3a + dnamCellPC1-3`, so much of what those PCs would adjust for is already in
M3a. That is a note on how to read M6d, not a reason to change it.

This is a bulk-tissue correlation between a methylation PC and an RNA-derived
proportion estimate. It establishes collinearity and **not** a cell type of
origin, which AGENTS.md §2.3 forbids inferring.

```diff
--- a/config/covariates.yml
+++ b/config/covariates.yml
@@ -39,7 +39,20 @@ primary_meqtl:
-  cell_composition: sensitivity_only
+  # NOT sensitivity_only under the locked M3a. methPC1 is 72% explained by this
+  # region's RNA MuSiC cell proportions (R2 = 0.721; Oligo rho = +0.766,
+  # p = 9.2e-31), so the locked primary model carries a substantial
+  # cell-composition adjustment through its latent factors. Recorded so that no
+  # downstream text describes the primary scan as cell-composition-free. This is
+  # collinearity between a methylation PC and an RNA-derived proportion estimate,
+  # NOT a cell type of origin (AGENTS.md 2.3).
+  cell_composition: implicit_in_primary_via_latent_factors
+  cell_composition_explicit_models: sensitivity_only   # M5 and M6d below
   # Phase 1 lock sensitivity grid (see PHASE1_COVARIATE_SENSITIVITY_PLAN.md)
   sensitivity_models:
     M0: baseline (agedeath+sex+primarydx+snpPC1-5)
     M1: M0 + pmi + ph
     M2: M0 + snpPC6-10
     M3a: M0 + methPC1-5 (PCA on M0 residuals; PEER unavailable) [LOCKED PRIMARY]
     M3b: M0 + methPC1-10
     M3c: M0 + methPC1-15
     M4: M1 + best methPC k (after pilot)
     M5: M0 + RNA MuSiC cellPC1-3 (legacy sensitivity only)
-    M6d: M3a + DNAm scMD dnamCellPC1-3 (gated sensitivity only)
+    M6d: M3a + DNAm scMD dnamCellPC1-3 (gated sensitivity only; its baseline M3a
+      already carries cell-composition structure through methPC1, so the M6d
+      contrast is an INCREMENT over that, not an unadjusted-vs-adjusted contrast)
```

If the PI prefers the key to keep its current value, the alternative is to leave
line 39 alone and add `cell_composition_note:` with the same text. What must not
happen is for `sensitivity_only` to stand unqualified, because it is the sentence
a Methods paragraph would be written from.

### Not blocked on either edit

The rerun can proceed without them: `recipe_source: code_default` is honest, and
`03_vmr_burden.R` already writes the cell-composition constraint into every run's
`results/interpretation-constraints.txt`. Edit 1 is what moves the recipe from
this agent's resolution to a PI decision; Edit 2 is what stops the config
asserting something the locked model does not do.
