# 02_local_genetic_variance — primary quantitative endpoint

Produces `h2_en_calibrated`:

> simulation-calibrated elastic-net estimate of local SNP-explained methylation
> variance

This is the manuscript's primary quantitative endpoint (AGENTS.md §7.2). Local
genetic control is treated as **continuous**; any quantile or extreme-group
presentation is secondary and must use a prespecified threshold.

**Status: migrated, adapter re-pointed, no accepted production run yet.**

## Relationship to `calibrated-simulation-analysis/`

This module *is* `calibrated-simulation-analysis/`, moved under the v2 numbering
and re-pointed at corrected VMRs. The estimator is unchanged and deliberately so.

What was **not** touched: the estimator, alpha grid, folds, repeats, lambda rule,
screen, and raw metric. Changing any of them requires recalibrating the entire
simulation grid (AGENTS.md §7.2), which is out of scope here.

What changed:

1. **The observed adapter reads a v2 VMR catalog.**
   `04_estimate_observed_vmr.R` hard-coded
   `vmr-analysis/all_individuals/<region>/_m/{vmr,plink_format}` and fell back to
   a second copy in the lead author's tree. It now accepts `--vmr_run_dir`
   pointing at an accepted `01_vmr_catalog` run, and when that is set there is
   **no** cross-repo fallback — a v2 run reads its own VMR set or fails. It also
   verifies that the run's cohort and region match what was requested, and
   carries `upstream_vmr_run_id` and `vmr_set_id` onto every output row.
2. **The calibration model is checksum-pinned.** `04_estimate_observed_vmr.R`
   asserts SHA-256
   `bbe9f9f3e897b19c536078c20e6bd50a2f5ea385ab1c1258039974ced855e389`
   before use — the model from `ajhg-calibration-v4-independent-validation`,
   frozen at `_m/calibration_frozen/elastic-net-calibration.rds`.
3. **V12.** `step_5_estimate_observed_vmr.sh:96` used
   `module load plink/2.0-alpha-3.3`; it now uses
   `/projects/p32505/opt/bin/plink2`.

## The estimator

Nested cross-fitting produces held-out diagnostics (`r2_oof`, `rho2_oof`,
`covariance_ratio_oof`, `score_variance_ratio_oof`), standardized within
N/SNP-count/LD strata. A pooled forward regression removes level-specific
attenuation against the known simulated h2 grid, and the result is blended 75/25
with an unclipped Haseman–Elston control variate.

- `h2_en_calibrated_unbounded` — raw hybrid, retained for audit
- `h2_en_calibrated` — clipped only at the upper end (max simulated h2 = 0.6),
  flagged by `h2_upper_boundary_hit`
- Negative values are **retained deliberately**

Calibration limits are *simulation-reference limits*, not confidence intervals,
and only within-domain loci may be interpreted numerically.

## Acceptance criteria

Fail-closed: observed jobs are refused unless every criterion is present and
passing. From `config/acceptance-criteria.tsv`, as met by
`ajhg-calibration-v4-independent-validation`:

| Metric | Rule | Threshold | Observed |
|---|---|---|---|
| absolute_mean_bias | ≤ | 0.05 | 0.0167 |
| rmse | ≤ | 0.15 | 0.1341 |
| calibration_interval_coverage | ≥ | 0.90 | 0.9449 |
| null_type1_error | ≤ | 0.06 | 0.0198 |
| spearman_truth_estimate | ≥ | 0.50 | 0.7351 |
| null_mean_estimated_h2 | ≤ | 0.05 | 0.0452 |
| max_absolute_h2_level_bias | ≤ | 0.10 | 0.0789 |

Additional v2 gate: the fraction of loci falling outside the calibration domain
must stay under `config/thresholds.yml` `gates.max_outside_calibration_domain`
(0.10). Too many out-of-domain loci is a §14 stop condition.

## Running

Both arms, with AA primary:

```bash
cd 02_local_genetic_variance/_m
CAL_H2_VMR_RUN_DIR=../../01_vmr_catalog/_m/runs/<ACCEPTED_RUN_ID> \
CAL_H2_COHORT=AA REGION=caudate POPULATION=AA \
  sbatch ../_h/step_5_estimate_observed_vmr.sh
```

`01_vmr_catalog` must have a **passing acceptance gate and a recorded run ID**
first (AGENTS.md §6). The adapter will not synthesize a VMR set.

## Terminology

Call it "local SNP-explained methylation variance". Call it "local SNP
heritability" only if the acceptance criteria pass. Do not describe the
calibration interval as a confidence interval.

## Accepted runs

| Run ID | Arm | Region | Upstream `vmr_set_id` | Acceptance | Notes |
|---|---|---|---|---|---|
| _(none)_ | | | | | |

### Verification not yet done

- No observed run against a v2 catalog — `01_vmr_catalog` has no accepted run.
- The acceptance-criteria replay and the ~50-VMR smoke described in the plan are
  blocked on that.
- The annotation sensitivity must be rerun on corrected VMRs before any
  §7.4 claim.

## `_m/` contents

`_m/runs/` and `_m/observed/` are gitignored and stay on Quest.

`_m/calibration_frozen/` holds the accepted calibration model plus its
acceptance results and validation metadata. The 4 MB
`elastic-net-calibration.rds` is **not** tracked — it stays on Quest like every
other data artifact — but `acceptance-results.tsv` and
`validation-metadata.tsv` are, so the expected checksum and the seven acceptance
numbers are readable from a clone. The authoritative copy is
`calibrated-simulation-analysis/_m/runs/ajhg-calibration-v4-independent-validation/calibration/elastic-net-calibration.rds`;
restore from there if the frozen copy is lost, and the SHA-256 assertion in
`04_estimate_observed_vmr.R` will confirm it is the right one.
