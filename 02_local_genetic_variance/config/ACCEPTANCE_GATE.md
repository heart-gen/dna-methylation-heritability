# Historical EN acceptance strategy — low-h² calibration and VMR ordering

This document records the superseded elastic-net calibration gate. The final
joint-estimator experiment and the 2026-08-21 interpretation decision replaced
it for active Module 02 work. The associated executables are preserved under
`../_h/archive/2026-08-21_estimator_development/`; this file is not an active
workflow specification. See `RELATIVE_LOCAL_GENETIC_CONTROL_STRATEGY.md` for
the manuscript-facing endpoint and the module README for the active `00..06`
pipeline.

## 1. Biological Claim

`h2_en_calibrated` can support continuous comparisons of local genetic control
across brain VMRs when it does not inflate null/low-signal loci and preserves
their biological ordering.

## 2. Primary Hypothesis

Across an independent simulation split, the estimator has small bias at true
h² ≤ 0.10, controlled null positive calls, acceptable level-specific bias, and
a strong monotonic relationship between true and estimated h².

## 3. Statistical Unit and Design

The unit is one independently seeded simulated locus. Calibration-only
simulations select the estimator; an untouched evaluation split determines
acceptance. Architecture, sample size, SNP count, and LD strata remain balanced.

## 4. Primary Comparison

Compare estimated with true local SNP-explained variance across the locked grid,
with true h² ≤ 0.10 prespecified as the biologically consequential low-signal
range.

## 5. Covariates and Exclusions

No post hoc strata are excluded. Numerical interpretation remains limited to
the simulated design domain, and computational failures remain disqualifying.

## 6. Expected Effect Direction

Null and low-signal estimates should be centered near truth, positive-signal
calls should be controlled under h²=0, and estimated VMR ranks should increase
with true h².

## 7. Minimal Primary Analysis

Hard gates cover null type-I error, null mean, pooled low-h² bias, maximum
low-h² level bias, maximum whole-grid level bias, interval coverage, overall
mean bias, and Spearman ordering. Aggregate RMSE ≤0.20 is a guardrail only.

## 8. Sensitivity Analyses

Report bias separately at every true-h² level and by sample-size, SNP-count,
LD, and architecture strata. Keep high-dimensional and high-h² dispersion
visible even when it does not determine acceptance.

## 9. Orthogonal Validation

Haseman–Elston remains the control variate. Downstream CpG cis-meQTL burden,
repeat/repressive enrichment, and transcription/splicing coupling test whether
the accepted VMR ordering has biological support.

## 10. Main Figure-Worthy Result

A truth-versus-estimate panel should emphasize level means and rank ordering,
with null/low-h² levels enlarged and the aggregate RMSE shown as a secondary
precision descriptor.

## 11. Reviewer Objections and Responses

- A pooled RMSE can be dominated by high-h²/high-p dispersion: it is reported
  against a 0.20 target but cannot alone veto an otherwise calibrated ranking.
- Bias can cancel across levels: pooled low-h² and per-level bias are both hard
  gates.
- A monotone but numerically noisy estimator may still support enrichment
  analyses: Spearman ≥0.70 is a hard gate, while numerical claims retain the
  calibration-domain restrictions.
- A revised gate could grandfather prior results: model and gate versions must
  match, and a fresh independent validation split is required.

## 12. Implementation Notes

The historical thresholds and roles are in `acceptance-criteria.tsv`;
`analysis.tsv` locks `low_h2_max=0.10`. Archived
`03_evaluate_calibration.R` wrote the low-h² and by-level metrics, and archived
`06_check_acceptance.R` failed only on hard criteria while recording guardrail
misses. Calibration fitters stamped
`acceptance_gate_version` into the model so historical models cannot pass the
new gate silently.
