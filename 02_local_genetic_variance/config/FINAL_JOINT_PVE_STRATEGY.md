# Final joint local-cis PVE estimator experiment

**Decision lock:** 2026-08-20, before fitting the joint estimator or generating
the independent validation simulations.

**Terminal status:** this is the final Module 02 experiment for absolute
locus-level PVE. No further estimator family, knot, penalty, weighting, or
acceptance-threshold changes are permitted after independent validation is
opened.

## 1. Biological Claim

Adult-brain VMRs can be placed on an absolute, simulation-calibrated scale of
local cis-SNP-explained methylation variance when complementary BSLMM,
Haseman--Elston, and nested out-of-fold prediction diagnostics are combined.

## 2. Primary Hypothesis

A single constrained calibration model using raw BSLMM PVE, unclipped HE,
`rho2_oof`, `r2_oof`, sample size, SNP count, genotype effective rank, and LD
will estimate true PVE with controlled null/low-PVE behavior, acceptable
whole-grid bias and RMSE, and preserved ordering on a completely independent
expanded-grid validation set.

## 3. Statistical Unit and Design

The statistical unit is one independently seeded simulated cis locus. The
locked grid crosses n = 117, 118, 153, 173, 177, 282; p = 100, 500, 1500,
3000, 5000, 12000; AR(1) LD = 0.2, 0.5, 0.8; sparse, oligogenic, and polygenic
architectures; and true PVE = 0, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 0.95.

The completed `bslmm-validation-20260817` simulations are reclassified as
development data. Their raw BSLMM results may be joined to newly regenerated
EN, HE, effective-rank, and LD features using the identical scenario seeds.
Replicates 1--3 within every architecture-by-PVE cell fit the regression;
replicates 4--5 calibrate its simulation-reference interval and null cutoff.
Those data can never again determine acceptance.

The sole acceptance set uses seed offset 500,000,000 and 120 simulations per
n/p/LD stratum, balanced as five replicates per architecture-by-PVE cell. It
is not used for fitting, interval construction, cutoff construction, model
selection, or model-family revision.

## 4. Primary Comparison

Compare `pve_cis_joint_calibrated` with known true PVE over all 12,960
independent validation loci. The prespecified low-PVE subset is exactly true
PVE 0, 0.05, and 0.10.

## 5. Covariates and Exclusions

The phenotype covariate model and nested EN settings remain identical to the
expanded simulation design. No architecture, n, p, LD, or PVE stratum may be
excluded. A missing BSLMM, HE, EN, effective-rank, or LD feature is a
computational failure and fails acceptance. Numerical prediction is bounded to
[0,1], while the unbounded linear predictor is retained for audit.

`p_eff` is the effective rank of the standardized genotype GRM:
`trace(K)^2 / trace(K^2)`, where `K = ZZ'/p`. LD is the existing median
adjacent-SNP r-squared diagnostic. These definitions are fixed for simulated
and observed loci.

## 6. Expected Effect Direction

Conditional on n, p, `p_eff`, and LD, calibrated PVE must be monotonically
nondecreasing in each of raw BSLMM PVE, HE, `rho2_oof`, and `r2_oof`.
Validation means must increase from true PVE 0 to 0.05 to 0.10, null estimates
must remain near zero, and ordering must remain positive across the full grid.

## 7. Minimal Primary Analysis

The one allowed model family is a weighted, nonnegative, monotone hinge-basis
ridge regression with a bounded output:

1. Clip raw BSLMM PVE and `rho2_oof` to [0,1], HE to [-1,2], and `r2_oof` to
   [-1,1]. Missing or nonfinite values are failures, not imputed.
2. Expand each signal diagnostic into fixed hinge bases at the knots recorded
   in `joint-pve-20260820.tsv`.
3. Include centered/scaled additive design terms for log(n), log(p),
   log(`p_eff`), and LD.
4. Fit Gaussian ridge regression with fixed lambda 0.0001. Signal-basis
   coefficients have lower bound zero; design coefficients are unconstrained.
   There is no hyperparameter search or candidate comparison.
5. Weight true PVE 0, 0.05, and 0.10 by 8, 6, and 4, respectively; all other
   levels receive weight 1.
6. Bound the reported point estimate to [0,1] and retain the unbounded value.
7. Use the held-out development-calibration subset to construct a symmetric
   95% split-conformal simulation-reference interval and a finite-sample 95th
   percentile null cutoff. These are simulation-reference quantities, not
   confidence intervals or locus-level tests.

Every hard criterion in `joint-pve-acceptance-criteria.tsv` must pass. In this
final experiment RMSE <= 0.20 is a hard requirement, not a guardrail.

## 8. Sensitivity Analyses

After the primary decision is frozen, report bias and RMSE by true-PVE level,
architecture, n, p, and LD; boundary-hit frequency; and coefficients grouped
by BSLMM, HE, `rho2_oof`, `r2_oof`, and design terms. Prespecified diagnostic
ablations (BSLMM omitted; EN diagnostics omitted; HE omitted) may quantify
incremental information but cannot replace the primary estimator or alter the
decision.

## 9. Orthogonal Validation

The independent simulations are the only validation of absolute PVE. If the
estimator passes, later CpG cis-meQTL burden, repeat/repressive architecture,
and transcription/splicing coupling provide biological convergence for the
continuous local-genetic-control axis; they do not validate numerical PVE.

## 10. Main Figure-Worthy Result

A validation-only truth-versus-estimate panel will show individual estimates,
mean and RMSE at every true-PVE level, an enlarged 0/0.05/0.10 inset, the
identity line, interval coverage, and the prespecified pass/fail criteria.

## 11. Reviewer Objections and Responses

- **The model was invented after seeing estimator failures.** The family,
  knots, weights, penalty, split, seed block, and gates are locked here before
  joint fitting or validation generation. The prior BSLMM validation is
  explicitly demoted to development data.
- **Calibration can hide a nonidentifiable estimand.** Independent RMSE, bias,
  low-PVE, coverage, null, ordering, reconciliation, and failure gates all have
  to pass. Calibration cannot be accepted on correlation alone.
- **The methods are not independent.** They are not claimed to be independent
  studies; they contribute partially distinct diagnostics from a common locus.
  The model tests whether their joint information is sufficient.
- **Boundary clipping creates artificial success.** Unbounded predictions and
  boundary-hit rates are retained; bias and RMSE are evaluated after the
  declared [0,1] estimand boundary, including true PVE 0 and 0.95.
- **A favorable subset drove the result.** The grid is balanced and no stratum
  can be removed. All criteria are evaluated on the full untouched validation.

## 12. Implementation Notes

The training run ID is `lgv-joint-pve-train-20260820`; the validation run ID is
`lgv-joint-pve-validate-20260820`. The model checksum, configuration checksum,
source development run, seed manifests, software versions, SLURM IDs, scenario
reconciliation, and output checksums must be recorded.

The terminal rule is binary:

- **PASS:** promote `pve_cis_joint_calibrated` as *simulation-calibrated local
  cis-SNP PVE*, then adapt the frozen feature definitions to corrected observed
  VMRs.
- **FAIL:** conclude that absolute locus-level PVE is not identifiable with
  sufficient reliability at these sample sizes. End estimator development and
  move the manuscript to a continuous relative/local-genetic-control axis.

