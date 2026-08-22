# Relative local-genetic-control strategy

**Interpretation lock:** 2026-08-21, after the absolute-PVE acceptance decision
and before any observed-score or downstream biological analysis.

**Analysis intent:** We will test whether VMRs with a higher within-cell rank of
the frozen simulation-calibrated joint estimator differ in repeat/repressive
annotation, cis-meQTL burden, transcriptional coupling, and schizophrenia-risk
support, supporting the claim that local genetic variation organizes a
continuous gradient of methylation control across adult human brain VMRs.

## 1. Biological Claim

Brain VMRs vary continuously in their relative local SNP contribution. Stronger
relative local genetic control can be tested for association with repeat-rich
and repressive genomic compartments and with orthogonal regulatory evidence.

This is not a claim that the estimator provides an exact locus-level percentage
of variance explained.

## 2. Primary Hypothesis

Within each cohort-by-region analysis cell, increasing local SNP contribution
score is associated with greater H3K9me3, quiescent-chromatin, and LINE/L1
overlap and with greater CpG cis-meQTL burden. Transcription/splicing and
schizophrenia-risk coupling are prespecified downstream extensions.

## 3. Statistical Unit and Design

The statistical unit is one eligible corrected VMR. The score is constructed
separately within each cohort-by-region cell from the frozen joint estimator.
Region-specific models are primary; cross-region evidence compares directions
and locked inference rather than comparing raw score levels.

Eligibility requires complete BSLMM, HE, and nested-prediction diagnostics; no
computational failure; a finite frozen-model estimate; and an observed locus
inside the expanded simulation design domain. Excluded loci remain in the
output with explicit reasons.

## 4. Primary Comparison

The manuscript-facing endpoint is `local_snp_contribution_score`, the empirical
within-cell midrank percentile of the internal
`pve_cis_joint_calibrated` estimate. Primary models use
`local_snp_contribution_score_z`, its within-cell standardized value.

This transformation uses the estimator property supported by validation:
relative ordering. It deliberately discards unsupported absolute-PVE spacing.

## 5. Covariates and Exclusions

Downstream models retain their module-specific covariates, including VMR
geometry, CpG density, methylation properties, WGBS coverage, local SNP
opportunity, mappability, segmental duplication/problematic-region overlap,
broad genomic annotation, and cell-composition variables where applicable.

No locus is excluded based on whether its score is favorable. Boundary hits,
ties, and design-domain exclusions are reported.

## 6. Expected Effect Direction

The biology-forward primary expectation is a positive gradient: higher local
SNP contribution scores should accompany greater repeat/repressive overlap and
greater cis-meQTL burden. Regionally heterogeneous LINE/L1 effects remain
allowed under the Module 04 interpretation gate.

## 7. Minimal Primary Analysis

Fit one adjusted continuous model per cohort-by-region cell using
`local_snp_contribution_score_z`. Report effect direction, uncertainty,
multiple-testing-adjusted inference, denominators, exclusions, and the exact
source run and `vmr_set_id`.

The prespecified absolute-PVE result remains
`FAIL_PIVOT_TO_RELATIVE_LOCAL_GENETIC_CONTROL`. The manuscript-use result is
`PASS_RELATIVE_GENETIC_CONTROL / FAIL_ABSOLUTE_LOCUS_PVE`. This is an
interpretation pivot, not a retroactive change to any validation gate.

## 8. Sensitivity Analyses

- top versus bottom quartile is secondary and uses the locked 25th/75th
  percentile boundaries; neither group is called genetically controlled or
  uncontrolled;
- show score deciles for visualization only while retaining continuous-model
  inference;
- audit tied ranks and lower/upper frozen-model boundary hits;
- repeat analyses after high-mappability, SNP-proximity, segmental-duplication,
  and cell-composition restrictions;
- where useful, fit a spline in score percentile to expose nonlinearity;
- report raw-estimate distributions only as descriptive, simulation-calibrated
  contribution estimates with an explicit absolute-calibration caveat.

## 9. Orthogonal Validation

CpG cis-meQTL burden and transcription/splicing coupling are convergent
biological evidence. Held-out local SNP prediction is a separate translational
endpoint. External eQTL support in the schizophrenia application is supportive
but does not establish mediation or causality.

## 10. Main Figure-Worthy Result

Show adjusted annotation or burden trends across score deciles for visual
interpretation, accompanied by the continuous per-SD model estimate. The figure
must label the x-axis “local SNP contribution score,” not PVE or heritability.

## 11. Reviewer Objections and Responses

- **Post hoc salvage:** absolute-PVE failure remains unchanged and visible; the
  original terminal rule explicitly required a pivot to a relative axis.
- **Null inflation:** no zero-versus-positive or 0.10 threshold is used; the
  bottom quartile is only relatively lower within the eligible analysis cell.
- **Cross-region comparability:** scores are ranked within cell and region
  results are compared by direction/inference, not by raw score means.
- **Simulation transportability:** primary inference is restricted to the
  expanded design domain and boundary/extrapolation counts are reported.
- **Nonlinearity:** percentile ranks, decile displays, and spline sensitivity
  avoid treating the failed absolute scale as interval-calibrated.

## 12. Implementation Notes

The executable lock is `config/local_genetic_control.yml`. Internal outputs may
retain `pve_cis_joint_calibrated` for reproducibility, but active downstream
models must consume `local_snp_contribution_score_z`. Code must fail if an
absolute PVE threshold, `positive_signal`, `h2_unscaled`, or `r_squared_cv` is
used to define biological groups. The accepted observed Module 02 run must be
listed in the README before Modules 04--09 can run in production.
