# Methods summary (generated)

Run `fig-all-20261009-b`. Parameters are read from the locked configuration snapshotted into this run (`code/config/`); run IDs from the module acceptance tables. Do not edit.

## Cohort and catalog
- Primary donor group `AA`; locked donor counts: Caudate n = 153, DLPFC n = 118, Hippocampus n = 117.
- VMR sets: `vmrset-AA-caudate-937a41979978`, `vmrset-AA-dlpfc-856067dfe289`, `vmrset-AA-hippocampus-2d907b892215` [tables/exclusions-and-denominators.tsv].

## Local genetic control (Module 02)
- cis window +/-500,000 bp; SNP MAF >= 0.05, missingness <= 0.05; >= 100 cis variants.
- Endpoint: `local_snp_contribution_score`, within_cell_empirical_midrank_percentile within cohort_by_region; model predictor `local_snp_contribution_score_z`. Frozen model SHA-256 `9f26c3273746fda85d9bbf21e224857db9a1ad79a521582a12f241854c03223a`; validation over 12,960 simulations: relative ordering pass, absolute PVE fail (`PASS_RELATIVE_GENETIC_CONTROL_FAIL_ABSOLUTE_LOCUS_PVE`).
- cis-GREML sensitivity (Module 02c): the same window, SNP QC and covariates; one GRM per VMR; GCTA 1.94.1 modes `fusion_unconstrained` (primary: --reml-no-constrain --reml-lrt 1); `fisher_unconstrained` (convergence_sensitivity: --reml-no-constrain --reml-alg 1 --reml-lrt 1). Existence = mean estimate over converged eligible VMRs, ordering = Spearman with `local_snp_contribution_score`, each with a delete-one-chromosome jackknife CI. A GCTA error after a runaway iteration trace, or a fit with h2 SE 0 or undefined, is an estimator outcome, not a converged fit.

## Held-out prediction (Module 03)
- `end_to_end_oof`: 5 outer x 5 inner folds, 5 repeats; elastic-net alpha grid {0.1, 0.5, 0.9, 1}, `lambda.min`; fold-internal screen `he_permutation` (1000 permutations, alpha 0.05); failed screens get `training_mean`.

## Axis inference (Modules 09b, 10)
- Debiased outcomes; variance = donor bootstrap + delete-one-chromosome block jackknife; relative magnitudes gated on Fieller denominator stability (`relative_magnitude_gate`), which gates the magnitude, never the test.

## Design constraint
- Caudate is sequencing batch 3 and DLPFC/hippocampus are batches 1-2, so region and batch are perfectly confounded. Single-region results are not exposed; caudate-vs-other differences are descriptive.

## Accepted upstream runs
- `01_vmr_catalog`: `vmrcat-AA-caudate-20260816`, `vmrcat-AA-dlpfc-20260816`, `vmrcat-AA-hippocampus-20260816`, `vmrcat-all_individuals-caudate-20260816`, `vmrcat-all_individuals-dlpfc-20260816`, `vmrcat-all_individuals-hippocampus-20260816` -- `ALL_FIVE_CRITERIA_PASS`
- `01b_estimation_cells`: `estcell-all_individuals.AA-caudate-20260910`, `estcell-all_individuals.AA-dlpfc-20260910`, `estcell-all_individuals.AA-hippocampus-20260910`, `estcell-all_individuals.EA-caudate-20260910`, `estcell-all_individuals.EA-dlpfc-20260910`, `estcell-all_individuals.EA-hippocampus-20260910`, `estcell-AA.n118r1-caudate-20260918`, `estcell-AA.n118r2-caudate-20260918`, `estcell-AA.n118r3-caudate-20260918` -- `ALL_FIVE_CRITERIA_PASS`
- `02_local_genetic_variance`: `lgv-AA-caudate-rescore-20260913`, `lgv-AA-dlpfc-rescore-20260913`, `lgv-AA-hippocampus-rescore-20260913`, `lgv-all_individuals-caudate-20260823`, `lgv-all_individuals-dlpfc-20260823`, `lgv-all_individuals-hippocampus-20260823`, `lgv-all_individuals.AA-caudate-20260913`, `lgv-all_individuals.AA-dlpfc-20260913`, `lgv-all_individuals.AA-hippocampus-20260913`, `lgv-all_individuals.EA-caudate-20260917`, `lgv-all_individuals.EA-dlpfc-20260917`, `lgv-all_individuals.EA-hippocampus-20260917`, `lgv-AA.n118r1-caudate-20260918`, `lgv-AA.n118r2-caudate-20260918`, `lgv-AA.n118r3-caudate-20260918` -- `PASS_RELATIVE_SCORE_OBSERVED_QC`
- `02b_greml_simulation_benchmark`: `greml-sim-ar1-20261003`, `greml-AA-caudate-20261003-a`, `greml-AA-dlpfc-20261003-a`, `greml-AA-hippocampus-20261003-a` -- `PASS_GREML_BENCHMARK_QC`
- `02c_cis_greml_sensitivity`: `cgs-AA-caudate-20261006-a`, `cgs-AA-dlpfc-20261006-a`, `cgs-AA-hippocampus-20261006-a` -- `PASS_CIS_GREML_SENSITIVITY_QC`
- `03_local_snp_prediction`: `lsp-AA-caudate-20260925-a`, `lsp-AA-dlpfc-20260925-a`, `lsp-AA-hippocampus-20260925-a`, `lsp-all_individuals.AA-caudate-20260918`, `lsp-all_individuals.AA-dlpfc-20260918`, `lsp-all_individuals.AA-hippocampus-20260918`, `lsp-all_individuals.EA-caudate-20260918`, `lsp-all_individuals.EA-dlpfc-20260918`, `lsp-all_individuals.EA-hippocampus-20260918`, `lsp-AA.n118r1-caudate-20260918`, `lsp-AA.n118r2-caudate-20260918`, `lsp-AA.n118r3-caudate-20260918` -- `PASS_OOF_PREDICTION_QC`
- `04_repeat_repressive_architecture`: `rra-AA-caudate-20260925-a`, `rra-AA-dlpfc-20260925-a`, `rra-AA-hippocampus-20260925-a` -- `GATES_APPLIED_2_OF_3_OUTCOMES_SUPPORTED`
  - non-gating: `rra-AA-crossregion-20261007` (Accepted non-gating sensitivity runs), `rra-AA-crossregion-20261007-a` (Accepted non-gating sensitivity runs)
- `05_cpg_meqtl_burden`: `cmb-AA-caudate-20260924`, `cmb-AA-dlpfc-20260924`, `cmb-AA-hippocampus-20260924`, `cmb-AA-crossregion-20261007` -- `PASS_CPG_MEQTL_BURDEN_QC / SLOPE_DIFFERENCE_CI_INCLUDES_ZERO`
- `06_partitioned_heritability`: `sldsc-AA-caudate-20261008`, `sldsc-AA-dlpfc-20261008`, `sldsc-AA-hippocampus-20261008` -- `PASS_PARTITIONED_H2_QC`
  - non-gating: `sldsc-AA-external-20261008` (Accepted non-gating positive control)
- `07_transcription_splicing_coupling`: `tsc-AA-caudate-20260925-b`, `tsc-AA-dlpfc-20260925-b`, `tsc-AA-hippocampus-20260925-b` -- `PASS_TX_COUPLING_QC`
- `08_region_donor_generalization`: `rdg-AA-crossregion-20261008` -- `PASS_REGION_DONOR_GENERALIZATION_QC`
- `09_schizophrenia_risk_application`: `scz-AA-caudate-20261008`, `scz-AA-dlpfc-20261008`, `scz-AA-hippocampus-20261008` -- `PASS_SCZ_APPLICATION_QC`
  - non-gating: `scz-all_individuals.EA-crossregion-20261008` (Accepted non-gating check)
- `09b_aging_application`: `age-AA-caudate-20261003`, `age-AA-dlpfc-20261003`, `age-AA-hippocampus-20261003` -- `PASS_AGING_AXIS_COVERAGE`
- `10_environmental_exploratory`: `env-AA-caudate-20261003`, `env-AA-dlpfc-20261003`, `env-AA-hippocampus-20261003` -- `PASS_EXPLORATORY_COVERAGE`
