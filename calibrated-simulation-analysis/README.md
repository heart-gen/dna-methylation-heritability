# Simulation-calibrated elastic-net local SNP variance

This directory implements a leakage-free, simulation-calibrated elastic-net
analysis for local SNP-explained methylation variance. It does not modify or
replace the existing results under `heritability/elastic_net_model/`.

The primary reported term is **simulation-calibrated elastic-net estimate of
local SNP-explained methylation variance** (`h2_en_calibrated`). The result
should only be shortened to *local SNP heritability* when the calibration
acceptance criteria pass, and it must remain qualified as elastic-net-based and
simulation-calibrated.

## Directory structure

The workflow follows the repository convention:

```text
calibrated-simulation-analysis/
├── _h/                         # R analysis and SLURM helper scripts
├── _m/                         # Generated manifests, models, results, logs, figures
├── config/                     # Prespecified grids and acceptance criteria
├── tests/                      # Deterministic unit and end-to-end smoke tests
├── environment.yml             # Pinned major R dependency versions
└── README.md
```

Generated `_m` results are excluded by the analysis-local `.gitignore`; code,
configuration, tests, and `_m/README.md` remain versionable.

## Estimand and safeguards

For each simulated or observed VMR, the analysis:

1. creates repeated outer donor folds;
2. fits the covariate adjustment using outer-training donors only;
3. imputes and filters SNPs using outer-training genotypes only;
4. performs optional phenotype-based SNP screening inside the outer training
   fold only;
5. tunes elastic-net alpha and lambda by inner cross-validation;
6. predicts outer held-out donors;
7. calculates held-out `r2_oof`, `rho2_oof`, `covariance_ratio_oof`, and
   `score_variance_ratio_oof`;
8. standardizes five cross-fit diagnostics within each N/SNP-count/LD stratum;
9. fits a pooled forward regression across design strata and removes its
   level-specific attenuation using the known simulated h2 grid;
10. blends that forward estimate with an unclipped Haseman–Elston
    method-of-moments control variate; and
11. constrains only the upper point estimate to the maximum simulated h2 (0.6),
    while retaining negative estimates to avoid reintroducing null bias.

The same loci receive an unclipped Haseman–Elston estimate (`he_h2`, with
standard error and p-value). The final hybrid uses it as a 25% unbiased control
variate and gives the forward elastic-net calibration 75% weight. That weight is
selected without using evaluation simulations: odd calibration replicates fit
the forward model, even replicates tune the weight, and the tuning null mean's
one-sided 95% upper confidence bound must satisfy the locked null criterion.
The final components are then refit on all calibration simulations.

The observed-data adapter refuses estimator-setting overrides by default,
because a different fold count, repeat count, alpha grid, lambda rule, or
screening limit requires recalibration.

The columns `h2_calibration_lower` and `h2_calibration_upper` are empirical
simulation-reference limits. They are not classical frequentist confidence
intervals and depend on the simulated architecture mixture.

`h2_en_calibrated_unbounded` retains the raw hybrid value for audit, while
`h2_en_calibrated` applies the upper calibration-domain limit. The boolean
`h2_upper_boundary_hit` identifies loci affected by that limit. A boundary hit
does not justify extrapolating beyond 0.6; a wider h2 grid would be required.

## Simulation design

[`config/analysis.tsv`](config/analysis.tsv) contains the locked primary grid.
It targets the AA discovery sample sizes and spans local SNP count, LD, sparse
through polygenic architectures, and true conditional local variance from zero
to 0.6. Sparse, oligogenic, and polygenic simulations are pooled when fitting
the primary calibration curve, preventing an unobservable causal architecture
from being selected after inspecting a real VMR.

Calibration and evaluation replicates receive different deterministic seeds.
The held-out evaluation simulations are never used to fit the forward model,
choose its control-variate weight, or set null thresholds. Null positive-signal
thresholds use a finite-sample
split-conformal upper order statistic at alpha 0.05 rather than an interpolated
sample quantile; with 30 calibration nulls per stratum, the attainable marginal
rate is conservatively 1/31 (3.23%).

Observed loci are assigned to the closest simulated sample-size/SNP-count/LD
stratum. `calibration_status` flags raw-statistic extrapolation or a locus outside
the design domain. Such estimates should not be interpreted numerically until
the simulation grid is expanded and rerun.

## Environment

Create and verify the dedicated environment once:

```bash
bash calibrated-simulation-analysis/_h/setup_environment.sh
```

All production and test scripts default to the dedicated
`/projects/p32505/opt/envs/calibrated-local-h2` environment. They do not fall
back to the broader epigenomics environment. `CAL_H2_ENV` is only an explicit
override, and submission stops before creating jobs if the selected environment
does not contain an executable `Rscript`. Every calibration and evaluation run
records R session information.

## Tests

Run the deterministic function test and small end-to-end workflow:

```bash
bash calibrated-simulation-analysis/tests/run_smoke.sh
bash calibrated-simulation-analysis/tests/test_slurm_spool_paths.sh
```

The smoke grid verifies execution only. It is too small to establish scientific
calibration. The spool-path regression test executes a copied SLURM stage script
from a temporary directory and verifies that it still finds the run-scoped R
code rather than looking beneath `/var/spool/slurmd/`.

## Full SLURM calibration workflow

From the repository root:

```bash
export SBATCH_ACCOUNT=p32505
export MAX_CONCURRENT=200
bash calibrated-simulation-analysis/_h/submit_simulation_workflow.sh ajhg-calibration-v1
```

The submission helper:

1. creates an immutable `_m/runs/{RUN_ID}/` directory and refuses to mix with
   an existing run;
2. snapshots the analysis configuration, acceptance criteria, environment
   specification, Git commit, exact Conda package specification, and SHA-256
   checksums of every submitted analysis script; the copied scripts are run
   from `_m/runs/{RUN_ID}/code/_h/`, so later repository edits cannot change an
   active run and SLURM spool paths cannot break relative script discovery;
3. generates `_m/runs/{RUN_ID}/config/scenarios.tsv`;
4. submits one simulation and nested-cross-fit array task per manifest row;
5. fits the forward hybrid calibration only after every array task succeeds;
6. evaluates and plots only after calibration succeeds, with failure of any
   locked scientific acceptance criterion returning a failed final SLURM job.

Important generated outputs are:

- `_m/runs/{RUN_ID}/calibration/elastic-net-calibration.rds`
- `_m/runs/{RUN_ID}/calibration/calibration-manifest.tsv`
- `_m/runs/{RUN_ID}/evaluation/calibrated-evaluation-estimates.tsv`
- `_m/runs/{RUN_ID}/evaluation/calibration-performance-overall.tsv`
- `_m/runs/{RUN_ID}/evaluation/calibration-performance-by-design.tsv`
- `_m/runs/{RUN_ID}/evaluation/acceptance-results.tsv`
- `_m/runs/{RUN_ID}/figures/calibration-truth-versus-estimate.pdf`
- `_m/runs/{RUN_ID}/provenance/submitted-jobs.tsv`

Prepare and validate a submission without calling `sbatch` with:

```bash
SUBMIT_CAL_H2_DRY_RUN=TRUE \
bash calibrated-simulation-analysis/_h/submit_simulation_workflow.sh dry-run-v1
```

The primary estimator is acceptable only if all prespecified thresholds in
[`config/acceptance-criteria.tsv`](config/acceptance-criteria.tsv) pass. Failure
means the manuscript should retain **local genetic predictability** and should
not report the calibrated result as local heritability.

Acceptance includes conditional calibration checks, not only performance
averaged across the simulation prior: mean estimated h2 at true h2=0 must be
≤0.05, and absolute mean bias at every simulated h2 level must be ≤0.10. These
guards prevent a regression-to-the-simulation-mean mapping from appearing
calibrated solely because positive and negative biases cancel across the grid.

If calibration/evaluation raw simulations completed but a later calibration
implementation is corrected, preserve the original run and derive new outputs
without resimulating:

```bash
bash calibrated-simulation-analysis/_h/reanalyze_completed_run.sh \
  SOURCE_RUN_ID DERIVED_RUN_ID
```

The derived run checks that all source scenarios have outputs, symlinks the
immutable raw simulations, snapshots the revised code and environment, and
records the derivation reason without overwriting the source run.

After method development, validate a frozen model on a completely new seed
split without refitting it:

```bash
bash calibrated-simulation-analysis/_h/submit_fresh_validation_workflow.sh \
  MODEL_RUN_ID VALIDATION_RUN_ID SEED_OFFSET
```

The validation workflow generates only fresh evaluation simulations, copies
the frozen model and its SHA-256 into a new immutable run, and applies the same
fail-closed seven-criterion gate. The accepted v4 validation used seed offset
`200000000` and passed all seven criteria across 2,430 new simulations.

## Applying the accepted calibration to observed VMRs

Do not submit observed jobs until the full simulation results pass every
aggregate and conditional calibration criterion. The run-scoped submission
entry point applies the accepted model to AA VMRs in all three regions:

```bash
bash calibrated-simulation-analysis/_h/submit_observed_calibrated_workflow.sh \
  CALIBRATION_RUN_ID OBSERVED_RUN_ID AA
```

The submission refuses to create an observed run if any current acceptance
criterion is missing or failed. It snapshots the accepted calibration model and
code, submits separate caudate/DLPFC/hippocampus arrays, records expected versus
analyzed/excluded/QC-failed/computationally-failed task counts, and combines
results only after every region array completes. Terminal task categories are
mutually exclusive:

- `summary/`: an estimate was produced;
- `excluded/`: a prespecified analysis exclusion, such as a sex-chromosome VMR;
- `qc_failures/`: either no SNP exists in the exact ±500 kb window or fewer
  than two SNPs remain after MAF and missingness QC;
- `failures/`: a computational or missing-input failure that requires recovery.

The configured failure-rate tolerance applies only to `qc_failures/`.
Computational failures have zero tolerance and cannot be accepted by the final
QC job, even when their numerical rate is small.

The adapter reads the existing region-specific `vmr.bed`, local PLINK files,
methylation phenotypes, and region-specific covariate files. Even when an
upstream PLINK file was extracted with a wider window, it filters the attached
map to the prespecified ±500 kb cis window, MAF ≥0.05, and SNP missingness ≤5%
before fitting. It never uses covariates from another region.

Autosomal VMRs are analyzed by default; X and Y VMRs are recorded under an
`excluded/` directory unless a sex-chromosome-specific calibration is developed.
Per-VMR output is written below the observed run directory:

```text
_m/observed-runs/{OBSERVED_RUN_ID}/results/{region}/AA/summary/
_m/observed-runs/{OBSERVED_RUN_ID}/results/{region}/AA/excluded/
_m/observed-runs/{OBSERVED_RUN_ID}/results/{region}/AA/qc_failures/
_m/observed-runs/{OBSERVED_RUN_ID}/results/{region}/AA/failures/
_m/observed-runs/{OBSERVED_RUN_ID}/results/combined/
```

Recover a failed immutable observed run into a new run directory with:

```bash
bash calibrated-simulation-analysis/_h/recover_observed_workflow.sh \
  SOURCE_OBSERVED_RUN_ID RECOVERY_OBSERVED_RUN_ID
```

Recovery copies successful terminal records, retries only computational
failures, and re-evaluates legacy zero-variant exclusions under the current QC
policy. If a per-VMR PLINK input is absent because its upstream extraction did
not complete, the retry reconstructs only the exact ±500 kb window from the
recorded cohort PGEN into the new run's `recovered-inputs/` directory. It does
not alter the source run or the shared genotype inputs. The recovery provenance
records the PGEN source, phenotype fallback, task manifests, code snapshot,
software environment, checksums, submitted jobs, and classification policy.

The automatically generated combined table retains raw prediction metrics,
Haseman–Elston estimates, forward and unbounded hybrid estimates, the
upper-domain-constrained calibrated estimate,
empirical calibration limits, null-threshold call, nearest calibration stratum,
and calibration-domain status. It does not overwrite legacy elastic-net files.

## Interpretation limits

- GREML/REML failure at small sample sizes motivates this estimator but does
  not remove the limited information in the data.
- Elastic net imposes shrinkage and sparsity assumptions; calibration is only
  transportable to architectures represented in the simulation grid.
- A positive-signal call is based on a stratum-specific simulated null
  threshold, not a conventional variance-component likelihood-ratio test.
- Numerical estimates outside the calibration domain must be withheld.
- CpG cis-meQTL burden and external brain meQTL support remain necessary
  orthogonal biological validation.
