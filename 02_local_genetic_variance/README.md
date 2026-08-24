# 02_local_genetic_variance — relative local genetic control

Module 02 supplies the continuous local-genetic-control axis for downstream
biology.

**Decision locked 2026-08-21:**
`PASS_RELATIVE_GENETIC_CONTROL / FAIL_ABSOLUTE_LOCUS_PVE`.

The final zero-overlap validation contained 12,960/12,960 complete simulations,
zero computational failures, and zero seed overlap against 22 prior
seed-bearing manifests. The frozen joint model preserved strong ordering
(Spearman 0.796) but failed 6/14 absolute-PVE gates, including bias, RMSE, null
mean, and level-specific calibration. Estimator development has ended.

## Active endpoint

The internal `pve_cis_joint_unbounded` estimate is converted within each
cohort-by-region cell to:

- `local_snp_contribution_score`: empirical midrank percentile,
  `(midrank - 0.5) / n_eligible`;
- `local_snp_contribution_score_z`: standardized score used in primary models;
- `local_snp_contribution_quartile`: secondary top/bottom-quartile contrast;
- `local_snp_contribution_score_basis`: the column actually ranked.

Ranking uses the pre-clip estimate. `pve_cis_joint_calibrated =
max(0, unbounded)` is retained as the bounded descriptive output but is not
the ranking basis: clipping collapsed 62.7% of loci to one tied value, and the
2026-08-22 observed-regime grid showed that ordering is recoverable rather
than noise. See the grid section below.

The lock is [../config/local_genetic_control.yml](../config/local_genetic_control.yml)
and the scientific strategy is
[config/RELATIVE_LOCAL_GENETIC_CONTROL_STRATEGY.md](config/RELATIVE_LOCAL_GENETIC_CONTROL_STRATEGY.md).

Allowed uses:

- adjusted continuous biological-association models;
- effect direction and relative gradients;
- score deciles for visualization;
- top-versus-bottom quartiles as secondary relative contrasts;
- raw-estimate distributions only with an absolute-calibration caveat.

Prohibited uses:

- exact locus-level PVE percentages;
- absolute PVE thresholds, including 0.10;
- heritable/nonheritable or genetically/nongenetically controlled groups;
- `positive_signal` as a biological class;
- exact percentage-point annotation contrasts;
- raw score-level comparisons across regions.

## Frozen model and terminal evidence

- training run: `lgv-joint-pve-train-20260820`;
- independent validation: `lgv-joint-pve-validate-20260821a`;
- terminal decision: `lgv-joint-pve-decision-20260821a`;
- model SHA-256:
  `9f26c3273746fda85d9bbf21e224857db9a1ad79a521582a12f241854c03223a`.

The initially evaluated `lgv-joint-pve-validate-20260820` is audit-only because
384 seeds overlapped an earlier settings screen. The replacement validation
changed only its seed block and run ID; the model, design, and gates were
unchanged.

Absolute-PVE failures were:

| Metric | Result | Required |
|---|---:|---:|
| Absolute mean bias | 0.0869 | <=0.05 |
| RMSE | 0.2127 | <=0.20 |
| Null mean estimated PVE | 0.0767 | <=0.05 |
| Low-PVE absolute mean bias | 0.0347 | <=0.025 |
| Maximum low-PVE level bias | 0.0767 | <=0.05 |
| Maximum PVE-level bias | 0.2049 | <=0.15 |

Passing properties included Spearman ordering 0.796, null type-I error 0.0537,
and simulation-reference coverage 0.9545. These authorize the relative score,
not absolute PVE.

## Active `00..06` observed-score pipeline

Only the following numbered entry points are active:

| Stage | Script | Contract |
|---:|---|---|
| 00 | `_h/00_prepare_observed_run.R` | Verify the accepted Module 01 run, frozen model checksum, settings, task universe, and run identity. |
| 01 | `_h/01_estimate_observed_joint_features.R` | For one VMR, calculate raw BSLMM, HE, nested `rho2_oof`/`r2_oof`, n, p, effective rank, and LD. Every task emits one terminal row. |
| 02 | `_h/02_combine_observed_joint_features.R` | Reconcile every expected, completed, excluded, QC-failed, missing, and computationally failed task. |
| 03 | `_h/03_apply_frozen_joint_model.R` | Check the model SHA, apply it once, and label expanded-design eligibility. |
| 04 | `_h/04_derive_local_snp_contribution_score.R` | Derive the within-cell midrank score, z-score, and secondary quartile labels. |
| 05 | `_h/05_check_observed_score.R` | Enforce reconciliation, zero computational failures, >=90% within-domain coverage among complete features, score integrity, and the absolute-PVE prohibition. |
| 06 | `_h/06_finalize_observed_run.R` | Record the decision, session, output checksums, and seal the immutable run. |

Shared code is limited to `_h/00_functions.R`,
`_h/bslmm_pilot_functions.R`, and `_h/joint_pve_functions.R`. The latter two
retain historical names because their fitted-model interface is checksum-bound;
renaming or rewriting the numerical implementation is unnecessary risk.

### Submission scripts

Each stage has its own launcher, following the Module 01 convention: one
`step_*.sh` per stage carrying its own `#SBATCH` resource block, and one
workflow script that chains them with SLURM dependencies.

| Stage | Step script | Resources |
|---:|---|---|
| 01 | `_h/step_01_observed_joint_features.sh` | array, 1 cpu, 8G, 4 h |
| 02 | `_h/step_02_combine_features.sh` | 1 cpu, 12G, 2 h |
| 03 | `_h/step_03_apply_joint_model.sh` | 1 cpu, 12G, 2 h |
| 04 | `_h/step_04_derive_score.sh` | 1 cpu, 8G, 1 h |
| 05 | `_h/step_05_check_score.sh` | 1 cpu, 8G, 1 h |
| 06 | `_h/step_06_finalize_run.sh` | 1 cpu, 8G, 1 h |

Every step script reads `LGV_RUN_DIR`, `LGV_H_DIR`, and `CAL_H2_ENV` from the
environment and takes no positional arguments, so any single stage can be
resubmitted by hand against an existing run without editing the launcher:

```bash
LGV_RUN_DIR=<run dir> LGV_H_DIR=<repo>/02_local_genetic_variance/_h \
  sbatch _h/step_04_derive_score.sh
```

The workflow script is `_h/submit_observed_local_control.sh`. It runs Stage 00
on the submit host -- the chunk count Stage 00 writes is what sizes the Stage 01
array -- then submits step 01 as a throttled array and chains steps 02--06.
Step 02 depends with `afterany` so scheduler cancellations are reconciled rather
than silently blocking the audit; steps 03--06 use `afterok`. The submitted job
IDs and their step scripts are recorded in `{RUN_DIR}/submitted-jobs.tsv`.

Example smoke preparation:

```bash
cd 02_local_genetic_variance
DRY_RUN=TRUE LGV_SMOKE_N=2 LGV_VMRS_PER_CHUNK=1 \
  _h/submit_observed_local_control.sh \
  lgv-AA-caudate-20260821a AA caudate vmrcat-AA-caudate-20260816
```

Example production submission:

```bash
cd 02_local_genetic_variance
_h/submit_observed_local_control.sh \
  lgv-AA-caudate-20260821 AA caudate vmrcat-AA-caudate-20260816
```

Do not submit production until the intended run IDs and available Quest
allocation are reviewed. A completed SLURM job is not acceptance; Stage 05 must
pass and the accepted run must be entered below.

## Observed-regime diagnostic grid (`07..09`)

The frozen joint model was calibrated on AR(1)-simulated genotypes. Their
effective dimension is far higher than any real cis-window at n = 153:
simulated `p_eff` runs 24.3--274.4 (median 116.4) while observed `p_eff` runs
3.5--81.2 (median 37.2). The AR(1) grid therefore never covered the regime the
production run actually occupies, and `joint_pve_domain_status` does not test
`p_eff`, so it reported `within_domain` for all 11,239 eligible loci while
56.4% of them produced an unbounded estimate below the global minimum
(-0.0584) of all 12,960 validation simulations.

`07..09` close that gap without touching the estimator. Each scenario keeps a
real VMR's genotype, covariates and donor alignment and replaces only the
methylation phenotype with one simulated at a known true PVE, so `n`,
`num_snps`, LD and `p_eff` are the observed values by construction.
`_h/observed_locus_io.R` is the single locus reader shared with Stage 01, so
the two paths cannot drift.

| Script | Role |
| --- | --- |
| `07_make_observed_regime_manifest.R` | Opens the run; samples 192 real loci across 4 x 4 `num_snps` x `p_eff` quartile strata; crosses them with the frozen 8-PVE x 3-architecture grid x 2 replicates (9,216 scenarios) |
| `08_run_observed_regime_scenario.R` | Computes the four raw joint features per scenario chunk, caching the locus genotype |
| `09_summarize_observed_regime.R` | Reconciles, applies the frozen model unmodified, and writes the boundary-rate, rank-information and verdict tables |

Launch with `_h/submit_observed_regime_grid.sh`; config is
`config/observed-regime-20260822.tsv` (locked 2026-08-22,
`model_changed=FALSE`, `criteria_changed=FALSE`, fresh seed block 900000000).

The grid answers two questions the AR(1) grid could not:

1. Is the observed 62.7% lower-boundary rate reproduced on real cis-window
   genotypes at the observed sample size? In simulation it is not: the maximum
   over any n x PVE x architecture cell is 8.9%.
2. Below the boundary, does `pve_cis_joint_unbounded` still order loci by true
   PVE? On the AR(1) grid it does not (Spearman -0.0025, bootstrap 95% CI
   [-0.133, 0.131], AUC 0.487).

This grid is diagnostic. It authorises no scoring change on its own; Stage 04
still ranks on `pve_cis_joint_calibrated` and must not be changed before
`results/combined/regime-verdict.tsv` is read.

Stage 01 was refactored onto the shared reader on 2026-08-22. Rerunning task 3
of `lgv-AA-caudate-20260822` reproduced its stored row on every scientific
field; only `bslmm_elapsed_sec` (wall clock) differed.

### Grid result (2026-08-22) and the scoring change it authorised

`lgv-observed-regime-20260822` completed 9,216 of 9,216 scenarios with zero
failures and zero unaccounted. Both questions resolved in favour of the
unbounded estimate, reversing what the AR(1) grid had indicated.

The support mismatch closes. Regime estimates span [-0.246, 0.982] against
observed [-0.267, 0.980]; 0.07% of observed loci fall below the regime
minimum, versus 56.4% below the AR(1) minimum.

Q1. The lower-boundary rate is an expected resolution limit, not a defect.
On real cis-window genotypes at n = 153 it is 0.998 at true PVE 0 and 0.05,
0.988 at 0.10, 0.912 at 0.20, 0.148 at 0.40 and ~0 above; pooled over
PVE <= 0.1 it is 0.995. The observed 62.7% sits inside that curve, which is
what a catalog mixing low- and moderate-PVE loci should give. The rate is
flat across `num_snps` strata (0.500--0.517) and `p_eff` strata
(0.489--0.518); the inverted `num_snps` trend seen against the AR(1) grid was
itself an AR(1) artifact.

Q2. Below the boundary the unbounded estimate carries real rank information:
Spearman 0.584, Kendall 0.452, bootstrap 95% CI [0.563, 0.603], AUC 0.741,
Wilcoxon p 8e-134, across 4,661 clipped scenarios. On the AR(1) grid the same
statistics were -0.0025 and [-0.133, 0.131]. Clipping costs about 80% of the
recoverable ordering where most real VMRs sit: Spearman at true PVE <= 0.1 is
0.312 unbounded against 0.059 clipped (overall 0.939 against 0.922).

Two changes follow, both applied 2026-08-22:

- Stage 04 ranks `local_snp_contribution_score` on `pve_cis_joint_unbounded`
  and records the basis in `local_snp_contribution_score_basis`. Stage 05
  counts ties on that same column. `pve_cis_joint_calibrated` is unchanged and
  remains the bounded descriptive output; `absolute_pve_interpretation_allowed`
  stays FALSE on every row.
- Stage 03 gates `num_snps`, `p_eff` and `ld_metric` against
  `config/joint-pve-characterized-support.tsv`, the union of the AR(1)
  training grid and the observed-regime grid, regenerated by
  `_h/10_write_characterized_support.R`. The previous gate bounded `p_eff`
  only by [1, n] and so flagged nothing; the union support (p_eff
  7.578--274.458) flags 0.92% of observed loci, well inside the 10%
  `max_outside_calibration_domain` threshold.

The grid covers caudate AA only. The other five cohort x region cells should
be characterised before their runs are interpreted.


### Do not edit `_h/` while an array is in flight

Four Stage 01 tasks of `lgv-AA-caudate-20260822` failed with
`cannot open file '01_estimate_observed_joint_features.R'` because a git
operation rewrote the script underneath the running array.

### Per-cell characterization (`11..13` + `submit_observed_regime_cell.sh`)

The grid stratifies loci on `num_snps` and `p_eff`, which for caudate AA came
from its completed production run. The other five cells have no production
run, and requiring one first would invert the dependency: Stage 03 needs the
characterized support, which needs the grid, which would need the run.

`p_eff` and `ld_metric` are pure genotype geometry, so `11..13` compute them
straight from the accepted Module 01 catalog at about 0.5s per locus -- no
phenotype, no elastic net, no BSLMM. `07` accepts either that scan or a
production feature table as its locus universe.

Submitted 2026-08-22, one chain per cell (geometry -> combine -> manifest ->
9,216 scenarios -> summary):

| Cohort | Region | n | Geometry run | Regime run |
| --- | --- | --- | --- | --- |
| AA | dlpfc | 118 | `lgv-geometry-AA-dlpfc-20260822` | `lgv-observed-regime-AA-dlpfc-20260822` |
| AA | hippocampus | 117 | `lgv-geometry-AA-hippocampus-20260822` | `lgv-observed-regime-AA-hippocampus-20260822` |
| all_individuals | caudate | 282 | `lgv-geometry-all_individuals-caudate-20260822` | `lgv-observed-regime-all_individuals-caudate-20260822` |
| all_individuals | dlpfc | 173 | `lgv-geometry-all_individuals-dlpfc-20260822` | `lgv-observed-regime-all_individuals-dlpfc-20260822` |
| all_individuals | hippocampus | 177 | `lgv-geometry-all_individuals-hippocampus-20260822` | `lgv-observed-regime-all_individuals-hippocampus-20260822` |

All five `n` values are already in the AR(1) training grid, so the `n` arm of
the domain gate does not reject them. `p_eff` is the arm that will move: the
AA dlpfc geometry smoke returned `p_eff` around 18--20 against an AR(1)
minimum of 24.9, so the mismatch is not specific to caudate.

### Characterized support widened over all six grids (2026-08-23)

All six observed-regime grids completed at 9,216/9,216 scenarios with zero
unaccounted scenarios and zero computational failures, and
`sub_boundary_rank_informative` is TRUE in every cell. The clipping cost
reproduces cohort-wide -- Spearman at true h2 <= 0.1, unbounded versus clipped:

| cell | overall (unbnd/clip) | low h2 (unbnd/clip) |
| --- | --- | --- |
| AA caudate | 0.939 / 0.922 | 0.312 / 0.059 |
| AA dlpfc | 0.924 / 0.910 | 0.248 / 0.077 |
| AA hippocampus | 0.923 / 0.908 | 0.241 / 0.030 |
| all_individuals caudate | 0.962 / 0.942 | 0.442 / 0.044 |
| all_individuals dlpfc | 0.946 / 0.933 | 0.328 / 0.050 |
| all_individuals hippocampus | 0.949 / 0.933 | 0.350 / 0.070 |

`boundary_rate_reproduced` is FALSE for the five new cells only because
`observed_boundary_rate` is NA there -- those cells had no production run when
the grid was summarized, so the comparison had no left-hand side. It becomes
evaluable once the 2026-08-23 production runs land. AA caudate, the only cell
with a production run at the time, reports TRUE.

`_h/10_write_characterized_support.R` now takes a comma-separated
`--regime-run-id` list and unions the AR(1) training grid with all six
observed-regime grids. The support widened to:

| feature | before | after |
| --- | --- | --- |
| num_snps | [100, 12000] | [100, 12000] |
| p_eff | [7.578, 274.458] | [2.058, 274.458] |
| ld_metric | [0.01183, 0.41251] | [0.01183, 1.0] |

`allowed_n` is now `117,118,153,173,177,282`.

Both new extremes come from one pair of loci at chr17:45.58-45.59 Mb -- the
17q21.31 MAPT inversion, where recombination suppression makes the cis-window
effectively two haplotypes (p_eff ~ 2, adjacent LD 1.0). These are real
genotype geometry, not a QC artifact: both loci converged on every estimator
and the frozen model ranks them against true PVE at Spearman 0.984 and 0.973.
The support is therefore the raw union, untrimmed.

## Production runs 2026-08-23

Because run directories are permission-locked immutable, the caudate run was
not re-scored in place. All six cells were submitted fresh as
`lgv-<cohort>-<region>-20260823`, so every cell carries identical provenance
and pins the same widened support
(`config_characterized_support_sha256 = f52026944b0346fa928a5e5cb9293dc742d2d4f5c8bacfea9913ec8854d2d9f7`).

Under the widened support, AA caudate's domain exclusions fall from 111 to 2
and its within-domain rate rises from 0.9902 to 0.99982; the other five cells'
locus geometry covers the low-`p_eff` and high-LD tail that the caudate-only
support had excluded.

`lgv-AA-caudate-20260822` is superseded and must not be entered in the
accepted-runs table.

## Production runs 2026-08-23: all six cells pass

All six cells return `PASS_RELATIVE_SCORE_OBSERVED_QC` on all six criteria,
with zero unaccounted tasks and zero computational failures, and all six run
directories are finalized and permission-locked.

| cell | expected | eligible | within-domain | boundary rate | max tie |
| --- | --- | --- | --- | --- | --- |
| AA caudate | 11,530 | 11,341 | 0.99982 | 0.6263 | 0.0001 |
| AA dlpfc | 9,572 | 9,341 | 0.99936 | 0.6225 | 0.0001 |
| AA hippocampus | 9,497 | 9,265 | 0.99925 | 0.5907 | 0.0001 |
| all_individuals caudate | 11,463 | 11,238 | 0.99956 | 0.5707 | 0.0001 |
| all_individuals dlpfc | 9,374 | 9,125 | 0.99912 | 0.5613 | 0.0001 |
| all_individuals hippocampus | 9,365 | 9,120 | 0.99923 | 0.5512 | 0.0001 |

`score_basis` is `pve_cis_joint_unbounded` in every cell.

`boundary_rate_reproduced` now evaluates TRUE in all six: each observed
boundary rate falls inside its own cell's simulated curve, between the
regime-wide overall rate and the maximum cell rate. The observed rates
(0.551--0.626) sit well below the low-h2 plateau (~0.99), which is what an
expected resolution limit looks like rather than a simulation-to-data mismatch.

### Seed-overflow defect found and fixed

The first submission failed two cells on `zero_computational_failures` only.
`stable_seed()` reduces modulo 2147483629, only 18 below
`.Machine$integer.max`, and `00_functions.R` derived per-repeat and per-fold
seeds by addition (`seed + repeat_id * 1009L`, `+ fold * 9173L`). A base seed
in the top band overflowed to NA and `set.seed(NA)` threw. Three loci died this
way: chr3:131957079-131957863 and chr9:136552488-136552675 in AA caudate,
chr11:59101497-59102234 in all_individuals dlpfc. The failure is deterministic,
so a plain resubmission would have reproduced it.

`offset_seed()` now does the arithmetic in double precision and takes one
modulo. Verified across all 60,801 loci and every repeat/fold offset in use:
363 previously-NA derived seeds become valid and zero previously-working seeds
change, so the four already-passing cells are bit-identical under the fix.

## Sequence still to run

1. Hand-enter the six 2026-08-23 runs in the accepted-runs table below.
   AGENTS.md section 6 makes this a human step, and no downstream module may
   read Module 02 until it is done. `lgv-AA-caudate-20260822` is superseded and
   must not be entered.

Any run's domain decision must be read against the
`config_characterized_support_sha256` recorded in its own manifest, not against
the working-tree file, which widens as cells are added.

## Scheduler-safe reorganization

Before the 2026-08-21 reorganization, `squeue -u $USER` showed no running or
pending Module 02 job. Other modules had pending/running work and were not
touched. The move changed source locations only; immutable `_m/runs/` outputs
were not modified.

## Archived scripts

Sixty-one historical estimator-development, validation, recovery, and obsolete
observed-EN entry points now live in
`_h/archive/2026-08-21_estimator_development/`. They are retained for forensic
review and later removal, not as an alternative production pipeline. See its
`README.md` and `ARCHIVE_MANIFEST.tsv` for the removal gate.

Completed historical run directories retain their own executed code/config
snapshots. No archived path may be called by an active launcher or downstream
module.

## Acceptance gate

For each cohort-by-region cell, acceptance requires:

1. an accepted, non-smoke Module 01 run and identical `vmr_set_id`;
2. the exact frozen joint-model checksum above;
3. deterministic complete feature calculation and task reconciliation;
4. zero unexplained computational failures;
5. at least 90% within-domain among complete-feature loci;
6. a finite, nondegenerate within-cell score;
7. `absolute_pve_interpretation_allowed = FALSE` on every row;
8. Stage 05 decision `PASS_RELATIVE_SCORE_OBSERVED_QC`;
9. immutable Stage 06 checksums and a manual README acceptance record.

Poor domain coverage or widespread boundary/tie concentration is reported, not
hidden by pooling cells or relaxing the score definition.

## Accepted runs

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| lgv-AA-caudate-20260823 | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-08-23 | Kynon J.M. Benjamin | PASS_RELATIVE_SCORE_OBSERVED_QC | Relative rank only; absolute PVE prohibited |
| lgv-AA-dlpfc-20260823 | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-08-23 | Kynon J.M. Benjamin | PASS_RELATIVE_SCORE_OBSERVED_QC | Relative rank only; absolute PVE prohibited |
| lgv-AA-hippocampus-20260823 | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-08-23 | Kynon J.M. Benjamin | PASS_RELATIVE_SCORE_OBSERVED_QC | Relative rank only; absolute PVE prohibited |
| lgv-all_individuals-caudate-20260823 | all_individuals | caudate | vmrset-all_individuals-caudate-cb5519d7d2ad | 2026-08-23 | Kynon J.M. Benjamin | PASS_RELATIVE_SCORE_OBSERVED_QC | Relative rank only; absolute PVE prohibited |
| lgv-all_individuals-dlpfc-20260823 | all_individuals | dlpfc | vmrset-all_individuals-dlpfc-e88f46904afb | 2026-08-23 | Kynon J.M. Benjamin | PASS_RELATIVE_SCORE_OBSERVED_QC | Relative rank only; absolute PVE prohibited |
| lgv-all_individuals-hippocampus-20260823 | all_individuals | hippocampus | vmrset-all_individuals-hippocampus-809f8de0db2d | 2026-08-23 | Kynon J.M. Benjamin | PASS_RELATIVE_SCORE_OBSERVED_QC | Relative rank only; absolute PVE prohibited |

No downstream production module may consume Module 02 until its cell-specific
run appears in this table.

`lgv-AA-caudate-20260822` is superseded by `lgv-AA-caudate-20260823` and is
deliberately absent: it was scored under the caudate-only characterized
support, which excluded 111 loci where the six-grid support excludes 2.

Every accepted run carries `absolute_pve_interpretation_allowed = FALSE` on
every row and the terminal decision
`PASS_RELATIVE_GENETIC_CONTROL_FAIL_ABSOLUTE_LOCUS_PVE`. Downstream modules may
consume `local_snp_contribution_score` and its percentile rank, never the PVE
magnitude.

## Verification status

- Scheduler audit: no active Module 02 jobs at reorganization time.
- Frozen simulation validation: pass for relative ordering, fail for absolute
  PVE.
- Active scripts: R and shell syntax checks pass; active-pipeline integrity,
  Stage 01 chunk-wrapper, frozen joint-estimator, and score/loader tests pass.
- Archive reconciliation: 61 manifest rows, 61 retained scripts, zero missing
  files, and no retired-path calls from active launchers.
- Real-data smoke: Stage 01 completed nested EN, HE, BSLMM, effective-rank, and
  LD features for two 2,672--3,009-SNP caudate VMRs.
- End-to-end smoke: `lgv-AA-caudate-20260822b`, a 10-locus evenly spaced task
  universe, ran Stages 00--06 through the scheduler as six chained step jobs.
  Reconciliation was 10 expected / 9 completed / 1 qc_failed / 0 unaccounted /
  0 computational failures; Stage 05 returned
  `PASS_SMOKE_ONLY_NOT_ACCEPTABLE` and Stage 06 sealed the run. Five of the
  nine scored loci sat at the lower boundary (`max_tie_fraction` and
  `boundary_rate` both 0.556). At smoke scale that is uninformative, but
  boundary and tie concentration is the property to watch in production,
  because heavy zero mass compresses the within-cell rank score.
- Fail-closed smoke: an intentionally incomplete 507-task universe returned
  `FAIL_OBSERVED_RELATIVE_SCORE_QC` with 505 unaccounted/computational failures.
- Accepted observed production score runs: none.

## `_m/` contents

Generated runs live only under immutable `_m/runs/{RUN_ID}/`. Each active run
contains its manifest, frozen configuration, task/chunk manifests, per-task
terminal rows, reconciliation, combined features and score, QC decision,
session record, and SHA-256 output manifest. See `_m/README.md`.
