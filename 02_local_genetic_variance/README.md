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

The internal `pve_cis_joint_calibrated` estimate is converted within each
cohort-by-region cell to:

- `local_snp_contribution_score`: empirical midrank percentile,
  `(midrank - 0.5) / n_eligible`;
- `local_snp_contribution_score_z`: standardized score used in primary models;
- `local_snp_contribution_quartile`: secondary top/bottom-quartile contrast.

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
| _(none)_ | | | | | | | |

No downstream production module may consume Module 02 until its cell-specific
run appears in this table.

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
- End-to-end smoke: Stages 00--06 completed and correctly returned
  `PASS_SMOKE_ONLY_NOT_ACCEPTABLE` when run on the two-locus task universe.
- Fail-closed smoke: an intentionally incomplete 507-task universe returned
  `FAIL_OBSERVED_RELATIVE_SCORE_QC` with 505 unaccounted/computational failures.
- Accepted observed production score runs: none.

## `_m/` contents

Generated runs live only under immutable `_m/runs/{RUN_ID}/`. Each active run
contains its manifest, frozen configuration, task/chunk manifests, per-task
terminal rows, reconciliation, combined features and score, QC decision,
session record, and SHA-256 output manifest. See `_m/README.md`.
