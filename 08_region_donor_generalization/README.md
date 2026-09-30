# 08_region_donor_generalization — cross-region robustness and identifiable regional heterogeneity

Establishes what **reproduces** across brain regions, which regional difference
is actually **identified**, and what the donor-group and matched-subset
contrasts can support.

**Status: accepted (AA, 2026-09-18), and its upstreams have since moved — see
"The accepted run now rests on superseded upstreams" below.**
`rdg-AA-crossregion-20260918` passed the gate and is recorded under **Accepted
runs**; criterion 10 was added the next day and it passes retrospectively. The
blocking upstreams were satisfied when it opened —
`04_repeat_repressive_architecture` (`rra-AA-*-20260906`) and
`05_cpg_meqtl_burden` (`cmb-AA-*-20260825`) both recorded passing acceptance
gates (AGENTS.md §6) — and the donor-group axis is unblocked by the six accepted
cell runs in each of `01b_estimation_cells`, `02_local_genetic_variance` and
`03_local_snp_prediction`.

Tier 3 is the one axis that needed new upstream compute rather than assembly:
`config/cohorts.yml` declares `AA.n118r{1,2,3}`, and their 01b → 02 → 03 chains
had to be sealed **and accepted** before this module would open a run. All nine
are accepted (2026-09-18).

## The accepted run now rests on superseded upstreams

**Added 2026-09-27.** Four of `rdg-AA-crossregion-20260918`'s upstreams were
superseded by the 2026-09-24/25 acceptances, so its manifest points at runs that
are no longer the ones a new production run may consume:

| upstream | in the accepted 08 run | now accepted |
|---|---|---|
| `03_local_snp_prediction` | `lsp-AA-*-20260825` | `lsp-AA-*-20260925-a` |
| `04_repeat_repressive_architecture` | `rra-AA-*-20260906` | `rra-AA-*-20260925-a` |
| `05_cpg_meqtl_burden` | `cmb-AA-*-20260825` | `cmb-AA-*-20260924` |
| `07_transcription_splicing_coupling` | `tsc-AA-*-20260902` | `tsc-AA-*-20260925-b` |

AGENTS.md §6 is not retroactive: the 2026-09-18 acceptance was valid when it was
made and stays in the table. What follows is narrower and is the operative rule
here — **no new production run may consume `rdg-AA-crossregion-20260918`**, and
its tier counts may not be quoted beside numbers from the 2026-09-25 upstreams.

### A rerun exists, passed its gate, and was deliberately not accepted

`rdg-AA-crossregion-20260925` (sealed 2026-09-25 at commit `e93319890`,
`git_dirty = false`) consumes all four replacements and returned
`PASS_REGION_DONOR_GENERALIZATION_QC` on 10 of 10 criteria. It is **not** in the
Accepted runs table, because four accounting defects were found in the tier
counts after it sealed. The gate did not catch them and could not: none of its
ten criteria inspects claim-family composition or the tier-2 denominator, which
is consistent with §7.7 — the gate certifies tiering and interpretation
constraints, not a positive finding.

The run reports *15 of 16 claim-family tests replicate, 9 strict* and *39
identified DLPFC-vs-hippocampus differences*, against 12/16, 12 strict and 1
difference in the superseded run. Both directions of that movement are artifacts:

1. **`expression_abc` enters the claim family.** *Repaired 2026-09-30 (T20).*
   `harvest()` defaults `outcome_role` to `prespecified_family` and the Module 07
   harvest spec does not map the field, so Module 07's `in_fdr_family = FALSE` —
   set by the locked per-modality power floor — is ignored here. Two of three ABC
   tests then count as replicated and strict, which is the whole of the 12 → 15
   rise. One carries `p = 0` exactly. The spec now maps the flag through a new
   `fdr_family_flag` key onto a `power_excluded` role, which `in_claim_family`
   does not admit; the flag's type is checked rather than trusted, because the
   NA-safe `%in% TRUE` idiom answers FALSE for a character `"TRUE"` and would
   empty the claim family silently.
2. **A one-region arm makes strict replication unsatisfiable.**
   `adjust_cell_composition_scmd` is correctly caudate-only (the scMD integration
   gate passes only there), so `n_regions = 1`, while `complete_across_regions`
   requires all three. Every one of the six repeat tests therefore reads
   `replicated_strict = FALSE`. That is the whole of the 12 → 9 fall, and it is an
   artifact of a conjunction over an arm that cannot be fitted, not lost
   robustness.
3. **A null estimate's sign vetoed a replication.** *Repaired 2026-09-30 (T21).*
   `line_l1_frac × score_z × exclude_segdups` flipped TRUE → FALSE only because
   the caudate estimate moved from +0.0165 (p = 0.745) to −0.0004 (p = 0.994),
   while DLPFC holds at p = 5.1e-11 and hippocampus at p = 4.1e-14.
   `direction_consistent := pmax(n_up, n_down) == n_regions` counts the sign of a
   non-significant estimate — there the SE is 140× the estimate — so a meaningless
   sign can veto a replication. Direction is now judged on the nominally
   supported regions only: a null estimate carries no direction and neither
   supports nor contradicts one, a region that contradicts **with** support still
   vetoes, and consistency among zero supported regions is FALSE rather than
   vacuously TRUE. The superseded rule is retained beside it as
   `direction_consistent_all_regions` so the two verdicts stay comparable.
4. **Tier 2 applies no `outcome_role` filter and no primary-arm filter.** Of its
   39 "genuine regional heterogeneity" differences, 38 are the `atac_*` controls
   registered under T8/T11, 28 rest on `r2_pred_oof_z` — which Module 04's own
   claims table forbids citing as corroboration — and each sensitivity refit is
   counted as a separate difference (9 primary + 30 refits), the inflation `_h/01`
   refuses by design. Primary-arm, non-control differences: **1**, unchanged from
   the superseded run.

Under all three tier-1 fixes the corrected reading is **13 of 13 replicate, 13
strict** — cleaner than either reported number, which is why nothing is lost by
leaving the run unaccepted until the code is fixed.

**Fixes 1 and 3 are implemented (2026-09-30).** They are defect repairs in
`_h/01_cross_region_replication.R` that change no configuration, so they carry no
PI decision. Fix 2 and the tier-2 filters change reported counts and touch
`config/region_donor_generalization.yml`, which is `pi_locked`, so they remain PI
decisions under AGENTS.md §12 and are not made here.

Re-running tier 1 against the same pinned upstreams with fixes 1 and 3 in place —
outside `_m/`, since `rdg-AA-crossregion-20260925` is sealed — gives **13 of 13
claim-family tests replicate, 7 strict**, and moves exactly the rows it should:

- 9 rows change `outcome_role` (the three `expression_abc` tests × three regions)
  and nothing else does, so `n_claim_tests` falls 16 → 13 and
  `n_power_excluded_tests` reads 3;
- 7 rows change `direction_consistent`, of which **one changes a verdict**:
  `exclude_segdups × line_l1_frac × score_z` becomes `replicated`. The other six
  lose a vacuous TRUE at zero nominal support and were already not replicating;
- **no row loses `replicated` and none loses `replicated_strict`.** The
  `cross-region-tests.tsv` table differs only in those 9 role cells, and
  `cross-region-rank-agreement.tsv` is byte-identical.

The strict count reads 7 rather than 13 because fix 2 is still open: the two
`expression_abc` strict artifacts leave the claim family (9 → 7) while the six
repeat tests stay blocked by the caudate-only `adjust_cell_composition_scmd` arm.
Fix 3's contribution is visible in `n_sensitivity_replicated` for the `line_l1_frac
× score_z` family, which rises 3 → 4 of 5; the unsatisfiable fifth is fix 2's, so
**strict reaches 13 only once the PI closes fix 2.**

Regression test: `tests/test_tier1_role_and_direction.R`. It evaluates the shipped
`harvest()` body against all three real Module 07 tables and the shipped
`direction_consistent` expressions against a five-case fixture plus the sealed
run's own rows, so a future edit that reintroduces either defect fails there
rather than in a production count.

## Accepted runs

Machine-readable, in the schema `00_shared/gates.R::read_accepted_runs()` parses.
A run of this module spans all three regions, so `region` is the literal
`crossregion` and the per-region `vmr_set_id`s are recorded in the manifest as
`vmr_set_id_{region}` rather than in this table.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| rdg-AA-crossregion-20260918 | AA | crossregion | see manifest vmr_set_id_{caudate,dlpfc,hippocampus} | 2026-09-18 | Kynon J.M. Benjamin | PASS_REGION_DONOR_GENERALIZATION_QC | 9/9 gate criteria; 14 outputs; built at 6c2a930e4. Tier 1: 12 of 16 prespecified claim-family tests replicate in all 3 regions (12 strict, all 154 tests complete), both specificity controls run opposite the claim family. Tier 2: 1 of 154 testable dlpfc-minus-hippocampus differences survives strict conjunction; 154 caudate rows retained descriptive-only. Tier 3: primary within-caudate paired delta on 11335 shared loci, A = 0.1426 (block-jackknife 95% CI 0.131-0.155, reported not gated), all 3 replicates same direction, gap_closed 0.984; reading donor_count_is_a_plausible_major_contributor; lower-boundary mass rises 0.6247 to 0.6444 (+2.0 pts), reported alongside and excluded from the 0.10 threshold. Donor-group axis: concordance only, rho 0.815/0.753/0.759 = 85.3%/88.9%/88.0% of the reliability ceiling; no ancestry effect claim. Caudate remains batch-confounded; residual excess may NOT be called biological. |

## Pipeline

| stage | tier | writes |
|---|---|---|
| `00_new_run.R` | — | `results/tiers.tsv`; gates all 30 region-axis + 18 cell + 9 tier-3 upstreams and pins their run IDs |
| `01_cross_region_replication.R` | 1 | `cross-region-{tests,replication,rank-agreement,summary}.tsv` |
| `02_identified_difference.R` | 2 and 4 | `identified-difference{,-summary}.tsv`, `descriptive-confounded-regions.tsv` |
| `03_donor_group_concordance.R` | donor group | `donor-group-concordance.tsv` |
| `04_caudate_downsampling.R` | 3 | `caudate-downsampling-{replicates,summary}.tsv` |
| `05_apply_gates.R` | — | `region-donor-generalization-{qc,decision}.tsv` |
| `06_finalize_run.R` | — | `interpretation-constraints.txt`, seals |

Run it with `COHORT=AA _h/submit_region_donor_generalization.sh`. Every stage is
cheap: the module assembles accepted upstream results and refits nothing.

## Acceptance gate

`config/region_donor_generalization.yml:gate` locks ten criteria and the
terminal decision `PASS_REGION_DONOR_GENERALIZATION_QC`:

1. `every_output_carries_exactly_one_tier`
2. `cross_region_replication_computed_for_all_three_regions`
3. `identified_difference_restricted_to_dlpfc_hippocampus`
4. `caudate_downsampling_replicates_complete`
5. `donor_group_policy_is_concordance_only`
6. `no_pooled_rank_or_pooled_r2_emitted`
7. `no_cross_region_raw_score_comparison`
8. `confounded_caudate_reported_separately_not_dropped`
9. `reliability_ceiling_attached_to_every_concordance_figure`
10. `cross_region_completeness_nonvacuous` — added 2026-09-19. Fails when no
    test, or no claim-family test, is observed in every region, and when the
    completeness flag disagrees with a row-by-row re-derivation. A null
    replication passes; a replication over zero complete tests does not. The
    accepted `rdg-AA-crossregion-20260918` was gated under nine criteria and
    passes this one retrospectively (154 of 154 complete, 16 claim tests).

**These criteria certify tiering and interpretation constraints, not a positive
finding.** A null cross-region replication and a null region difference both
pass, exactly as `06_partitioned_heritability` passed with
`sldsc_supports_brain_enrichment = FALSE`. A gate that required a finding would
be a gate on the conclusion.

## Migrating from

Donor-group, cross-region, and downsampling modules: `meqtl-validation/04_cross_region_sharing/`, `05_donor_group_comparison/`, `10_downsampling_caudate/`, `13_vmr_universe_nmatched/`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Scope, set by the PI 2026-09-06 before implementation

This module is **not** a symmetric "shared versus context-dependent" analysis.
Brain region is perfectly confounded with sequencing batch (AGENTS.md §8.1;
`writing-notes/WGBS_BATCH_REGION_CONFOUNDING.md` carries the evidence and the
full design review), so replication and difference are not two equal findings
here. Every output carries exactly one tier:

| tier | analysis | licenses |
|---|---|---|
| **Strongest** | cross-region replication | robustness of genetic-control architecture across technical and regional contexts |
| **Identified difference** | DLPFC vs hippocampus, both within batches 1–2 | genuine regional heterogeneity |
| **Mechanistic sensitivity** | caudate downsampled to n=118 | whether donor count explains the caudate excess — and nothing more |
| **Descriptive only** | caudate vs other regions | nothing; reported because readers will ask |

1. **Cross-region replication is the primary deliverable.** Describe it as
   robustness across technical and regional contexts. Do not write that the
   confounding strengthens the result: the batch boundary makes a successful
   replication more compelling, and does nothing for the interpretability of a
   difference. Existing material that belongs here and is currently only a
   by-product elsewhere: Module 03's check C, Module 05's 3% beta agreement
   across cells, Module 07's expression coupling.
2. **DLPFC vs hippocampus is the only clean region-difference analysis.** Use
   the Module 04 template: primary claim on the identified contrast,
   strict-conjunction sensitivity gating, confounded cell reported separately.
   The contrast is a two-sample z that adds the per-region variances, and the two
   regions **share 115 of 118 donors** (Jaccard 0.96), so the estimates are
   positively correlated and the true variance of the difference is smaller than
   the one used. The test under-rejects, which is the safe direction for a tier
   licensed to claim heterogeneity; a comment in `_h/02_identified_difference.R`
   asserted independence until 2026-09-19 and was wrong. The paired donor
   bootstrap in `09b_aging_application/_h/05_cross_region_concordance.R` is the
   estimator that uses the covariance, and adopting it here is a scope change.
3. **Caudate downsampled to n=118** tests Module 03's untested attribution of
   the caudate excess to donor count (153 vs 118). Read it in two directions
   only — the excess largely disappears (donor count is a plausible major
   contributor) or it persists (donor count does not explain it). It **cannot**
   establish that a residual is biological; caudate stays batch-confounded
   whatever it shows. This limit belongs in the module's interpretation
   constraints, not in the reader's head.

   **Primary endpoint, PI 2026-09-18.** The primary tier-3 result is the
   **within-caudate** change on the *identical shared caudate locus set* — full
   n=153 against each n=118 replicate, paired on locus. Both terms are the same
   region and the same loci, so no cross-region reference enters it.
   `config/region_donor_generalization.yml:caudate_downsampling.primary_endpoint`
   locks it as `within_caudate_paired_delta_r2`.

   **Three locked criteria, PI 2026-09-18.** (i) and (ii) are the primary and
   are evaluated on the shared caudate loci alone; (iii) is separate:

   | # | criterion | config key | kind |
   |---|---|---|---|
   | i | all three replicates attenuate in the same direction | — | direction |
   | ii | mean relative attenuation `A >= 0.10`, where `A = mean_r (R²_full − R²_n118,r) / R²_full` | `primary_min_relative_attenuation` | **magnitude / effect size** |
   | iii | `gap_closed >= 0.50` — licenses calling donor count a plausible major contributor *to the caudate–DLPFC difference* | `major_contributor_gap_closed_min` | **interpretive** |
   | — | replicate spread in `A_r` within 0.10 | `replicate_agreement_max_range` | **agreement tolerance** |

   **None of these is a significance cutoff.** An effect-size criterion was
   chosen over a paired test deliberately: with ~11k paired loci a trivial
   attenuation would be "significant" and would still say nothing about whether
   donor count matters. Uncertainty on `A` is **reported, not gated** — a
   delete-one-chromosome weighted block jackknife (Busing et al. 1999, the
   construction `06_partitioned_heritability` uses for block standard errors),
   because loci within a chromosome are not independent.

   **The gap ratio is secondary and descriptive.** Writing it out,

   ```
   fraction_of_excess_closed_by_matching_n
       = (mean_r2_full - mean_r2_subset) / (mean_r2_full - dlpfc_reference)
   ```

   the numerator is the primary within-caudate attenuation and is clean — the
   DLPFC reference cancels. The **denominator is not**: caudate and DLPFC carry
   different `vmr_set_id`s, so no cross-region locus intersection exists. The
   caudate means are restricted to the loci shared by the full run and all three
   replicates, while the DLPFC mean is over all DLPFC-scored loci, unrestricted.
   Selection into the shared set is not random with respect to r², so the
   denominator carries an uncontrolled term. Use the ratio **only** to say how
   much of the region-level gap the primary attenuation would represent; it
   cannot on its own decide whether donor count explains the excess. Stage 04
   enforces this: the reading requires a replicate-consistent primary
   attenuation before the ratio is consulted at all.

   **Estimator resolution is retained as an explicit sample-size sensitivity.**
   Module 02's `boundary_rate` — the fraction of eligible loci whose unbounded
   estimate sits at the frozen model's output floor, i.e. loci with no
   detectable local genetic control — rises with the draw-down:

   | cell | n | boundary rate |
   |---|---|---|
   | `lgv-AA-caudate-20260823` | 153 | 0.6263 |
   | `AA.n118r1` | 118 | 0.6431 |
   | `AA.n118r2` | 118 | 0.6445 |
   | `AA.n118r3` | 118 | 0.6456 |

   In caudate this is **entirely the lower boundary** — zero upper-boundary hits
   in the arm or any replicate — so removing 35 donors pushes ~1.8% more loci
   below the floor, consistently across draws (spread 0.0025). **Any attenuation
   after downsampling therefore partly reflects statistical resolution rather
   than a biological change**, and must be reported that way. The per-replicate
   numbers ride on the tier-3 output rows, and
   `boundary_shift_is_estimator_resolution_not_biology` is carried beside them so
   the caveat cannot be separated from the number.

   It is reported **alongside** the attenuation and is **never folded into the
   0.10 threshold**: `A` is computed on prediction r² alone and carries no
   boundary-rate correction (`boundary_shift_excluded_from_attenuation_threshold`).
   That is what keeps overall prediction attenuation distinguishable from the
   accompanying loss of estimator resolution, instead of silently mixing the two
   into one number.
4. **Caudate-vs-other differences stay descriptive.** Reuse the Module 04
   mechanism — `interpretation.technically_confounded_regions` sets the cell
   aside from the claim while keeping the estimate fitted, written and surfaced
   in dedicated columns.

### Donor-group axis

**PI decision 2026-09-06:** follow the v1 approach — compare AA and EA local
genetic-control estimates **on the common pooled-discovery VMR set**, the
`all_individuals` catalog, rather than on cohort-specific catalogs. Discovery
happens once in the pooled sample; the donor groups are then two disjoint sets
of donors evaluated on one fixed, shared locus set, so neither group's VMR
calling can advantage it.

- Do **not** contrast `AA` against `all_individuals`. Those are nested —
  `all_individuals` contains the AA donors — and a set-versus-superset
  comparison is not a donor-group contrast.
- **The EA estimation cell now exists.** `01b_estimation_cells` materializes
  `all_individuals.AA` and `all_individuals.EA` from the sealed pooled Module 01
  runs (six accepted runs, 2026-09-10), Modules 02 and 03 have six accepted cell
  runs each (2026-09-17 / 2026-09-18), and both recombination stages have run
  for all three regions. This axis is unblocked; see AGENTS.md §7.7.
- This axis is **not** exposed to the batch confounding: donor groups interleave
  within each region's libraries.
- Expect unequal n. Report it, and match or downsample as a tier-3 sensitivity
  rather than adjusting for it post hoc.

**Inference policy, PI 2026-09-18.** `config/analysis_thresholds.yml:donor_group`
sets `donor_group_inference: concordance_only`,
`cross_group_raw_score_comparison: false` and
`ancestry_effect_claim_allowed: false`. It replaced
`require_interaction_for_ancestry_claim: true`, which named a test this design
cannot run — an interaction term needs a pooled group × genotype model, and §7.6
forbids comparing the Module 02 score across cells at any level — so the key was
unsatisfiable and guarded nothing. `00_shared/gates.R::donor_group_inference_policy()`
reads and enforces the three keys, and refuses a widened policy rather than
defaulting.

Concordance is read **against the analytic reliability ceiling** from
`02/_h/16_reliability_ceiling.R`. Without the ceiling an imperfect ρ reads as a
donor-group difference when most of it is input uncertainty, which is the single
most likely misreading of this axis.

## Priorities

Prioritize biological generalization of local variance, repeat enrichment, meQTL
burden, and effect direction. Cross-population predictor portability is optional
and secondary.

AGENTS.md §7.6 already forbids raw score-level comparison across regions, so the
Module 02 score cannot be contrasted across cells at all. The region axis was
always narrower than the migration sources imply.

Do not attribute differences to ancestry-specific biology without eliminating
sample size, MAF, LD, SNP availability, assay, covariate, and brain-region
explanations. Use donor-group or population language approved by the PI:
`AA` = "Black American", `EA` = "non-Hispanic white American".

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
