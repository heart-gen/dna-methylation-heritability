# 08_region_donor_generalization — cross-region robustness and identifiable regional heterogeneity

Establishes what **reproduces** across brain regions, which regional difference
is actually **identified**, and what the donor-group and matched-subset
contrasts can support.

**Status: accepted (AA, 2026-10-08: `rdg-AA-crossregion-20261008`, gate 10/10,
0 primary claim-family tier-2 differences once the meQTL-slope contrast reads the
donor-robust SE).** The history below records the 2026-09-18 acceptance and the
reruns that superseded it; the current run is the one in **Accepted runs**.
`rdg-AA-crossregion-20260918` passed the gate and is recorded under **Accepted
runs**; criterion 10 was added the next day and it passes retrospectively. The
blocking upstreams were satisfied when it opened —
`04_repeat_repressive_architecture` (`rra-AA-*-20260906`) and
`05_cpg_meqtl_burden` (`cmb-AA-*-20260825`) both recorded passing acceptance
gates — and the donor-group axis is unblocked by the six accepted
cell runs in each of `01b_estimation_cells`, `02_local_genetic_variance` and
`03_local_snp_prediction`.

Tier 3 is the one axis that needed new upstream compute rather than assembly:
`config/cohorts.yml` declares `AA.n118r{1,2,3}`, and their 01b → 02 → 03 chains
had to be sealed **and accepted** before this module would open a run. All nine
are accepted (2026-09-18).

## Tier 2's one difference rested on VMR-level SEs (2026-10-08)

The accepted run's only tier-2 difference is the `meqtl_burden` score_z slope:
DLPFC 2.094 against hippocampus 2.367, delta −0.273, z −3.63. Its SE is the
square root of the sum of the two regions' Module 05 SEs, and those are
VMR-level HC3. HC3 treats the ~9,000 VMRs in a region as independent, and it
ignores the 115 donors the two regions share.

`cmb-AA-crossregion-20261007` (Module 05, TASKS A3) refits the same model with a
paired delete-d donor jackknife plus a joint chromosome jackknife.
- Both region estimates reproduce exactly.
- Each region's SE is about 3× its HC3 SE.
- The difference's SE is 0.204 against 0.075.

On that SE the difference is z −1.34, p 0.18. The two regions' slopes agree.

**The fix:** `identified_difference.donor_robust_se` in
`config/region_donor_generalization.yml` names a source by module, region and
test.
- Stage 00 gates the source with `require_accepted_upstream()`, so production
  cannot use it until the PI accepts it. Only a smoke run may name an
  unaccepted source, through `V2_SMOKE_DONOR_ROBUST_RUN`.
- Stage 02 refuses a source that does not reproduce both stage-01 estimates to
  1e-9. When it accepts one, the source's SEs replace the HC3 ones in three
  places:
  - the delta;
  - the per-region magnitude check;
  - the both-nominal check.
- `identified-difference.tsv` keeps the old SE as `delta_se_hc3_independent`,
  and `delta_se_source` names where each row's SE came from.

Every other tier-2 row keeps its upstream SE. Those SEs carry the same
limitation, but no donor-robust refit exists for them.

**What this did to the accepted run.** It could not be changed in place. The
rerun `rdg-AA-crossregion-20261008`, accepted 2026-10-08, reports **0** primary
claim-family differences, with the counts the smoke run predicted (39 to 37 rows
across all arms, the ATAC row's q 0.0493 to 0.0506). Do not quote the
2026-09-30 run's difference as regional heterogeneity.

A smoke run (`V2_SMOKE_DONOR_ROBUST_RUN=cmb-AA-crossregion-20261007`, deleted
after inspection) confirms this:
- The meQTL row gets SE 0.204, p 0.18 and `difference_claimed = FALSE`.
- The primary claim-family count goes from 1 to 0.
- The all-rows count goes from 39 to 37. BH runs over all 314 testable pairs, so
  the meQTL row's higher p moves other rows' q. One secondary row
  (`repeat_architecture` × `high_mappability` × `atac_opc_frac` × `r2_pred_oof_z`)
  goes from q 0.0493 to 0.0506 and loses its claim.
- No other row changes, and the gate passes 10 of 10.

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

The locked analysis plan is not retroactive: the 2026-09-18 acceptance was valid when it was
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

**All four are now repaired in `_h/` (2026-09-30), and that does not make this run
acceptable.** Its tables were written by the defective code and `_m/` is immutable,
so the corrected counts below describe what a rerun produces. Reacceptance needs a
fresh run ID (T23).

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
2. **A one-region arm makes strict replication unsatisfiable.** *Repaired
   2026-09-30 (T22).* `adjust_cell_composition_scmd` is correctly caudate-only (the
   scMD integration gate passes only there), so `n_regions = 1`, while
   `complete_across_regions` requires all three. Every one of the six repeat tests
   therefore reads `replicated_strict = FALSE`. That is the whole of the 12 → 9
   fall, and it is an artifact of a conjunction over an arm that cannot be fitted,
   not lost robustness. The conjunction now ranges over arms that exist
   cross-region. The remedy is forced, not chosen: judging the arm on the region
   where it *was* fitted would let a caudate-only arm veto a cross-region
   replication claim, and caudate is perfectly confounded with sequencing batch
   (§8.1), so a batch effect could decide a cross-region verdict. Arms outside the
   conjunction are counted in `n_sensitivity_sets_single_region` and keep their own
   `replicated` verdict, and `strict_conjunction_vacuous` flags the opposite
   failure — a test that would pass having survived nothing.
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
4. **Tier 2 applies no `outcome_role` filter and no primary-arm filter.**
   *Repaired 2026-09-30 (T22).* Of its 39 "genuine regional heterogeneity"
   differences, 30 are sensitivity refits and the rest are secondary scales,
   cell-type breakdowns and independent-assay contrasts; each sensitivity refit was
   counted as a separate difference (9 primary + 30 refits), the inflation `_h/01`
   refuses by design. Primary-arm, non-control differences: **1** here and in
   `-20260930`, but **0** in the superseded `-20260918` -- whose single difference
   was the `expression_abc` artifact fix 1 excludes. The count is not carried over
   and the test is not the same one; see "### Superseded". `_h/02` now carries
   `outcome_role` through its cast and reports the primary claim-family count as
   the headline, with the unfiltered count kept as the auditable denominator. `difference_claimed` is deliberately
   unchanged — it answers a statistical question about one row, and both the gate
   (`_h/05_apply_gates.R:121`) and Module 11's QQ panel read it — so the filters
   are additive columns, not a redefinition.

### All four are repaired (2026-09-30)

None of the four needed a configuration change. An earlier version of this section
said fixes 2 and 4 "touch `config/region_donor_generalization.yml`, which is
`pi_locked`", and that was wrong on the facts: the config declares the tiers, their
licences, the contrast, `require_strict_conjunction`, `alpha` and `fdr_method`, and
says nothing about which `analysis_set`s or `outcome_role`s enter a conjunction or
a count. There is no key to edit. The locked analysis plan turns on whether a change is a
silent scientific decision, not on whether it edits a YAML file, and none of the
four is: fix 1 reads a flag the source module already publishes, fix 3 stops a
number with an SE 140× its magnitude from carrying a direction, fix 2's remedy is
forced by §8.1, and fix 4 makes tier 2 count the way tier 1 in the same module
already counts. What remains a human act is accepting the rerun, which §6 requires
regardless.

Re-running tiers 1 and 2 against the same pinned upstreams — outside `_m/`, since
`rdg-AA-crossregion-20260925` is sealed — gives the corrected reading:

| | as sealed | corrected |
|---|---|---|
| claim-family tests | 16 | **13** |
| replicate | 15 | **13** |
| replicate strictly | 9 | **13** |
| tier-2 differences | 39 | **1** |

and moves only the rows it should:

- 9 rows change `outcome_role` (the three `expression_abc` tests × three regions)
  and nothing else does, so `n_claim_tests` falls 16 → 13 and
  `n_power_excluded_tests` reads 3;
- 7 rows change `direction_consistent`, of which **one changes a verdict**:
  `exclude_segdups × line_l1_frac × score_z` becomes `replicated`. The other six
  lose a vacuous TRUE at zero nominal support and were already not replicating;
- 6 rows gain `replicated_strict`, all of them primary claim-family repeat tests,
  each with 5 sensitivity sets of which 4 are cross-region and all 4 replicate;
- **no row loses `replicated` and none loses `replicated_strict`.**
  `cross-region-rank-agreement.tsv` is byte-identical, and in tier 2 both
  `difference_claimed` and `delta_z` are unchanged for all 374 rows.

The 39 → 1 fall is accounting, not lost evidence: all 39 rows remain in
`identified-difference.tsv` with their own verdicts, and 1 is the number this
section already reported as the true count of primary-arm non-control differences.

Gate criterion 10 `cross_region_completeness_nonvacuous` still holds on the
corrected tables (314 of 374 complete by both the flag and the re-derivation, 13
claim-family tests complete), and `n_claim_strict_conjunction_vacuous` is 0.

Regression test: `tests/test_tier1_role_and_direction.R`, 17 checks. It binds to
the shipped code rather than a transcription of it — it evaluates the real
`harvest()` body against all three real Module 07 tables, and the real
`direction_consistent`, `replicated_strict`, `strict_conjunction_vacuous` and
tier-2 filter expressions against fixtures plus the sealed run's own rows. It also
asserts that `difference_claimed` still reads `n_sensitivities_passed`, so a future
edit cannot redefine the column the gate depends on.

**This run stays unaccepted.** The fixes are in `_h/`; nothing in `_m/` was
touched, so the numbers above describe what a rerun will produce, not a sealed
result. Reacceptance needs a fresh run ID (T23), after which Modules 09, 09b and 11
rerun on it.

## Accepted runs

Machine-readable, in the schema `00_shared/gates.R::read_accepted_runs()` parses.
A run of this module spans all three regions, so `region` is the literal
`crossregion` and the per-region `vmr_set_id`s are recorded in the manifest as
`vmr_set_id_{region}` rather than in this table.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| rdg-AA-crossregion-20261008 | AA | crossregion | see manifest vmr_set_id_{caudate,dlpfc,hippocampus} | 2026-10-08 | Kynon J. Benjamin | PASS_REGION_DONOR_GENERALIZATION_QC | 10/10 gate criteria; 14 outputs; built at 03a4de5c1, git_dirty false, smoke_run FALSE. Same upstreams as rdg-AA-crossregion-20260930 plus the tier-2 donor-robust SE source cmb-AA-crossregion-20261007 (Module 05, accepted 2026-10-08). Tier 1, tier 3 and tier 4 outputs are byte-identical to the 2026-09-30 run. Tier 1: 13 of 13 prespecified claim-family tests replicate in all 3 regions, 13 strict; 374 tests, 74 primary, 13 in the claim family; both specificity controls run opposite the claim family. Tier 2: **0** primary claim-family dlpfc-minus-hippocampus differences. The meqtl_burden score_z delta -0.273 now carries the donor-robust SE 0.204 (HC3 independence form 0.075), p 0.18, q 0.36, so the 2026-09-30 run's one difference is withdrawn; 37 rows survive the statistical conjunction across all arms (8 in the primary arm, all ATAC outcomes outside the claim family), and all remain in the table. Every other tier-2 row keeps its upstream SE, which ignores the 115 shared donors and between-VMR correlation; no donor-robust refit exists for them. Tier 3: within-caudate paired delta A = 0.1432 (block-jackknife 95% CI 0.1317-0.1548), 3 replicates agree, gap_closed 0.9923; reading donor_count_is_a_plausible_major_contributor. Donor-group axis: concordance only, rho 0.8146/0.7527/0.7589 = 85.3%/88.7%/87.9% of the reliability ceiling (0.9546/0.8483/0.8635, now from lgv-all_individuals.EA-*-20260917); no ancestry effect claim. Caudate remains batch-confounded; residual excess may NOT be called biological. |

### Superseded

`rdg-AA-crossregion-20260930` (accepted 2026-09-30, superseded 2026-10-08). Its
one tier-2 difference (meqtl_burden, z -3.63, q 0.0038) rested on Module 05's
VMR-level HC3 SEs combined as if the two regions were independent. On the
donor-robust SE from `cmb-AA-crossregion-20261007` the same delta is p 0.18, so
**do not quote a DLPFC-hippocampus difference in meQTL-burden slope**. Its
reliability ceilings for DLPFC and hippocampus were computed from retired EA
runs (fixed in #127); its rho values stand. Tiers 1, 3 and 4 are byte-identical in
the replacement. `rdg-AA-crossregion-20261007`, built between the two, was never
accepted.

`rdg-AA-crossregion-20260918` (accepted 2026-09-18, superseded 2026-09-30).
Superseded on two independent counts, and its tier-1 and tier-2 counts should not
be quoted again.

First, its upstreams moved. It consumed `rra-AA-*-20260906`, `cmb-AA-*-20260825`
and `tsc-AA-*-20260902`, all three of which were replaced by the 2026-09-24/25
acceptances. The locked analysis plan is not retroactive, so the acceptance was sound when it
was made; it is the forward use that was already prohibited above.

Second, and the reason a rerun was needed rather than a pointer, its tier counts
were produced by code with four accounting defects, documented in full earlier in
this README and repaired on 2026-09-30. It reported *12 of 16 claim-family tests
replicate, 12 strict*. Three of those 16 were Module 07 arms the source module had
excluded from its own FDR family, and its strict conjunction ranged over a
caudate-only arm that could never satisfy it, so neither the numerator nor the
denominator means what it appears to.

**Its one tier-2 difference does not survive the corrected accounting, and the
accepted run's is a different test.** This is the one place where a count that
looks stable is not. Scored the way this module now scores, the superseded run has
**zero** primary claim-family differences, not one: its single claimed difference
was `expression_coupling / expression_abc / any_meqtl_support` (DLPFC 2.21 against
hippocampus 19.58, delta -17.37, q 2.4e-09), and `expression_abc` is exactly the
underpowered arm -- 250 and 243 VMRs against Module 07's locked floor of 500 --
that fix 1 excludes from the claim family. An estimate of 19.6 on 243 VMRs is the
kind of number that filter exists for.

The accepted run's one difference is `meqtl_burden` score_z, and it is new rather
than carried over. In the superseded run that same test read DLPFC 2.449 against
hippocampus 2.495, delta -0.046, p 0.64 -- not remotely significant. It reaches
DLPFC 2.094 against hippocampus 2.367, delta -0.273, p 2.9e-04, q 0.0038 only on
the `cmb-AA-*-20260924` burden estimates. So tier 2's answer moved from one
spurious difference to one real one, for two independent reasons: fix 1 removed the
artifact, and the Module 05 reacceptance created the finding.

Do not read the unfiltered counts (1 against 39) as a change of verdict either.
The test universe grew from 154 pairs to 374 when Module 04 registered its ATAC
outcomes on 2026-09-25 (T8/T11), and most of that 39 is `atac_*` rows outside the
claim family.

What did **not** change across the rerun is worth recording, because it is what
makes Module 09's decision 1 safe: tier 3 and the donor-group concordance tables
are byte-identical between `-20260925` and `-20260930`. The tier-3 reading has been
`donor_count_is_a_plausible_major_contributor` in every run, with relative
attenuation 0.1426 / 0.1432 / 0.1432.

`rdg-AA-crossregion-20260925` was never accepted. It is the run whose sealing
exposed the four defects; see the section above. Nothing in it is citable.

## Pipeline

| stage | tier | writes |
|---|---|---|
| `00_new_run.R` | — | `results/tiers.tsv`; gates all 18 region-axis (Modules 01, 02, 03, 04, 05, 07 x 3 regions) + 18 cell + 9 tier-3 upstreams, plus each tier-2 donor-robust SE source, and pins their run IDs |
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

The legacy v1 trees were retired on 2026-10-06; their tracked files are
recoverable from the annotated tag `v1-legacy-final`.

## Scope, set by the PI 2026-09-06 before implementation

This module is **not** a symmetric "shared versus context-dependent" analysis.
Brain region is perfectly confounded with sequencing batch, so replication and difference are not two equal findings
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
  for all three regions. This axis is unblocked; see the locked analysis plan.
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

The locked analysis plan already forbids raw score-level comparison across regions, so the
Module 02 score cannot be contrasted across cells at all. The region axis was
always narrower than the migration sources imply.

Do not attribute differences to ancestry-specific biology without eliminating
sample size, MAF, LD, SNP availability, assay, covariate, and brain-region
explanations. Use donor-group or population language approved by the PI:
`AA` = "Black American", `EA` = "non-Hispanic white American".

## Contract

This module follows the repository layout: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
