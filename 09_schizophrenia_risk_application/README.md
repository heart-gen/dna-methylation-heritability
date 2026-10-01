# 09_schizophrenia_risk_application — required disease application

Tests whether schizophrenia-risk variants regulate methylation within genetically anchored VMRs. Intended for the main text, conditional on surviving corrected VMRs and the new local-genetic-control axis.

**Status: accepted (AA, 2026-10-01); decision 1 resolved, decision 2 pending the stage 15 rerun.** The three per-region runs `scz-AA-{caudate,dlpfc,hippocampus}-20261001` are accepted; see **Accepted runs**, and "### Superseded" for the `-20260918` runs they replace. Every upstream is current (see **Upstream currency**), including `rdg-AA-crossregion-20260930` for the Module 08 dependency. Decision 1 is `CAUDATE_MAGNITUDE_CLAIM_NOT_SUPPORTED`. **Decision 2 lives only in `_m/combined/`, written by stage 15, which has not been rerun** -- `_m/combined/` still holds the 2026-09-19 tables built on the superseded runs, and they are not citable beside the accepted per-region numbers; the superseded answer was `RETAIN_MAIN_TEXT`. The same applies to stages 17 and 18: `config/gwas_negative_controls.yml` was PI-locked on 2026-09-20 and both stages ran against the `-20260918` runs, so their conclusions stand -- they read the axis, which is byte-identical across the rerun -- but every number in them must be regenerated before it is quoted again. Stages 17 and 18 are post hoc sensitivities that **qualify** Module 09's claim and change neither decision (AGENTS.md §7.8).

## Migrating from

`meqtl-validation/08_schizophrenia_risk_application/`, whose Phase 7 decision file records `retain_main_text_proof_of_application`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## That decision does not carry forward

The existing Phase 7 result (31 caudate loci, 361 pairs, 38 VMRs, eight
TX-coupled VMRs) was conditioned on the legacy predictability metric and
pre-repair VMR sets. Per AGENTS.md §8 it is a hypothesis to retest. Hero loci
`rs8048039` and `rs13331198` remain **candidates**: the prioritization rule in
`config/schizophrenia.yml` names no locus, and `_h/10_prioritize_loci.R`
contains none.

## Preserved design principles

PGC schizophrenia loci are defined independently of methylation results —
`_h/01_define_scz_loci.R` reads only the published intervals, the published
index-SNP table and the public summary statistics, and records
`methylation_used = FALSE`. Risk-variant×CpG tests keep their own FDR family,
corrected once over all autosomes in `_h/03b_combine_risk_variant_tests.R`.
Association is tested against the relative local SNP contribution score, not
absolute PVE or legacy predictability. GTEx eQTL evidence is support, not proof
of mediation. At most five illustrative loci are prioritized by a prespecified
rule.

The integration analysis asks whether schizophrenia-linked VMRs are enriched
for LINE/L1, H3K9me3, quiescent chromatin, high-mappability repeat intervals,
and expression/splicing coupling. Positive connects the disease application to
the repeat/repressive architecture; null presents Phase 7 as a separate proof of
disease relevance. Both are results.

## Colocalization

`config/analysis_thresholds.yml:phase7_scz` listed coloc as `deferred` while no
ancestry-matched QTL resource was wired up. PI decision 2026-09-09 authorised
it, and AGENTS.md §7.8 permits the claim once the analysis "has been run with
adequate ancestry-matched LD and passes its own gate". Two arms:

| Arm | Pair | LD | Status |
| --- | --- | --- | --- |
| `gtex_eqtl`, `gtex_sqtl` | PGC3 European GWAS × GTEx v11 European brain | matched | gate-bearing; may support a claim |
| `meqtl` | PGC3 European GWAS × this cohort's African-American CpG meQTL | **not matched** | exploratory only; never claimable |

The meQTL arm is the biologically interesting pair, but coloc assumes both
studies share an LD structure and these do not. Every row of that arm carries
`ld_ancestry_matched = FALSE` and `arm_status = CROSS_ANCESTRY_LD_UNMATCHED`,
is excluded from the gate, and cannot set `claimable_colocalization` — a
condition the gate re-checks rather than trusting. It becomes gate-eligible
without a code change once the `all_individuals` cohort carries EA donors across
all three regions: set `colocalization.arms.meqtl.qtl_ancestry: EUR` and
`gate_eligible: true`.

`coloc.abf` is primary; `coloc.susie` runs as a sensitivity check on the
single-causal-variant assumption over prioritized loci only, and never
overturns the primary result.

## Pipeline

| Stage | Script | Purpose |
| --- | --- | --- |
| 00 | `00_new_run.R` | Open the run; gate upstreams 01/02/04/05/06/07; assert one shared `vmr_set_id`; freeze the coloc arms |
| 01 | `01_define_scz_loci.R` | PGC3 intervals + index SNPs, lifted to hg38; GWAS sliced to locus windows (via `slice_gwas_sumstats.sh`) |
| 02 | `02_link_loci_to_vmrs.R` | Link loci to corrected VMRs; freeze the tested universe and the background set |
| 03 | `03_risk_variant_cpg_meqtl.py` | Per autosome: index SNPs + LD proxies × member CpGs, from Module 05's nominal pairs |
| 03b | `03b_combine_risk_variant_tests.R` | Pool the autosomes; apply the module's own FDR family once |
| 04 | `04_architecture_axis.R` | Enrichment of SCZ-linked VMRs along `local_snp_contribution_score_z` |
| 05 | `05_integration_enrichment.R` | The five prespecified repeat/repressive annotations, with per-region constraints |
| 06 | `06_transcriptional_coupling.R` | Project Module 07's accepted coupling onto loci |
| 07 | `07_gtex_support.py` | External shared-variant support from GTEx v11 brain eQTL/sQTL |
| 08 | `08_prepare_coloc_regions.py` | Per autosome: harmonise GWAS and QTL statistics onto common variants |
| 09 | `09_run_coloc.R` | Per autosome: `coloc.abf` per region, per arm |
| 09b | `09b_combine_coloc.R` | Pool coloc; roll up to loci, gate-eligible arms kept separate |
| 10 | `10_prioritize_loci.R` | Apply `ranked_composite_v1`; ≤5 loci, every ranking key recorded |
| 11 | `11_coloc_susie_sensitivity.R` | `coloc.susie` credible-set sensitivity on prioritized loci |
| 12 | `12_apply_gates.R` | Conduct gate + coloc gate; retention criteria; interpretation constraints |
| 13 | `13_plot_locus_panels.py` | Stacked regional panels for the prioritized loci |
| 14 | `14_finalize_run.R` | Seal: session info, manifest fields, checksums, read-only |
| 15 | `15_cross_region_axis_concordance.R` | **Module-level, not per-run.** Cross-region axis concordance over the three accepted per-region runs; resolves decision 2; measures donor overlap and upstream freshness. Writes `_m/combined/` |
| 16 | `16_axis_downsampling_sensitivity.R` | **Module-level.** Repeats the axis contrast in the three `AA.n118r*` caudate cells: does the contrast survive matching n to 118? Writes `_m/combined/` |

## Two independent decisions

PI 2026-09-18. This module's framing question is *does regional variation in
local genetic control of methylation intersect schizophrenia-relevant regulatory
biology, and is that relationship shared or region-dependent?* That has two
separable answers, and Module 08 tier 3 speaks to only one of them. They are
recorded as separate fields in `scz-decision.tsv` and
`_m/combined/scz-application-decisions-{cohort}.tsv`:

| decision | question | gated on | scope |
|---|---|---|---|
| `caudate_magnitude_claim` | Is caudate's *stronger* signal more than its larger donor count? | Module 08 tier 3 reading, and nothing else | caudate runs only |
| `scz_application_retention` | Is there a disorder-related pattern across all three regions? | H1 concordance (stage 15) + the per-region criteria, **excluding** `caudate_not_sample_size_artifact` | module-level |

**H1, as locked:** SCZ-linked VMRs show *lower* local genetic-control scores in
all three regions, with statistical support in at least two, including at least
one non-caudate region. A negative-but-nonsignificant third region does **not**
falsify H1.

Three constraints that the code enforces and the manuscript must respect:

- **The three regions are not independent replicates.** 100 donors appear in all
  three; DLPFC and hippocampus share 115 of 118 (Jaccard 0.96); the union is 168
  donors. Stage 15 computes this and writes it out. Report *concordance*, never a
  pooled p-value — the fixed-effect estimate stage 15 emits is labelled
  descriptive-only and its SE is anticonservative by construction.
- **Caudate may contribute to concordance but never carry it.** It is perfectly
  confounded with sequencing batch (AGENTS.md §8.1), hence
  `axis_requires_support_outside: caudate`.
- **Attenuation is not bias.** Module 08 tier 3 shows donor count is a plausible
  major contributor to caudate's larger *magnitude*. It does **not** show the
  caudate estimate is biased. Write "caudate magnitude attenuates after
  n-matching", never "n-inflated".
- **No blanket repressive-chromatin claim.** The direct linked-vs-background
  annotation test exists (stage 05), but claimability is narrower than
  significance. Stage 15 reports, per annotation, the regions where the depletion
  is claimable; cite those, and do not generalize to "repressive chromatin".

Submit one cell:

```
./_h/submit_schizophrenia.sh AA caudate      # DRY_RUN=1 to print the job graph
```

Stages 01–02 run on the submit host; 03 and 08–09 are 22-way arrays.
Colocalization runs in the `coloc` conda env (`00_shared/slurm.sh:run_r_coloc`);
`coloc` and `arrow` are not in `epigenomics`.

## Acceptance gate

`_h/12_apply_gates.R` writes `results/scz-decision.tsv`. It checks that the
analysis was **conducted** correctly, not that it produced a positive result: a
null Phase 7 is a reportable finding (AGENTS.md §7.8 provides for it
explicitly), and Modules 06 and 07 set the precedent that a reportable null
seals.

1. Locus definition recorded `methylation_used = FALSE`.
2. Every stage produced its output table.
3. The risk-variant FDR family is the one config names, corrected once.
4. At least `min_loci_tested` loci, `min_vmrs_linked` VMRs and
   `min_pairs_tested` pairs entered the analysis.
5. At least one architecture test produced a finite estimate.
6. Colocalization ran on at least `min_loci_coloc_tested` ancestry-matched loci
   with at least `min_variants_shared` harmonised variants each.
7. No claimable colocalization originates from an ancestry-unmatched arm.

Decision codes: `PASS_SCZ_APPLICATION_QC`, `PASS_SMOKE_ONLY_NOT_ACCEPTABLE`,
`FAIL_SCZ_APPLICATION_QC:<reasons>`.

### Main-text retention is reported separately

The five prespecified criteria in `config/schizophrenia.yml:retention_criteria`
are evaluated into `results/retention-criteria.tsv`, each `PASS*`, `FAIL_*`,
`NOT_APPLICABLE_NON_CAUDATE_REGION`, or `INDETERMINATE_MODULE_08`.

`caudate_not_sample_size_artifact` is read off Module 08 tier 3
(`caudate-downsampling-summary.tsv:reading`) since
`gates.require_module_08_downsampling: true` (2026-09-18), and is evaluated in
caudate runs only. It feeds `caudate_magnitude_claim` (decision 1) and is
**excluded** from `scz_application_retention` (decision 2). A per-region run
cannot resolve decision 2 — it writes `PENDING_CROSS_REGION` — so the retention
answer lives only in `_m/combined/scz-application-decisions-{cohort}.tsv`, from
stage 15.

**Current state (accepted runs `scz-AA-{region}-20261001`): decisions 1 and 2
are not symmetric, and only decision 1 is resolved.**

- decision 1 `CAUDATE_MAGNITUDE_CLAIM_NOT_SUPPORTED`, from the caudate run's own
  `retention-criteria.tsv`: tier 3 of `rdg-AA-crossregion-20260930` reads
  `donor_count_is_a_plausible_major_contributor`.
- decision 2 **pending**. It lives only in `_m/combined/`, written by stage 15,
  and stages 15 through 18 have not been rerun on the accepted runs. `_m/combined/`
  still holds the 2026-09-19 tables built on `scz-AA-*-20260918`, **which are not
  citable beside the accepted per-region numbers**. The superseded answer was
  `RETAIN_MAIN_TEXT`; the inputs stage 15 reads -- the per-region axis tables --
  are byte-identical across the rerun, so it is expected to repeat, but the
  accepted answer is whatever stage 15 writes.
- stages 17 and 18 (`step_9_negative_controls.sh`) are likewise not rerun. Every
  number in "### Negative-control traits" and "### Locus architecture" currently
  rests on the superseded runs. The trait-general conclusion and the three §7.8
  writing rules are not in doubt -- they read the axis, which did not move -- but
  the tables must be regenerated before any of their numbers are quoted again.

Note also that Module 06's accepted S-LDSC result is **null**
(`sldsc_supports_brain_enrichment = FALSE`), which
`config/analysis_thresholds.yml` lists as an `omit_or_supplement_if` condition.
The decision file carries that upstream value.

### Negative-control traits (stage 17, post hoc sensitivity)

The axis model adjusts for technical covariates only. Module 04 shows the
low-control end of the axis is gene-proximal and active, and GWAS loci of most
traits are gene-dense. So, as it stands, "SCZ-linked VMRs are low-control"
cannot be told apart from "VMRs near any GWAS locus are gene-proximal, and
gene-proximal VMRs are low-control". Stage 17 asks that question two ways.
It was requested by the PI on 2026-09-19, after the runs above were accepted,
and is recorded as a post hoc sensitivity in `config/gwas_negative_controls.yml`.
It changes neither decision.

1. **Trait distribution.** `_h/17a_extract_gwas_leads.sh` pulls the
   genome-wide significant rows from every trait in the harmonized hg38
   collection (`/projects/b1213/resources/gwas/imputed_gwas_hg38_1.1`, 114
   traits, one format) and from PGC3 schizophrenia (lifted). `_h/17_negative_control_traits.R`
   then applies **one locus rule to every trait, schizophrenia included**:
   p < 5e-8, greedy distance clumping at 500 kb, extended MHC excluded
   (the PGC3 fine-mapped table carries no MHC locus), then the Module 09
   linkage window (±500 kb) and the Module 09 logistic model. SCZ's coefficient
   is read against the other traits' distribution: overall, within category,
   and within traits with a similar lead count (×/÷ 2). `pgc.scz2` and the UKB
   self-reported schizophrenia trait are positive controls, excluded from the
   null distribution. Traits with fewer than 10 leads are fitted and shown but
   excluded from the distribution.
2. **Genomic-context adjustment.** Every contrast is fitted twice: model A with
   the Module 09 covariates, model B adding `broad_genomic_annotation`. The
   attenuation of the published-interval SCZ estimate under model B is reported.

3. **Ancestry sensitivity** (added 2026-09-19 with the results below). The
   collection is European or European-dominated and this cohort is admixed
   African American, so the config's `focal_sumstats` also carries the two
   multi-ancestry releases that include African-American cohorts: PGC3
   schizophrenia "primary" and PGC bipolar 2024 multi-ancestry. An
   ancestry-**matched** test is impossible, not skipped: the African-ancestry
   PGC3 schizophrenia and PGC bipolar releases each yield zero genome-wide
   significant variants at this rule.

Two SCZ rows exist per region: `PGC3_SCZ` (uniform rule; the comparable row)
and `PGC3_SCZ_published` (the accepted linkage; used only for the context arm).
Outputs: `_m/combined/scz-negative-control-{traits,summary,by-category,ancestry,leads}-AA.tsv`.

#### Result (2026-09-19): the depletion is trait-general

| region | SCZ, uniform rule | null median, 63 traits | SCZ percentile | traits with negative q<0.05 |
|---|---|---|---|---|
| caudate | −0.173 (p 8e-6) | −0.085 | 25th | 25 / 63 |
| dlpfc | −0.139 (p 1e-3) | −0.053 | 17th | 11 / 63 |
| hippocampus | −0.098 (p 0.02) | −0.049 | 33rd | 12 / 63 |

Schizophrenia sits at the 10th-15th percentile among the 20 traits of similar
lead count. Psychiatric traits are indistinguishable as a category (Wilcoxon
p 0.18-0.96). Anthropometric traits are the most depleted category in all three
regions; immune traits show no depletion and trend positive outside caudate.
Adding broad genomic context moves the SCZ estimate by 1-3%.

Both multi-ancestry releases give a **stronger** depletion than the European
schizophrenia file, so the ancestry mismatch is not producing the result:

| release | leads | caudate | dlpfc | hippocampus |
|---|---|---|---|---|
| PGC3 schizophrenia, European | 197 | −0.173 | −0.139 | −0.098 |
| PGC3 schizophrenia, multi-ancestry | 280 | −0.223 | −0.152 | −0.162 |
| PGC bipolar 2024, multi-ancestry | 91 | −0.191 | −0.111 | −0.127 |

**The three writing rules this imposes on the manuscript** are in AGENTS.md §7.8
and machine-readably under `interpretation.writing_rules` in the config: write
the depletion as trait-general with schizophrenia as a typical example; name the
European ancestry of the locus set wherever the locus set is described, locating
that limitation on the locus definition and not on the axis; and describe stage
17 as a qualification of Module 09 rather than a validation of it.

A parsing defect was found and fixed here on 2026-09-19. Some rows of the PGC3
"primary" release are whitespace-separated with trailing empty tab fields, so a
plain tab split left the p-value empty and awk's `"" + 0 < 5e-8` admitted every
such row as a spurious genome-wide hit. `_h/17a_extract_gwas_leads.sh` now
requires a numeric p-value, re-splits only rows that need it, and reports refused
rows. The European extraction is byte-identical before and after, so no accepted
number moved; 16 spurious leads were removed from the multi-ancestry release.

### Locus architecture (stage 18, descriptive)

Stage 17 says schizophrenia is typical. Stage 18 asks the biological question it
leaves open: what separates the traits that *are* depleted from the ones that are
not? `_h/18_locus_architecture.R` fits the same logistic model with one Module 04
annotation at a time as the predictor of trait linkage, each with and without the
control score, plus a high-mappability arm (≥0.9) for the repeat and
heterochromatin annotations as AGENTS.md §7.4 requires. Spec:
`locus_architecture` in the config. Outputs
`_m/combined/scz-locus-architecture{,-axis-link,-by-category}-AA.tsv`.

The axis contrast is largely reporting whether a trait's loci sit in active or in
quiescent sequence. Spearman correlation across the 58-62 distribution traits
between a trait's annotation enrichment and its axis estimate, score-adjusted:

| annotation | caudate | dlpfc | hippocampus | direction |
|---|---|---|---|---|
| accessible | −0.60 | −0.54 | −0.52 | marks depleted traits |
| H3K27ac | −0.48 | −0.59 | −0.44 | marks depleted traits |
| quiescent | +0.52 | +0.47 | +0.38 | marks non-depleted traits |
| H3K27me3 | +0.31 | +0.39 | +0.28 | marks non-depleted traits |
| H3K9me3 | +0.22 | +0.54 | +0.30 | marks non-depleted traits |

All of the above reach q < 0.05 in at least two regions except H3K9me3 in
caudate. Between 17 and 24 traits per region are individually enriched for
accessible chromatin at q < 0.05.

**No LINE/L1 statement is licensed.** Not one trait reaches q < 0.05 for LINE/L1
enrichment in any region or either arm, the axis-link correlation disagrees in
sign between caudate and the other two regions, and it collapses under the
high-mappability restriction. The immune category is **not** heterochromatic: its
loci are significantly *depleted* of quiescent chromatin, H3K9me3 and segmental
duplications where significant, never enriched. Immune traits are simply the
least enriched for accessible chromatin among the well-represented categories,
which is consistent with their lack of axis depletion. Note also that the
extended MHC is excluded for every trait, which removes the dominant immune
locus, so immune traits' remaining loci are a non-representative subset of their
architecture.

## Upstream currency

**Resolved 2026-10-01.** The accepted runs consume the current acceptance of
every upstream: `vmrcat-AA-*-20260816`, `lgv-AA-*-rescore-20260913`,
`rra-AA-*-20260925-a`, `cmb-AA-*-20260924`, `sldsc-AA-*-20260925`,
`tsc-AA-*-20260925-b`, and `rdg-AA-crossregion-20260930` for the Module 08
dependency that `gates.require_module_08_downsampling: true` makes mandatory.

The four-upstream supersession recorded here from 2026-09-27, and the Module 08
blocker that made the rerun impossible, are discharged. Module 08's four
tier-accounting defects were repaired and `rdg-AA-crossregion-20260930` accepted
on 2026-09-30; the Module 09 rerun followed on 2026-10-01. The prediction made
when that blocker was recorded held: neither decision moved, and the
architecture-axis tables are byte-identical across the rerun.

### Stage 18's central finding has no independent-data-source check

Stage 18's load-bearing result is that a trait's axis depletion tracks its
enrichment for **accessible chromatin and H3K27ac** (ρ −0.44 to −0.60). Both
predictors, `accessible_any` and `h3k27ac_any`, are Roadmap calls on one
reference epigenome per region — one consortium, one build, one pipeline. Module
04 registered `atac_union_frac` (BrainScope ATAC CRE union, 562,098 hg38 peaks,
an independent assay, cohort and pipeline) on 2026-09-25 for exactly this
robustness question, and `config/gwas_negative_controls.yml:locus_architecture.indicators`
does not list it. Adding it is a `pi_locked` config amendment, so it is a PI
decision under AGENTS.md §12 and is not made here. The two limitations recorded
with that track travel with it: it is the same intervals in all three regions, so
it cannot support a region-specific statement, and it is a different assay in
every respect, so read it as "does the correlation survive a change of data
source", not as a second estimate of the same quantity.

### The meQTL colocalization arm is still not gate-eligible

Its condition is that the `all_individuals` cohort carries EA donors across all
three regions. `01b_estimation_cells`, Module 02 and Module 03 each have six
accepted `all_individuals.{AA,EA}` cell runs, but **Module 05 has none** —
`05_cpg_meqtl_burden/_m/runs/` holds only `cmb-AA-*` runs, and the meQTL arm
needs an EA CpG meQTL map, not an EA score. The arm therefore stays
`CROSS_ANCESTRY_LD_UNMATCHED` and exploratory, and flipping
`colocalization.arms.meqtl.{qtl_ancestry,gate_eligible}` would assert a matching
that does not exist.

## Accepted runs

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| scz-AA-caudate-20261001 | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-10-01 | Kynon J. Benjamin | PASS_SCZ_APPLICATION_QC | Built at 644054e01, git_dirty false, smoke_run FALSE, on vmrcat-AA-caudate-20260816, lgv-AA-caudate-rescore-20260913, rra-AA-caudate-20260925-a, cmb-AA-caudate-20260924, sldsc-AA-caudate-20260925, tsc-AA-caudate-20260925-b and rdg-AA-crossregion-20260930; 106,067 output files, failures "". 2,654 VMRs linked to 522 of 612 published loci; 363,164 pairs tested in the single `scz_risk_variant_cpg_pairs_per_region` family, 27,470 significant; **94** loci with CpG-meQTL support. Axis significant, `lower_in_scz_linked`: Wilcoxon mean score_z difference -0.283 (q 8.0e-37), adjusted logistic log-odds **-0.226** per SD (SE 0.026, q 4.2e-18, adjusted for vmr_length, cpg_count, gc_content, mappability) -- **byte-identical to the superseded run**, since the axis reads Modules 01 and 02 only. 145 loci with transcriptional coupling. 70 ancestry-matched loci evaluated for colocalization, **267** claimable, coloc_claim_permitted TRUE, 19,046 cross-ancestry regions retained exploratory-only. 5 prioritized loci (rule ranked_composite_v1), 4 with GTEx support; the set is loci 73/276/432/163/198 and 3 of the superseded 5 are retained -- see "### The five prioritized loci are illustrative, and churn". Integration: 1 of 5 annotations claimable (quiescent chromatin -0.088, q 2.2e-20), **0 enriched**, all claimable annotations DEPLETED, so integration_supports_repressive_architecture_link FALSE. Decision 1 **CAUDATE_MAGNITUDE_CLAIM_NOT_SUPPORTED** -- criterion `caudate_not_sample_size_artifact` reads FAIL_SAMPLE_SIZE_ARTIFACT off rdg-AA-crossregion-20260930 tier 3 (`donor_count_is_a_plausible_major_contributor`); retention-criteria.tsv is byte-identical to the superseded run. Decision 2 PENDING_CROSS_REGION -- a per-region run cannot resolve it. Upstream `sldsc_supports_brain_enrichment = FALSE` carried from Module 06. Caudate is batch-confounded (AGENTS.md 8.1); the residual excess may not be called biological, and permitted wording is "caudate magnitude attenuates after n-matching", never "inflated" or "biased". |
| scz-AA-dlpfc-20261001 | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-10-01 | Kynon J. Benjamin | PASS_SCZ_APPLICATION_QC | Built at 644054e01, git_dirty false, smoke_run FALSE, on vmrcat-AA-dlpfc-20260816, lgv-AA-dlpfc-rescore-20260913, rra-AA-dlpfc-20260925-a, cmb-AA-dlpfc-20260924, sldsc-AA-dlpfc-20260925, tsc-AA-dlpfc-20260925-b; 81,569 output files, failures "". 2,063 VMRs linked to 565 of 612 published loci; 225,527 pairs tested, 12,547 significant; **64** loci with CpG-meQTL support. Axis significant, `lower_in_scz_linked`: Wilcoxon -0.166 (q 6.6e-11), adjusted logistic log-odds **-0.161** per SD (SE 0.026, q 1.3e-09) -- byte-identical to the superseded run. 100 loci with transcriptional coupling. 42 ancestry-matched loci evaluated, **190** claimable, 8,945 cross-ancestry exploratory. 5 prioritized loci, 4 with GTEx support; loci 490/493/544/283/136, 3 of the superseded 5 retained. **`rs13331198` (locus 544) is prioritized here at rank 3** -- the only region where either legacy hero locus passes corrected prioritization (AGENTS.md §8). Integration: 3 of 5 claimable -- quiescent chromatin -0.086 (q 6.0e-14), LINE/L1 -0.017 (q 3.2e-04), H3K9me3 -0.012 (q 0.032) -- **0 enriched**, all three DEPLETED, integration_supports_repressive_architecture_link FALSE. This is the only region where LINE/L1 is claimable, and the sign is DEPLETION in SCZ-linked VMRs -- a different quantity from Module 04's LINE/L1 enrichment along the control axis, which this does not speak to either way. Nothing here licenses a repeat statement about schizophrenia, and a single claimable region could not support one in any case. Decision 1 NOT_APPLICABLE_NON_CAUDATE_REGION. Decision 2 PENDING_CROSS_REGION. Upstream `sldsc_supports_brain_enrichment = FALSE`. |
| scz-AA-hippocampus-20261001 | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-10-01 | Kynon J. Benjamin | PASS_SCZ_APPLICATION_QC | Built at 644054e01, git_dirty false, smoke_run FALSE, on vmrcat-AA-hippocampus-20260816, lgv-AA-hippocampus-rescore-20260913, rra-AA-hippocampus-20260925-a, cmb-AA-hippocampus-20260924, sldsc-AA-hippocampus-20260925, tsc-AA-hippocampus-20260925-b; 71,159 output files, failures "". 2,017 VMRs linked to 560 of 612 published loci; 216,837 pairs tested, **10,056** significant -- **down from 12,204 in the superseded run, the only region where the meQTL evidence weakened**; 64 loci with CpG-meQTL support, unchanged. Axis significant, `lower_in_scz_linked`: Wilcoxon -0.164 (q 1.7e-10), adjusted logistic log-odds **-0.160** per SD (SE 0.027, q 2.0e-09) -- byte-identical to the superseded run. 78 loci with transcriptional coupling, up from 66. 40 ancestry-matched loci evaluated, **192** claimable, 9,509 cross-ancestry exploratory. 5 prioritized loci, **4 with GTEx support, down from 5 of 5**: loci 441/264/118/136/198, 4 of the superseded 5 retained. The new rank-1 locus 441 (rs139139) carries no GTEx support and max PP4 0.015, and entered because the Module 05 hard filter admitted it while the old rank-1 locus 164 (rs11076631) left the eligible pool; read it as an illustrative locus, not a lead. Integration: 2 of 5 claimable -- quiescent chromatin -0.057 (q 5.1e-07), H3K9me3 -0.015 (q 0.021) -- 0 enriched, both DEPLETED, integration_supports_repressive_architecture_link FALSE; LINE/L1 is not claimable here (-0.008, q 0.080). Decision 1 NOT_APPLICABLE_NON_CAUDATE_REGION. Decision 2 PENDING_CROSS_REGION. Upstream `sldsc_supports_brain_enrichment = FALSE`. |

### Superseded

`scz-AA-{caudate,dlpfc,hippocampus}-20260918` (accepted 2026-09-19, superseded
2026-10-01). Superseded because four upstreams moved, which is the condition
the section above this table described and which the rerun has now discharged.
AGENTS.md §6 is not retroactive: the 2026-09-19 acceptance was sound when it was
made.

**What is safe to carry across, because it is byte-identical in all three
regions:** the architecture-axis tests and descriptives, the PGC3 locus
definition and index SNPs, the locus-to-VMR links and the VMR linkage table.
Decision 1's `retention-criteria.tsv` is byte-identical in caudate. So neither
decision changed, and the two decisions were always the thing at risk.

**What must not be quoted from the superseded runs:** every count that reads
Module 05 or Module 07 -- pairs significant, loci with CpG-meQTL support, loci
with transcriptional coupling, ancestry-matched loci evaluated for coloc,
claimable colocalizations -- and the identity of the five prioritized loci. The
superseded caudate run also carried `git_dirty true`, so its code state is not
exactly recoverable from its commit; the accepted runs are clean.

One criterion was tested more weakly in the superseded runs than it is now.
`config/schizophrenia.yml` renamed
`enrichment_along_local_control_axis_in_primary_region` to
`..._in_every_region` in `9ba7f1aae` on 2026-09-19, after those runs sealed. The
status token `PASS_INTERPRETABLE_CONTRAST_DEPLETION` is identical in both, but it
certifies the every-region form only in the accepted runs. It passes because the
axis is significant with `lower_in_scz_linked` in 3 of 3 regions.

### The five prioritized loci are illustrative, and churn

Retained from the superseded set: caudate 3 of 5, DLPFC 3 of 5, hippocampus 4 of
5. This is expected and is recorded here so that no draft treats the five as a
stable result.

`prioritization.rule` applies a hard filter on
`has_significant_risk_variant_cpg_meqtl`, which is Module 05-dependent, so the
eligible pool changes whenever Module 05 is reaccepted. It then ranks on
`max_abs_local_snp_contribution_score_z`, and the top of that ordering is a
near-tie: hippocampus ranks 1-6 are 1.7318, 1.7170, 1.7140, 1.7136, 1.7114,
1.7053, so the five-locus cut falls inside a spread of 0.026. Both mechanisms
operated -- hippocampus's old rank-1 locus left the pool, and the surviving order
is decided at the third decimal place.

AGENTS.md §8's hero-loci retest is answered by these runs. **`rs13331198`
survives corrected prioritization in DLPFC only** (rank 3; rank 15 caudate, rank
12 hippocampus). **`rs8048039` survives in no region** (rank 48 / 53 / 47),
though it carries meQTL support, GTEx support and max PP4 0.966 in all three; it
should not be presented as a hero locus.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
