# 06_partitioned_heritability — genome-wide disease relevance

Tests whether common-variant heritability for brain-relevant traits concentrates
in the VMRs whose methylation is under stronger local genetic control, using
stratified LD score regression. The estimand is a **within-VMR gradient**: among
the tested VMR universe, do higher-scoring loci carry more heritability than
lower-scoring ones?

**Status: accepted (AA, 2026-09-25). The two-annotation estimand is implemented and
the result is an interpretable null.** The three 2026-09-08 acceptances answered a
different question and are superseded; see **Estimand** and **Accepted runs**.
Module 02's accepted runs for the AA cells are now the `lgv-AA-{region}-rescore-20260913`
rescores, which supersede `lgv-*-20260823`.

This module depends only on Module 02. Its position in 
dependency list is a total order, not a claim that it consumes Modules 03–05.

## Estimand: two annotations, not one

The S-LDSC model is baselineLD plus **two** annotations of ours, in this column
order:

| # | annotation | form | role |
|---|---|---|---|
| 1 | `VMR_TESTED` | binary | membership in the tested VMR universe |
| 2 | `LOCAL_SNP_CONTRIBUTION_Z` | continuous | the within-cell relative score |

`local_snp_contribution_score_z` is a standardized within-cell midrank
percentile, so it is symmetric about zero **by construction** and zero is the
value of a median-ranked VMR. A single-column annotation therefore gave a SNP in
no VMR and a SNP in a median-ranked VMR the same number, and with no membership
term in the model there was nothing to absorb a VMR-versus-genome difference:
tau blurred the within-VMR gradient with a membership effect it never meant to
test. That is why the one-annotation result was not merely negative but
uninterpretable, and it is why a null from it is not reportable.

With `VMR_TESTED` in the model, tau on the score is conditional on membership
and estimates the within-VMR gradient. tau on `VMR_TESTED` answers a second,
separately interesting question — are tested VMRs enriched at all — which the
one-annotation model could not ask either.

**Prespecified, before any two-annotation result was seen:** the FDR family stays
exactly the frozen traits on the **score** annotation, which carries the
hypothesis. The membership tau is written to `sldsc-membership-metrics.tsv` with
a nominal p and no q-value, so the family size does not change and no q-value
the primary hypothesis depends on is revised.

The contract lives in `_h/annotations.py` — names, column order, roles, and the
assertion that neither name is a prefix of the other, since stage 06 identifies
`.results` rows by name. Column order is positional in `.annot.gz`,
`.l2.ldscore.gz` and `.l2.M_5_50`, so every stage reads it from that one file;
row identification downstream is always by name and never by position.

Enrichment is reported for both annotations and is interpretable for **only** the
binary one. For a signed continuous annotation LDSC's denominator is a signed
sum, so the ratio is not a share; each emitted row carries
`enrichment_interpretable` rather than leaving that to a reader's memory. The
`m_5_50` column makes it visible: a SNP count for `VMR_TESTED`, a signed sum for
the score.

See `writing-notes/ISSUE_06_two_annotation_model.md` for the full argument.

### Not fixed by this, and still open

The annotation covers ~0.6% of SNPs and the module has **no positive control**, so
a null cannot be distinguished from a power null whatever the estimand. A
separate issue should run the same pipeline on an annotation of comparable
footprint known to be enriched for brain traits (Roadmap brain DHS or H3K4me3).
Correcting the estimand does not remove that limitation.

**Being addressed (2026-10-07):** stages 10-13 run a published brain WGBS
annotation of comparable footprint through this pipeline. See **Positive
control: external annotations** below.

## Why this module exists

`config/analysis_thresholds.yml` names `adds_nothing_beyond_nonsignificant_sldsc`
as an `omit_or_supplement_if` criterion for the schizophrenia application in
Module 09. Until this module produces a result, that criterion cannot be
evaluated. PI decision 2026-08-26: build S-LDSC as its own module rather than a
Module 09 sub-analysis, because it is SNP-level and genome-wide where every
other module is VMR-level.

## Migrating from

`local-snp-prediction/BA_only/tissue_comparison/clinical_enrichment/s-ldsc/`
(and the `all_individuals` counterpart), including `make_annot_continuous.py`,
`ldsc_wrapper.py`, `region_heritability.py`, `fdr_correction.py`, and
`interpreting_sldsc_results.md`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Endpoint discipline

The score annotation is **continuous**. Do not rebuild the retired v1 form, which
partitioned VMRs into heritable and non-heritable classes at a threshold on
`r_squared_cv` — banned. No threshold, no grouping, and no
absolute locus PVE enters the annotation.

`VMR_TESTED` is binary, and that is not a reintroduction of the banned form: it
indicates **membership in the tested universe**, not a class of local genetic
control. It is derived from the same overlap test as the score and from no
threshold on it. `_h/01_build_annotation.R` still fails closed if the score
itself looks like a partition rather than a gradient.

Report the metrics named in the legacy `interpreting_sldsc_results.md`:
enrichment, enrichment p, and the tau coefficient z-score. A significant
enrichment on a continuous annotation is a statement about where common-variant
heritability concentrates, not about any individual VMR.

## LD reference and ancestry

**Resolved 2026-09-02 (PI):** primary arm uses the standard EUR-derived LD
scores; an AFR LD-score sensitivity arm is built separately.

Standard S-LDSC LD scores are EUR-derived, while the primary cohort is Black
American donors and the annotation is defined on VMRs discovered in that cohort.
v1 carried separate `_AA` and `_EA` variants of `region_heritability.py`.

The position this module takes, which **must appear in Methods and not be left
implicit**:

> The annotation is a genomic feature — a map of where in the genome
> methylation variance is under strong local genetic control — rather than a
> property of the donors in whom it was measured. Partitioning EUR-ancestry GWAS
> heritability with EUR LD scores over that annotation therefore asks whether
> common-variant heritability concentrates in those genomic regions. It does not
> assume the donor groups share LD structure, and no statement about ancestry-
> specific genetic architecture follows from it.

`ld_reference_arm` in `config/partitioned_heritability.yml` selects the arm, so
the sensitivity analysis is a config switch rather than a forked pipeline.

**This closes task T15** (the S-LDSC LD reference for an AA primary arm, left open
by the 2026-08-26 gap review). The decision is the PI's, dated, recorded in
versioned config as `ld_reference_arm: eur_primary`, and stamped onto every run's
`partitioned-h2-decision.tsv`. What remains is not a decision but a writing
obligation: the indented paragraph above must appear in Methods verbatim, and it
is not yet in `content/`.

### AFR sensitivity arm (separate issue, not required for the primary result)

No AFR LD scores exist on Quest: `/projects/b1213/resources/ldsc/` holds EUR
panels only, and `1000G_Phase3_plinkfiles.tgz` contains EUR. They must be built
from `/projects/b1213/resources/1kGP/GRCh38_phased_vcf/` (subset AFR samples,
convert with plink2, lift to hg19 because this pipeline is hg19), then
`ldsc.py --l2`. Ancestry-matched sumstats exist for few of the frozen traits, so
the sensitivity arm will cover a reduced trait list and must be reported as a
targeted robustness check rather than a parallel screen.

## Pipeline

| Stage | Script | Purpose |
|---|---|---|
| 00 | `_h/00_new_run.R` | Mint the run ID; gate on the accepted Module 02 run; freeze the trait family. |
| 01 | `_h/01_build_annotation.R` | Build the continuous hg38 annotation BED. |
| 02 | `_h/02_liftover_annotation.py` | hg38 → hg19, with every dropped interval recorded; drops a VMR whose lifted width changes more than `liftover_max_span_ratio` (2×). |
| 03 | `_h/03_make_annot.py` | Map membership + score onto reference SNPs (thin-annot, two columns). |
| 05 | `_h/05_compute_ldscores.sh` | Array 1–22: annotation + LD scores per chromosome; `05a` checks the annotation set that reached them. |
| 04 | `_h/04_munge_sumstats.py` | Munge the frozen trait list to LDSC format. |
| 06 | `_h/06_partition_h2.py` | S-LDSC per trait; one metrics row per annotation. |
| 07 | `_h/07_fdr_and_gates.R` | FDR over the frozen family (score annotation only); membership reported separately; acceptance gate. |
| 08 | `_h/08_plot.py` | Figures. |
| 09 | `_h/09_finalize_run.R` | Seal the run (sealing is not acceptance). |

Submit one cell with
`_h/submit_partitioned_heritability.sh <cohort> <region>`; `DRY_RUN=1` prints
the job graph, `SMOKE_N=1` permits unaccepted upstreams, and `SMOKE_CHROMS`
restricts the LD-score array.

## Liftover span guard (added 2026-10-08)

Stage 02 lifts each VMR's two ends separately with pyliftover. If the ends land
in different chain blocks, the "lifted" interval covers everything between
them. Nothing checked for this until 2026-10-08, and the accepted 2026-09-25
cells carry these intervals:

| region | hg38 VMR | hg38 width | hg19 width |
|---|---|---:|---:|
| caudate | chr8:144270354-144271615 | 1,261 | 170,346 |
| caudate | chr8:141737779-141738438 | 659 | 71,818 |
| DLPFC | chr8:141737835-141738448 | 613 | 71,772 |
| DLPFC | chr1:228556705-228556858 | 153 | 4,615 |
| DLPFC | chr14:106345434-106350667 | 5,233 | 19,472 |
| hippocampus | chr1:148679673-148679757 | 84 | **24,596,960** |
| hippocampus | chr8:141737805-141738448 | 643 | 71,802 |
| hippocampus | chr1:228556520-228556858 | 338 | 4,800 |

The hippocampus chr1 VMR lifts across the centromere to chr1:120.6-145.2 Mb,
which holds 3,978 reference SNPs. All of them entered both of that cell's
annotations, the membership one and the score one.

Stage 02 now drops a VMR whose hg19 width is more than
`liftover_max_span_ratio` (2) times larger or smaller than its hg38 width. It
records each one in `excluded/liftover-span-changed.tsv` and counts them in
`liftover-summary.tsv`.
- **Where 2 sits:** in the accepted cells the largest within-block change is
  1.49× (an indel-sized 493 bp), and the smallest cross-block change is 3.72×.
- **Smoke test:** `tests/test-liftover-span-guard.py` reruns stage 02 on each
  accepted cell's hg38 annotation. Exactly the listed VMRs drop, and every other
  hg19 interval is identical.
- **Module 04 is not affected.** Its rtracklayer liftover keeps only VMRs that
  lift to a single interval.

### Rerun with the guard (2026-10-08, accepted 2026-10-08)

`sldsc-AA-{caudate,dlpfc,hippocampus}-20261008` were sealed 2026-10-08 at
`8a8630d2a`, with `git_dirty false`. All three return `PASS_PARTITIONED_H2_QC`,
with 8 of 8 traits completed and the same upstream `lgv-AA-{region}-rescore-20260913`.
The PI accepted them on 2026-10-08, replacing the 2026-09-25 cells (see Accepted
runs and Superseded runs below).

| region | VMRs dropped by the guard | membership SNPs (MAF ≥ 5%), before → after | brain FDR hits | control FDR hits |
|---|---:|---:|---:|---:|
| caudate | 2 | 41,875 → 41,871 | 0 / 6 | 0 / 2 |
| DLPFC | 3 | 30,239 → 30,216 | 0 / 6 | 0 / 2 |
| hippocampus | 3 | 31,822 → 30,511 | 0 / 6 | **1 / 2 (CAD)** |

- **Caudate and DLPFC barely move.** No trait's tau z changes by more than 0.03.
- **Hippocampus changes for every trait.** The single 24.6 Mb interval made up
  4.1% of the membership annotation's reference SNPs. With it removed:
  - SCZ tau z goes from 1.45 to 0.82;
  - asthma goes from −2.67 (q 0.06) to −1.00;
  - CAD goes from 0.18 to 3.05 (tau p 0.0023, q 0.018).

  So the 2026-09-25 hippocampus tau values are artifacts of one mis-lifted
  interval, not estimates with noise.
- **The module result is unchanged.** `sldsc_supports_brain_enrichment = FALSE`
  in every cell, and no brain trait reaches FDR anywhere.
- **The CAD result is a prespecified non-brain control reaching FDR in one cell.**
  It is not repeated in caudate (tau z 0.72) or DLPFC (0.81). The gate does not
  act on a control result, and this one is reported, not explained away. It
  weakens one reading in particular: a future brain-trait hit on this
  annotation could not be called brain-specific without beating the controls in
  the same cell.

## Negative controls

The frozen trait list carries two prespecified non-brain traits (asthma, CAD)
alongside the six brain-relevant ones. Without them "heritability for
brain-relevant traits concentrates in this annotation" is not a testable claim,
because a continuous annotation covering gene-rich, SNP-dense regions could
enrich for any polygenic trait. The controls are in the FDR family, not outside
it.

## Positive control: external annotations (stages 10-13, non-gating)

**Why.** The accepted null (`sldsc_supports_brain_enrichment = FALSE`) is "no
detectable enrichment at this footprint", and the module had no positive
control to separate it from a power null. The PI asked on 2026-10-07 for the
published neuronal CG-DMRs to be run "through exactly the same S-LDSC pipeline,
baseline model, reference LD, and GWAS summary statistics as the
genetic-control VMR annotation". DLPFC neuronal H3K27ac was named as an
optional secondary control. Its peaks (Girdhar 2018, Synapse syn9998643) are
under NRGR controlled access, so the PI dropped it on 2026-10-08.

**Annotations** (`config/sldsc_external_annotations.yml`, `pi_locked`):

| name | source | intervals | role |
|---|---|---:|---|
| `RIZZARDI_CGDMR_NEUNPOS_REGIONS` | Rizzardi 2019 `CG-DMRs.pos.bb`: CG-DMRs between brain regions in NeuN+ nuclei, WGBS | 13,074 (11.9 Mb) | primary positive control |
| `RIZZARDI_CGDMR_NEUNPOS_VS_NEG` | Rizzardi 2019 `CG-DMRs.pos_vs_neg.bb`: NeuN+ vs NeuN− CG-DMRs | 100,875 (70.0 Mb) | secondary |
| `VMR_TESTED_<REGION>` | the accepted run's own `annotation-hg19.bed` | per region | like-for-like comparator |

`inputs/supportfiles/_h/03_build_rizzardi_dmr_asset.py` downloads each bigBed,
checks its SHA-256 against the config, converts it to a merged autosomal hg19
BED, and checks the interval count.

**What is identical, and what is not.**
- **Identical, read rather than restated:** baselineLD v2.2, the LD reference,
  weights, frq files, print-snps and `ld_wind_cm`. They all come from
  `config/partitioned_heritability.yml` at the accepted `ld_reference_arm`.
- **Identical GWAS:** the 8 frozen traits. Stage 10 proves the three accepted
  cells' munged `.sumstats.gz` files have identical decompressed content, then
  copies them.
- **Different, by design:** each external annotation enters alone, as one
  binary annotation on top of baselineLD. The accepted VMR model carries two
  annotations (membership and score). So each region's VMR membership is
  refitted alone here as `VMR_TESTED_<REGION>`, the like-for-like row. The
  accepted two-annotation rows sit beside it unchanged.

**Reported per annotation × trait:** proportion of SNPs, proportion of h2,
enrichment, its SE and p, τ with SE and two-sided p, and τ\* (Gazal 2017:
τ·sd(a)·M/h2g over reference SNPs with MAF ≥ 5%) with SE. BH q is computed
across the 8 traits within each annotation, for description only.

| Stage | Script | Purpose |
|---|---|---|
| 10 | `_h/10_external_new_run.R` | Open `sldsc-{cohort}-external-{date}`; check the accepted runs and sumstats identity; stage the annotation BEDs. |
| 11 | `_h/11_external_ldscores.sh` (+ `11a_external_make_annot.py`) | Array 1–22: thin annot and LD scores for each annotation. |
| 12 | `_h/12_external_partition_h2.py` | One S-LDSC regression per annotation × trait; metrics including τ\*. |
| 13 | `_h/13_external_summarize.R` | Refuse a partial grid; add the accepted rows; descriptive BH; seal. |

Submit with `_h/submit_external_annotations.sh <cohort>`. `SMOKE_N=1`,
`SMOKE_CHROMS` and `SMOKE_TRAITS` restrict a smoke run, and `DRY_RUN=1` prints
the job graph.

**Reading, fixed before the production run:**
- If a Rizzardi annotation is enriched for the brain traits and the standalone
  VMR annotation is not, the pipeline can detect enrichment at this footprint,
  and the VMR null is not a pure power null.
- If neither is enriched, the VMR null stays "no detectable enrichment at this
  footprint".

Either way the accepted decision is unchanged. Every row carries
`gating = FALSE`.

**Result: `sldsc-AA-external-20261008`** (sealed 2026-10-08 16:17 at `c0daa121a`,
clean tree, 42/42 jobs, 40 regressions). It reads the span-guard cells
`sldsc-AA-{caudate,dlpfc,hippocampus}-20261008`. Munged sumstats were identical
across the three cells for all 8 traits. The staged Rizzardi BEDs are
byte-identical to the asset content and to the 2026-10-07 run's, and no
annotation is listed as not covered.

| annotation | % common SNPs | SCZ enrichment (SE; p) | SCZ τ\* (SE) | BIP enrichment (p) | smoking enrichment (p) | asthma | CAD |
|---|---:|---|---|---|---|---|---|
| NeuN+ between-region CG-DMRs | 0.41 | **16.3 (2.85; 3.4e-07)** | 0.894 (0.184) | **15.3 (5.4e-04)** | **7.6 (3.2e-05)** | −13.2 (0.23) | −6.6 (0.026) |
| NeuN+ vs NeuN− CG-DMRs | 2.45 | **6.23 (0.88; 2.1e-08)** | 0.665 (0.146) | **6.72 (2.8e-09)** | **6.11 (3.8e-08)** | 0.28 (0.87) | 2.29 (0.45) |
| VMR membership, caudate (standalone) | 0.70 | 1.25 (0.82; 0.76) | 0.097 (0.071) | 0.43 (0.59) | 0.29 (0.30) | −1.91 (0.58) | 0.80 (0.91) |
| VMR membership, DLPFC (standalone) | 0.51 | 0.01 (0.94; 0.29) | −0.015 (0.069) | −0.92 (0.11) | 1.00 (1.00) | −8.2 (0.12) | 0.33 (0.76) |
| VMR membership, hippocampus (standalone) | 0.51 | 0.13 (0.98; 0.38) | −0.010 (0.072) | −1.10 (0.091) | 1.31 (0.73) | −8.7 (0.11) | −0.07 (0.60) |

Every metric is in `results/external-annotation-metrics.tsv`: proportion of
SNPs, proportion of h2 and its SE, enrichment with SE and p, τ with SE and
two-sided p, τ\* with SE, and descriptive BH q. The accepted two-annotation VMR
rows are in `results/vmr-membership-accepted-metrics.tsv`. They agree with the
standalone ones; SCZ is 1.75 / −0.06 / −0.21.

**Reading.**
- **The pipeline detects enrichment at this footprint.** The NeuN+
  between-region CG-DMRs cover 0.41% of common SNPs, less than any VMR
  annotation. They carry 6.6% of SCZ h2 (16-fold) and 6.3% of BIP h2. Both
  q-values are ≤ 1.5e-3, and the result matches the original report. AD, PD
  and MDD are not enriched; for MDD the SE is too wide to say anything.
  Neither non-brain control is positively enriched.
- **So the VMR null is not a pure power null.** The standalone VMR SCZ
  enrichment's upper 95% bound is 2.9 / 1.8 / 2.0. That excludes anything near
  the neuronal CG-DMRs' 16-fold, and even their 6-fold NeuN+ vs NeuN− level.
  Write: "no detectable enrichment, and enrichment of the magnitude seen for
  neuronal CG-DMRs is excluded". Do not write "no enrichment".
- No VMR row is FDR-significant for any trait; the smallest within-annotation
  q is 0.45.
- The comparison is between annotations, not a test of a difference. The
  Rizzardi DMRs were defined from neuronal versus glial and between-region
  contrasts, a different selection rule from population variability.
- The accepted decision (`sldsc_supports_brain_enrichment = FALSE`) is
  unchanged.

**What changed from `sldsc-AA-external-20261007`.** That run read the
2026-09-25 cells, whose liftover had no span check:
- **The worst case was hippocampus `chr1:148679673-148679757`.** In hg38 it is
  84 bp; `_h/02_liftover_annotation.py`, lifting start and end separately, put
  it at chr1:120,612,168-145,209,128 in hg19, which is 24.6 Mb across the
  centromere and 3,978 reference SNPs. The others were chr8 VMRs of about
  650 bp that became 72-170 kb.
- `8a8630d2a` added `liftover_max_span_ratio = 2.0`, which drops and records
  such intervals. The three cells were rerun as `-20261008`.
- The cells dropped 2 / 3 / 3 VMRs, which removes 2 / 1 / 2 intervals from the
  merged annotations staged here. Hippocampus goes from
  31.3 Mb to 6.7 Mb, caudate from 8.66 to 8.42 Mb and DLPFC from 6.63 to
  6.53 Mb.
- The Rizzardi rows are identical to the last digit. The caudate and DLPFC VMR
  rows move in the third significant figure. Hippocampus moves more, as
  expected. Its SCZ enrichment goes from −0.61 to 0.13, and its nominal BIP
  p 0.011 becomes 0.091. Asthma goes from 11.7 to −8.7, both with wide SEs.
  No reading changes.
- The 2026-10-07 run also listed `PSYCHENCODE_NEUNPOS_H3K27AC_DLPFC` as not
  covered. That control has since been dropped from the config.

## Acceptance gate

For each cohort-by-region cell, acceptance requires:

1. an accepted Module 02 run for the same cell, and its `vmr_set_id` recorded;
2. a continuous, unthresholded, ungrouped score annotation carrying no absolute
   PVE;
3. every declared trait munged and analysed — a partial family is refused,
   because BH over a smaller family understates every q;
3b. both annotations present for every trait, so the reported tau is conditional
   on VMR membership. A one-annotation run fails QC rather than earning a
   caveat: its tau is not the within-VMR gradient this module reports, null or
   otherwise;
4. finite, positive tau standard errors for every trait on **both**
   annotations — they are fitted jointly, so a degenerate SE on either means the
   joint fit did not identify the model;
5. at least one trait whose total observed-scale h2 is distinguishable from zero
   (`min_total_h2_z`);
6. liftover loss below `max_annotation_missing_fraction`;
7. Stage 07 decision `PASS_PARTITIONED_H2_QC`;
8. immutable Stage 09 checksums and a manual README acceptance record.

## Accepted runs

`tau_conditional_on_vmr_membership = TRUE` in all three cells, so the reported tau
**is** the within-VMR gradient this module defines as its estimand. These cells
replace the 2026-09-25 acceptances, whose hippocampus annotation carried one
liftover-inflated interval (see "Liftover span guard").

The scientific result for brain traits is null, which is a legitimate outcome:
0 of 6 brain traits are FDR-significant on the score annotation in any region, and
`sldsc_supports_brain_enrichment = FALSE` everywhere. One prespecified non-brain
control is significant: CAD in hippocampus (tau z 3.05, q 0.018). It does not repeat
in caudate (0.72) or DLPFC (0.81). It is reported, not explained away, and it means
a future brain-trait hit on this annotation could not be called brain-specific
without beating the controls in the same cell. All 8 declared traits completed
Stage 06, and 7 of 8 have total observed-scale h2 distinguishable from zero. EUR LD
scores (`eur_primary`). The AFR sensitivity arm is not part of this acceptance.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| sldsc-AA-caudate-20261008 | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-10-08 | Kynon J.M. Benjamin | PASS_PARTITIONED_H2_QC | Span-guard rerun: 2 VMRs dropped (hg19 span changed more than 2x), membership SNPs 41,875 to 41,871. 0/6 brain and 0/2 control traits FDR-significant; every tau z within 0.03 of the 2026-09-25 cell. 8/8 traits completed, 7/8 with interpretable total h2; tau conditional on VMR membership |
| sldsc-AA-dlpfc-20261008 | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-10-08 | Kynon J.M. Benjamin | PASS_PARTITIONED_H2_QC | Span-guard rerun: 3 VMRs dropped, membership SNPs 30,239 to 30,216. 0/6 brain and 0/2 control traits FDR-significant; every tau z within 0.03 of the 2026-09-25 cell; tau conditional on VMR membership |
| sldsc-AA-hippocampus-20261008 | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-10-08 | Kynon J.M. Benjamin | PASS_PARTITIONED_H2_QC | Span-guard rerun: 3 VMRs dropped, including the 84 bp VMR lifted to 24.6 Mb; membership SNPs 31,822 to 30,511 (4.1%). 0/6 brain FDR-significant; 1/2 control: CAD tau z 3.05, q 0.018. Every trait's tau moved from the 2026-09-25 cell (SCZ z 1.45 to 0.82), which was an artifact of the one interval; tau conditional on VMR membership |

Provenance: `lgv-AA-{region}-rescore-20260913` -> this run, sealed
2026-10-08T07:54-07:56 at commit `8a8630d2a`, `git_dirty = false`,
`smoke_run = FALSE`, 102/102 SLURM jobs completed,
`absolute_pve_interpretation_allowed = FALSE`, `fdr_family = traits_within_cell`
over `LOCAL_SNP_CONTRIBUTION_Z` only, with the membership tau reported separately
and carrying no q-value. `liftover_max_span_ratio = 2.0`.

### What this null does and does not license

It **does** let Module 09's `adds_nothing_beyond_nonsignificant_sldsc` criterion
in `config/analysis_thresholds.yml` finally be evaluated against a matching
estimand, which the 2026-09-08 runs could not supply (T17). Note that the
criterion is `pi_judgment_not_gated`: no Module 09 decision is computed from it,
so this unblocks a judgement the PI can now make, not a gate that was failing.

It does **not** establish that common-variant heritability is unenriched in these
VMRs. The limitation recorded under **Not fixed by this, and still open** is
unchanged by the estimand repair: the annotation covers ~0.6% of SNPs and the
module still has **no positive control**, so a null cannot be distinguished from a
power null. Write it as "no detectable enrichment at this footprint", never as
"no enrichment". A separate issue should run the same pipeline on an annotation of
comparable footprint known to be enriched for brain traits.

### The estimand is now declared, not only implemented (2026-09-27)

The two-annotation set was a **code-level** contract in `_h/annotations.py`: the
names, their order, their roles, which one carries the primary hypothesis and
which one's enrichment is interpretable. Every stage read it from that file, and
the file is snapshotted into `runs/{RUN_ID}/code/_h`, so each run recorded the set
it used — but a PI reading `config/` could not see any of it.

The PI declared it in `config/partitioned_heritability.yml` on 2026-09-27, under
`annotation.annotations`, `annotation.primary_hypothesis_annotation` and
`annotation.fdr_family_annotation`. A declaration that nothing verifies would be
worse than none, so it is enforced in two places:
`_h/annotations.py::check_against_config()`, called by stage 03 before any
annotation is written, and a preflight block in `_h/00_new_run.R`, which stops a
drifted config before a run directory exists. Both refuse a reordering (the order
is positional in the `.annot.gz` and the LD-score columns), a changed role, a
changed `enrichment_interpretable`, and any move of the FDR family off
`LOCAL_SNP_CONTRIBUTION_Z` — that last would revise already-computed q-values,
which AGENTS.md §10.3 forbids.

Both checks **tolerate a config with the keys absent**, so the accepted
2026-09-25 runs stay reproducible from their own snapshotted config, which
predates the declaration. Those runs executed exactly what is now declared, and
they say so themselves: `partitioned-h2-decision.tsv` stamps
`annotations_in_model = LOCAL_SNP_CONTRIBUTION_Z,VMR_TESTED` and
`fdr_family_annotation = LOCAL_SNP_CONTRIBUTION_Z` from the observed metrics. The
edit does change the file's SHA-256, so the
`config_partitioned_heritability_sha256` in those three manifests now predates
the live config; no rerun is required for that reason alone, because no key the
model reads changed.

### Superseded runs

`sldsc-AA-{caudate,dlpfc,hippocampus}-20260925` (accepted 2026-09-25, superseded
2026-10-08). Stage 02 lifted each VMR's two ends separately and never checked the
lifted span, so 2-3 VMRs per cell became much wider in hg19. The worst was the
hippocampus VMR chr1:148679673-148679757 (84 bp), which became 24.6 Mb across the
centromere and put 3,978 reference SNPs into both annotations. Caudate and DLPFC
are numerically unchanged by the fix. Every hippocampus tau moved, so **do not cite
the 2026-09-25 hippocampus tau values**. The module decision
(`sldsc_supports_brain_enrichment = FALSE`) is the same in both sets.

`sldsc-AA-{caudate,dlpfc,hippocampus}-20260903` (accepted 2026-09-08, withdrawn
2026-09-23). Computationally sound, but produced by the **one-annotation** model:
with no `VMR_TESTED` term there was nothing to absorb a VMR-versus-genome
difference, so tau blurred the within-VMR gradient with a membership effect it
never meant to test. `sldsc_supports_brain_enrichment = FALSE` from these runs
must not be cited anywhere, including Module 09's
`adds_nothing_beyond_nonsignificant_sldsc` criterion: a null from an estimand that
does not match the question is not a null for that question. Stage 07 now fails QC
on one-annotation metrics, so a rerun cannot quietly reproduce the old estimand.
Superseded by `sldsc-AA-{region}-20260925`.

`sldsc-AA-{caudate,dlpfc,hippocampus}-20260902` failed at Stage 06 for all eight
traits and must not be used. The LD scores and the S-LDSC regressions themselves
were sound; the defect was in labelling. This LDSC build writes the LD-score
column of a single `--thin-annot` file as a bare `L2`, discarding the annotation
name carried in the `.annot.gz` header, so the annotation reached the `.results`
file as the positional label `L2_1` and Stage 06 refused to identify it by
position. Stage 05 now restores the name (`05a_label_ldscore_column.py`), which
makes each `.results` file self-describing and keeps no stage dependent on row
order. Because Stage 05 changed, those runs were replaced rather than resumed:
their code snapshots no longer describe the code that would produce them.
Superseded by `sldsc-AA-{region}-20260902-a`.

`sldsc-AA-{caudate,dlpfc,hippocampus}-20260902-a` confirmed that fix -- all eight
traits completed Stage 06 in all three cells and the `.results` rows are named
`LOCAL_SNP_CONTRIBUTION_ZL2_1` -- but then failed at Stage 07. `config` declares
`fdr_method: fdr_bh`, the name of the procedure, and Stage 07 passed that string
straight to `p.adjust()`, which knows the correction as `BH` and aborts on any
other spelling. Stage 07 now translates config's procedure name to `p.adjust()`'s
argument and aborts on an unrecognized one; the config keeps naming the procedure
rather than the R argument, and `sldsc-metrics.tsv` still records `fdr_bh`. No
S-LDSC output was wrong -- Stage 07 never produced any -- but Stage 07 changed,
so these runs are likewise replaced rather than resumed.
Superseded by `sldsc-AA-{region}-20260903`.

## Contract

This module follows: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
