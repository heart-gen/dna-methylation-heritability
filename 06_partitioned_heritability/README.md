# 06_partitioned_heritability — genome-wide disease relevance

Tests whether common-variant heritability for brain-relevant traits concentrates
in the VMRs whose methylation is under stronger local genetic control, using
stratified LD score regression. The estimand is a **within-VMR gradient**: among
the tested VMR universe, do higher-scoring loci carry more heritability than
lower-scoring ones?

**Status: no interpretable run. The three 2026-09-08 acceptances passed QC but
answer a different question than the one above; see Estimand and Accepted runs.**
Gated on `02_local_genetic_variance` acceptance ("No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID"). Module 02 has accepted runs
for all six cells (`lgv-*-20260823`). See **Accepted runs** below.

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
| 02 | `_h/02_liftover_annotation.py` | hg38 → hg19, with every dropped interval recorded. |
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

## Negative controls

The frozen trait list carries two prespecified non-brain traits (asthma, CAD)
alongside the six brain-relevant ones. Without them "heritability for
brain-relevant traits concentrates in this annotation" is not a testable claim,
because a continuous annotation covering gene-rich, SNP-dense regions could
enrich for any polygenic trait. The controls are in the FDR family, not outside
it.

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

**Estimand withdrawn 2026-09-23; a rerun is required.** The three runs below
passed QC and are computationally sound, but they were produced by the
one-annotation model described under **Estimand** above, so their tau is not the
within-VMR gradient this module reports. `sldsc_supports_brain_enrichment = FALSE`
from these runs **must not be cited**, in Module 09's
`adds_nothing_beyond_nonsignificant_sldsc` criterion or anywhere else: a null
from an estimand that does not match the question is not a null for that
question. The runs stay as they are — `_m/` is immutable — and stages 03, 05, 06,
07 and 08 must be rerun for all three cells before any acceptance record here is
valid again. Stage 07 now fails QC on one-annotation metrics, so a rerun cannot
quietly reproduce the old estimand.

QC passed with a null scientific result: `sldsc_supports_brain_enrichment = FALSE`
in all three cells (0/8 traits FDR-significant). EUR LD scores; annotation is a
genomic feature, not a donor-group LD claim. AFR sensitivity is not part of this
acceptance.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| sldsc-AA-caudate-20260903 | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-09-08 | Kynon J.M. Benjamin | PASS_PARTITIONED_H2_QC | Null: 0 brain and 0 control traits FDR-significant; 7/8 traits with interpretable total h2 |
| sldsc-AA-dlpfc-20260903 | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-09-08 | Kynon J.M. Benjamin | PASS_PARTITIONED_H2_QC | Null enrichment; same frozen 8-trait family |
| sldsc-AA-hippocampus-20260903 | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-09-08 | Kynon J.M. Benjamin | PASS_PARTITIONED_H2_QC | Null enrichment; same frozen 8-trait family |

### Superseded runs

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
