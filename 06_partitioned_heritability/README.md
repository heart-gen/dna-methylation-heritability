# 06_partitioned_heritability — genome-wide disease relevance

Tests whether common-variant heritability for brain-relevant traits is enriched
in a **continuous** annotation built from the relative local SNP contribution
score (`local_snp_contribution_score_z`, Module 02), using stratified LD score
regression.

**Status: implemented, smoke-verified, production submitted 2026-09-02.**
Gated on `02_local_genetic_variance` acceptance (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID"). Module 02 has accepted runs
for all six cells (`lgv-*-20260823`), so this module is unblocked. No run is
accepted yet — see **Accepted runs** below.

This module depends only on Module 02. Its position in the AGENTS.md §6
dependency list is a total order, not a claim that it consumes Modules 03–05.

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

The annotation is **continuous**. Do not rebuild the retired v1 form, which
partitioned VMRs into heritable and non-heritable classes at a threshold on
`r_squared_cv` — banned by AGENTS.md §3. No threshold, no grouping, and no
absolute locus PVE enters the annotation.

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
| 01 | `_h/01_build_annotation.R` | Build the continuous hg38 annotation BED and enforce AGENTS.md §3. |
| 02 | `_h/02_liftover_annotation.py` | hg38 → hg19, with every dropped interval recorded. |
| 03 | `_h/03_make_annot.py` | Map the score onto reference SNPs (thin-annot). |
| 05 | `_h/05_compute_ldscores.sh` | Array 1–22: annotation + LD scores per chromosome. |
| 04 | `_h/04_munge_sumstats.py` | Munge the frozen trait list to LDSC format. |
| 06 | `_h/06_partition_h2.py` | S-LDSC per trait; enrichment, enrichment p, tau z. |
| 07 | `_h/07_fdr_and_gates.R` | FDR over the frozen family; acceptance gate. |
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
2. a continuous, unthresholded, ungrouped annotation carrying no absolute PVE;
3. every declared trait munged and analysed — a partial family is refused,
   because BH over a smaller family understates every q;
4. finite, positive tau standard errors for every trait;
5. at least one trait whose total observed-scale h2 is distinguishable from zero
   (`min_total_h2_z`);
6. liftover loss below `max_annotation_missing_fraction`;
7. Stage 07 decision `PASS_PARTITIONED_H2_QC`;
8. immutable Stage 09 checksums and a manual README acceptance record.

## Accepted runs

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| _(none)_ | | | | | | | |

AGENTS.md §6 makes acceptance a human step: no row appears here until a
production run's gate stage passes and the PI records it. A smoke run is never
acceptable, and a completed SLURM job is not acceptance.

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

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
