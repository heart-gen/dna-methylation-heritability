# 10_integrated_manuscript_outputs — one source of truth

Consumes only accepted immutable upstream runs and produces every manuscript number, table, and figure.

**Status: Figures 1 and 2, Table 1, and the cohort QC panels implemented.**
Figures 1-2 last built as `fig-all-20260826-a`; Table 1 and QC are new and have
no accepted run yet. Figures 3-5,
the manuscript-number registry, and the consolidated tables remain gated on
acceptance of their upstream modules (AGENTS.md §6: "No downstream production
run may consume an upstream result until the upstream README records a passing
acceptance gate and immutable run ID"). Modules 03 and 05 have production runs
but no recorded acceptance gate; Module 04 has smoke runs only, so Figure 3
cannot be built.

## Implemented figures

| Output | Content | Upstream runs |
|---|---|---|
| `figure1_vmr_catalog[_all_individuals]` | cohort, catalog structure, off-array coverage vs 450K, genic context, distance to nearest gene, legacy turnover | `vmrcat-*-20260816`, `vmrcatqc-*-20260826-a` |
| `figure1_vmr_catalog[_all_individuals]_epic` | as above, EPIC as the stricter array comparator | same |
| `figure2_local_genetic_control[_all_individuals]` | estimator concordance vs locus geometry, held-out R² across rank deciles, cross-region rank concordance, genic context across the rank, denominators | `lgv-*-20260823`, `vmrcatqc-*-20260826-a` |
| `figureS_local_control_audit_unbounded[_all_individuals]` | unbounded joint estimate distribution, **audit only** | `lgv-*-20260823` |
| `table1_cohort` (`.tsv`, `.tex`) | donor demographics, both arms x three regions | `vmrcat-*-20260816` |
| `figureS_ancestry_pcs` | genotype PC1/2 with 1000 Genomes reference, donor-group confirmation | `vmrcat-*-20260816` |
| `figureS_sample_integrity` | cross-region donor concordance, sample-swap screen | `vmrcat-*-20260816` |

AA is the primary arm; `all_individuals` renders from the same builders as the
sensitivity supplement.

### Figure 2 constraint

Module 02's terminal decision is
`PASS_RELATIVE_GENETIC_CONTROL_FAIL_ABSOLUTE_LOCUS_PVE`, so only the relative
score is admissible: no absolute PVE, no thresholds, no heritable/non-heritable
groups. `02_figure2_local_control.R` asserts this at runtime and stops if a
retired column reappears upstream.

`local_snp_contribution_score` is a within-cell midrank percentile, so its
distribution is uniform by construction. The figure therefore shows what the
ranking *agrees with* -- independent estimators, held-out prediction, the other
regions -- rather than the distribution of the score itself.

## Table 1 and cohort QC (PI decision D3, 2026-08-26)

`04_table1_cohort.R` and `05_qc_sample_integrity.R` replace `sample_summary/`
and the parts of `qc_analysis/` that Module 01 does not cover.

The legacy Table 1 is not merely unmigrated, it reports the **wrong cohort**: it
derived its donor set from `vmr-analysis/all_individuals/{region}/_m/samples.txt`
(invalidated by V1) and honoured the retired sample blacklists. The replacement
reads the donor list from the accepted Module 01 run and hard-stops if the count
disagrees with the locked `design_n` in `config/cohorts.yml`.

### Open finding: the retired blacklist tracks a real QC signal

`config/cohorts.yml` states the legacy blacklists "were never a QC exclusion"
and existed only to reconcile a stale phenotype file. The cross-region
integrity screen does not support that reading. Of the 8 donors readmitted by
retiring the blacklists, **4 are flagged** by the screen (Br1249, Br1693,
Br1700, Br1927), against 11.7% of all tested donors -- Fisher exact
p = 0.0064, OR = 9.8.

This is a post-hoc test on an admittedly underpowered screen and is **not**
grounds for reinstating the blacklist. It is grounds for the PI to look at
those four donors before submission, because the current v2 position is that
their exclusion was purely clerical.

## Build

Figures 1 and 2 depend on support files and a Module 01 QC refresh that the
accepted `vmrcat-*-20260816` runs predate. Full reproduction order:

    # 1. Array probe universes (450K required, EPIC for the supplement).
    #    EPIC needs the annotation package in the epigenomics env first:
    #      BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    cd inputs/supportfiles/_m && mkdir -p logs
    PLATFORM=450K sbatch ../_h/step_1_build_array_universe.sh
    PLATFORM=EPIC sbatch ../_h/step_1_build_array_universe.sh
    # then add the printed rows to _m/annotation_asset_manifest.tsv

    # 2. Module 01 QC refresh -- array coverage and genomic context, on a new
    #    run ID (the accepted catalog runs are sealed and are not modified).
    cd 01_vmr_catalog/_m && mkdir -p logs
    ../_h/submit_qc_refresh.sh

    # 3. Figures. Mints a run ID, builds every figure, seals the run.
    cd 10_integrated_manuscript_outputs/_m && mkdir -p logs
    ../_h/submit_manuscript_figures.sh

Steps 1 and 2 are one-time: once the universes exist and a QC refresh run is
recorded, step 3 alone rebuilds the figures. If the QC refresh is re-run, update
`QC_RUN()` in `01_figure1_catalog.R` and `CTX_RUN()` in
`02_figure2_local_control.R` to the new run IDs.

Both submit wrappers accept `DRY_RUN=1` to print the plan without queueing.

`00_figure_theme.R` holds the shared theme, palette, `save_figure()`, and
`write_source_data()`. It replaces the `BASE_THEME`/`save_plot()` block that was
copy-pasted into ~40 scripts across the three legacy cohort trees; new panels
source it rather than redefining a theme.

## Migrating from

Manuscript figures and consolidated tables currently scattered across `meqtl-validation/12_supplementary_data/` and per-module figure directories.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Products

Manuscript-number registry; main and supplementary tables; main and supplementary
figures; figure source-data tables; analysis-to-claim matrix; exclusions and
denominator table; software and run manifest; manuscript-ready Methods and
Results summaries.

## Rule

Do not manually assemble final figures from files copied across old directories.
Every figure panel must record its source run ID, table, script, and filter.

Always report denominators, exclusions, brain region, donor group, VMR set, and
the exact metric used.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
