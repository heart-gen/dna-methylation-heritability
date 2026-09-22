# 11_integrated_manuscript_outputs — one source of truth

Consumes only accepted immutable upstream runs and produces every manuscript number, table, and figure.

**Status: Figures 1-5, the supplementary figures, Table 1, the cohort QC
panels and every AGENTS.md §7.11 product are implemented. Built as
`fig-all-20260920`; no run of this module is accepted yet.**

Figures 1-2 were previously built as `fig-all-20260826-a` on `lgv-AA-*-20260823`,
retired 2026-09-17. That could happen because the builders resolved upstream run
IDs from string templates, so nothing failed when the README of record moved on.
Every builder now resolves them through
`00_shared/gates.R::require_accepted_upstream()`, and
`submit_manuscript_figures.sh` refuses to queue if any cell it needs lacks an
accepted run. The one exception is Module 01's QC-refresh run, which re-runs QC
over an already accepted catalog and so has no acceptance row of its own; it is
named once, in `00_figure_theme.R::QC_REFRESH_RUN()`.

Two upstreams were unblocked on 2026-09-20 rather than worked around:

- **Module 09 stages 17/18.** Their outputs carried `-UNACCEPTED` only because
  they had been built with `--allow-unlocked`.
  `config/gwas_negative_controls.yml` was PI-locked on 2026-09-20, so both
  stages were re-run without the flag (`09/_h/step_9_negative_controls.sh`) and
  their tables are now citable. This is what makes Figure 5 possible.
- **Module 10.** Its three sealed runs were accepted on 2026-09-20 and
  `06_collate_regions.R` re-run without `--allow-unaccepted-runs`, so
  `_m/combined/` carries `citable = TRUE` and `figureS_environmental_axis` is
  wired into the build.

## Implemented figures

| Output | Content | Upstream runs |
|---|---|---|
| `figure1_vmr_catalog[_all_individuals][_epic]` | **a** cohort **b** VMRs per chromosome **c** VMR width and CpGs per VMR **d** off-array coverage **e** genomic compartment **f** distance to nearest gene | `vmrcat-*-20260816`, `vmrcatqc-*-20260826-a` |
| `figure2_local_genetic_control[_all_individuals]` | **a** estimator concordance vs locus geometry **b** held-out R² across rank deciles **c** cross-region rank concordance **d** genic context across the rank | `lgv-AA-*-rescore-20260913` (AA), `lgv-all_individuals-*-20260823` |
| `figure3_repeat_repressive_architecture` | **a** the BH family (quiescent, H3K9me3, LINE/L1) **b** complementary contrasts and the H3K27me3 specificity control **c** the five locked analysis sets **d** the continuous gradient | `rra-AA-*-20260906` |
| `figure4_meqtl_burden_coupling` | **a** meQTL-positive CpG fraction across the rank **b** burden model with distal-null λ **c** coupling by modality and predictor **d** coupled-VMR denominators | `cmb-AA-*-20260825`, `tsc-AA-*-20260902` |
| `figure5_gwas_architecture_axis` | **a** every trait's axis estimate by GWAS category **b** schizophrenia against its own null distribution **c** psychiatric vs other traits **d** what a trait's depletion tracks | `scz-AA-*-20260918` + stages 17/18 |
| `figureS_catalog_turnover[...]` | legacy-catalog turnover, **audit only** | as Figure 1 |
| `figureS_local_control_denominators[...]` | denominators and exclusion reasons | as Figure 2 |
| `figureS_local_control_audit_unbounded[...]` | unbounded joint estimate, **audit only** | as Figure 2 |
| `figureS_partitioned_heritability` | S-LDSC across the frozen 8-trait family. **A null**, reported as one | `sldsc-AA-*-20260903` |
| `figureS_aging_axis` | **a** primary age gradient **b** gating sensitivities, incl. the composition arm that removes it | `age-AA-*-20260919` |
| `figureS_environmental_axis` | stage B proportional gradients, with both acceptance caveats on the panel | `env-AA-*-20260920-a` |
| `figureS_schizophrenia_application` | **a** the axis depletion **b** locus evidence tiers **c** prioritized loci | `scz-AA-*-20260918` |
| `figure_region_donor_generalization` + `_sensitivity` | Module 08 tiers | `rdg-AA-crossregion-20260918` |
| `table1_cohort` (`.tsv`, `.tex`) | donor demographics, both arms x three regions | `vmrcat-*-20260816` |
| `figureS_ancestry_pcs`, `figureS_sample_integrity` | genotype PCs over 1000 Genomes; cross-region swap screen | `vmrcat-*-20260816` |
| `manuscript-number-registry.tsv` | every citable number -> panel, run, table, column, filter. **This is Supplementary Data 14** | all of the above |
| `analysis-to-claim-matrix.tsv` (`.tex`) | claim -> module -> accepted run -> decision token | module READMEs |
| `exclusions-and-denominators.tsv` | donors, VMRs called, scored, eligible, and why excluded | Modules 01, 02 |
| `supplementary-table-index.tsv` | the tracked `_m/combined/` deliverables | all modules |
| `software-and-run-manifest.tsv` | environment, git commit, upstream run IDs | this run |

AA is the primary arm; `all_individuals` renders from the same builders as the
sensitivity supplement.

### The environmental supplement

Built, and wired into the run. `10_environmental_exploratory` was accepted on
2026-09-20, so `figureS_environmental_axis` renders from the citable
`_m/combined/` tables. It stays a supplement whatever it shows: the module's
decision row carries `main_text_retention = NEVER_SUPPLEMENT_ONLY`, and
AGENTS.md §2.3 forbids exposure results from defining the title, abstract,
primary groups or main causal interpretation.

Both acceptance caveats are rendered on the panel rather than left to the
legend: no percentage is formally identifiable (`mean_omega_z` never reaches
1.96), and the donor bootstrap inflates the ratio's denominator 2.0×-38.5×.
A negative gradient is never evidence that exposure effects concentrate in
weakly controlled VMRs — that is the permanent variance-budget limitation
(AGENTS.md §7.10), and it belongs in the Discussion.

### Figure 5 is trait-general, not schizophrenia-specific

The axis depletion is a property of trait-associated loci in general:
schizophrenia sits at the 16th-33rd percentile of 63 GWAS traits by region and
psychiatric traits are indistinguishable as a category (Wilcoxon p 0.11-0.96).
Figure 5 therefore shows the **distribution**, with schizophrenia marked in
place as one trait among the rest, which is what AGENTS.md §7.8 rule 1 and
§11 require the text to say.

**PI decision, 2026-09-22: the GWAS collection is the main-text result and the
schizophrenia detail is supplemental.** This is the split the build already
implements, now recorded rather than inferred.

Module 09's decision 2 is `scz_application_retention = RETAIN_MAIN_TEXT`, and
that is honoured: schizophrenia appears in Figure 5 panels a, b and c, as one
trait among the rest. What sits in `figureS_schizophrenia_application` is the
locus-level detail — evidence tiers, prioritized loci, the per-region axis
contrast — which is specific to the one trait and would otherwise crowd out
the general result.

The two decisions are compatible and neither overrides the other:
`RETAIN_MAIN_TEXT` requires schizophrenia to *appear* in the main text, not to
*own* a figure. Nothing here withdraws the schizophrenia application, and the
Module 09 decision row is untouched — an agent must not edit a scientific
result in `_m/` (AGENTS.md §5.2), and this decision does not call for it.

Two constraints the builder asserts at runtime, because both are easy to get
wrong from memory:

- **No trait reaches q < 0.05 for LINE/L1 *enrichment*** in any region or arm
  (0 of 888 per-trait tests). The *axis-link* Spearman is a different table and
  does have nominally significant LINE/L1 rows in the all-VMR arm — in
  inconsistent directions across regions, collapsing under high mappability.
  Panel d shows both arms so that collapse is visible, and nothing here
  licenses a repeat statement.
- The GWAS collection is European or European-dominated while the cohort is
  admixed African American. The limitation attaches to the **locus definition**,
  not to the axis, which is a within-cohort rank.

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
    cd 11_integrated_manuscript_outputs/_m && mkdir -p logs
    ../_h/submit_manuscript_figures.sh

Steps 1 and 2 are one-time: once the universes exist and a QC refresh run is
recorded, step 3 alone rebuilds the figures. If the QC refresh is re-run, update
`QC_REFRESH_RUN()` in `00_figure_theme.R` — one place, not two — to the new
run ID. Every other upstream run ID is resolved through the acceptance gate and
needs no edit here when a module supersedes a run.

`step_1_figures.sh` builds in one order and seals last: both-arm figures, then
the AA-only Figures 3-5 and supplements, then Table 1 and the QC panels, then
`10_manuscript_tables.R` (which reads what the figures actually rendered), then
`03_close_figure_run.R`. The old `step_2_table1_qc.sh` ran *after* the seal, so
no sealed run ever contained a `tables/` directory; it has been folded in and
removed.

Both submit wrappers accept `DRY_RUN=1` to print the plan without queueing.

`00_figure_theme.R` holds the shared theme, palette, `save_figure()`,
`write_source_data()`, `sig_stars()`, `scale_fill_log2or()` and `fig_tags()`.
It replaces the `BASE_THEME`/`save_plot()` block that was copy-pasted into ~40
scripts across the three legacy cohort trees; new panels source it rather than
redefining a theme. It lives in `_h/` rather than `00_shared/` deliberately:
`_m/runs/*/code/` snapshots `_h/` and `config/` but not `00_shared/`, which is
how a shared-code defect reached three sealed Module 10 runs with no provenance
trail on 2026-09-20.

Four things it settles that the v1 tree did not:

- **Panel tags are lowercase.** `content/91.figure-legends.md` cites panels as
  **a.**, **b.**; v1 rendered them uppercase and relabelled by hand during
  manual assembly. `fig_tags()` emits them lowercase so that step is gone.
- **The PDF device is `cairo_pdf`.** Base `pdf()` writes a single-byte encoding
  and silently drops anything it cannot map — it was dropping the ρ in
  Figure 2's "Out-of-fold ρ²" axis label, in the file destined for the
  journal, while the PNG review copy rendered correctly.
- **Figures are capped at 9.5 in tall** and `save_figure()` stops above it. The
  first v2 drafts of Figures 1 and 2 were 11.4 and 10.2 in, which no journal
  page accommodates.
- **Every device gets `bg = "white"`**, and an SVG is written alongside the PDF
  and PNG for the manubot manuscript build.

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
