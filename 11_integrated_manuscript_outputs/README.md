# 11_integrated_manuscript_outputs — one source of truth

Consumes only accepted immutable upstream runs and produces every manuscript number, table, and figure.

**Status: Figures 1-5, the supplementary figures, Table 1, the cohort QC
panels and every AGENTS.md §7.11 product are implemented. Built as
`fig-all-20260922-c` and accepted 2026-09-22 (see Accepted runs below).**

That run is the first one a reader can reproduce from the run itself. It
carries a `code/` snapshot of `_h/` and `config/` (33 files), records
`git_dirty = false` against the commit that contains the builders, checksums
every file it holds, and seals all of them. Earlier builds did none of that:
`fig-all-20260920` recorded a dirty tree against a commit that did not contain
the builders, leaving an uncommitted working tree as the only record of what
produced the figures. `fig-all-20260922` and `-a`/`-b` were the intermediate
rebuilds that surfaced two seal defects in `00_shared/runid.R::close_run()` --
dotfiles escaping both the checksum manifest and the seal, and an unanchored
exclusion pattern dropping `tables/software-and-run-manifest.tsv` from the
checksums. All four earlier runs are superseded and recorded in
`DEPRECATED_RUNS.tsv`, whose first cleanup tranche deleted them on 2026-09-22.
Their `manifest.tsv` and `output_checksums.tsv` are kept under
`_deleted_run_provenance/`, so the superseded builds remain auditable and this
README's account of them stays checkable after the directories are gone.

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

## Accepted runs

Machine-readable, in the schema `00_shared/gates.R::read_accepted_runs()`
parses. Three things about this table differ from every other module's, and all
three are properties of Module 11 rather than oversights:

- **Nothing consumes it.** Module 11 is terminal, so this row unblocks no
  downstream gate. It records which figures the manuscript cites, and it is what
  flips Supplementary Data 14 from `pending_acceptance` to `ready`.
- **There is no gate script, so `decision` is not a computed token.**
  AGENTS.md §7.11 lists the products this module owes and sets no pass/fail
  criterion, so `ACCEPTED_MANUSCRIPT_OUTPUTS` records a judgement about
  completeness and provenance. It is deliberately not spelled `PASS_*`: every
  other `PASS_*` in this repository was emitted by a gate stage, and borrowing
  the prefix would imply a check that does not exist.
- **One run spans six cells** (two arms × three regions), so `cohort` is `all`
  and `region` is the literal `crossregion`, following Module 08. The six
  `vmr_set_id`s are recorded per arm × region in the run's own
  `tables/exclusions-and-denominators.tsv`, which is sealed and checksummed,
  rather than crushed into one cell here.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| fig-all-20260922-c | all | crossregion | see `tables/exclusions-and-denominators.tsv` (6 cells: arm × region) | 2026-09-22 | Kynon J.M. Benjamin | ACCEPTED_MANUSCRIPT_OUTPUTS | Built at `efaae6220` with `git_dirty = false` and a `code/` snapshot of `_h/` and `config/` (33 files), so the run is reproducible from itself. 23 figures × PDF/SVG/PNG, 74 panel source tables, 11 tables; 190 files, 188 checksummed, 0 writable (the 2 exclusions are the run's own `manifest.tsv` and `output_checksums.tsv`). 40 upstream runs cited: 34 accepted, 6 Module 01 QC-refresh (`vmrcatqc-*-20260826-a`, the documented exception with no acceptance row of their own), 0 unaccepted. Every PDF has matching panel source data; every figure ≤ 9.5 in tall at 7.09/5.51 in column widths; every PDF carries embedded fonts and a ToUnicode map. `manuscript-number-registry.tsv` has 74 rows and **is** Supplementary Data 14; `analysis-to-claim-matrix.tsv` has 64. **No gate script exists for this module**, so this decision certifies completeness, provenance and the claim constraints asserted at build time — not a computed pass. Figure 5 is trait-general with schizophrenia as a marked example; the SCZ locus detail is `figureS_schizophrenia_application` (PI decision 2026-09-22), which honours Module 09's `scz_application_retention = RETAIN_MAIN_TEXT` without giving one trait a main figure. Supersedes `fig-all-20260920` and `fig-all-20260922{,-a,-b}`, all recorded in `DEPRECATED_RUNS.tsv`. |

## Implemented figures

| Output | Content | Upstream runs |
|---|---|---|
| `figure1_vmr_catalog[_all_individuals][_epic]` | **a** cohort **b** VMRs per chromosome **c** VMR width and CpGs per VMR **d** off-array coverage **e** genomic compartment **f** distance to nearest gene | `vmrcat-*-20260816`, `vmrcatqc-*-20260826-a` |
| `figure2_local_genetic_control[_all_individuals]` | **a** estimator concordance vs locus geometry **b** held-out local SNP prediction (end-to-end OOF R²) across rank deciles **c** cross-region rank concordance **d** genic context across the rank | `lgv-AA-*-rescore-20260913` (AA), `lgv-all_individuals-*-20260823`, `lsp-AA-*-20260825` (panel b) |
| `figure3_repeat_repressive_architecture` | **a** the BH family (quiescent, H3K9me3, LINE/L1) **b** complementary contrasts and the H3K27me3 specificity control **c** the five locked analysis sets **d** the continuous gradient | `rra-AA-*-20260906` |
| `figure4_meqtl_burden_coupling` | **a** meQTL-positive CpG fraction across the rank **b** burden model with distal-null λ **c** coupling by modality and predictor **d** coupled-VMR denominators | `cmb-AA-*-20260825`, `tsc-AA-*-20260902` |
| `figure5_gwas_architecture_axis` | **a** every trait's axis estimate by GWAS category **b** schizophrenia against its own null distribution **c** psychiatric vs other traits **d** what a trait's depletion tracks | `scz-AA-*-20260918` + stages 17/18 |
| `figureS_catalog_turnover[...]` | legacy-catalog turnover, **audit only** | as Figure 1 |
| `figureS_local_control_denominators[...]` | denominators and exclusion reasons | as Figure 2 |
| `figureS_local_control_audit_unbounded[...]` | unbounded joint estimate, **audit only** | as Figure 2 |
| `figureS_partitioned_heritability` | S-LDSC across the frozen 8-trait family. **A null**, reported as one | `sldsc-AA-*-20260903` |
| `figureS_aging_axis` | **a** primary age gradient **b** gating sensitivities, with the verdict derived from the run (see below) | `age-AA-*-20260919` |
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

### The aging supplement's verdict is derived, not typed (corrected 2026-09-23)

`figureS_aging_axis` panel b used to carry its conclusion as a string literal:
"The VMR composition-sensitivity arm removes the gradient in every region,
which is why the cross-region token is NOT_SUPPORTED", with a matching claim in
the source table's `row_filter`. A per-region verdict written as prose goes
stale silently the next time the verdict changes, and this one did — the
Module 09b scMD-gate correction is projected to flip DLPFC's reading and move
the stage-05 token off `NOT_SUPPORTED`.

The caption is now built from the run: the token from
`_m/combined/aging-cross-region-decision-{cohort}.tsv`, the supported regions
from `region_supported`, and the failing arms from `fitted`/`survives` in each
run's `gating-sensitivities.tsv`, with the fitted denominator **counted** rather
than asserted as "every region". Where Module 09b supplies a `reason` for a
not-fitted arm, the caption quotes it. The verdict also ships as data on the
panel table (`cross_region_token`, `region_reading`, `region_supported`,
`caption_rendered`), so a reader can check the caption against the numbers it
was built from.

This shares a root cause with the Figure 2 panel b correction above: in both
cases Module 11 stated something its declared upstream run did not supply, and
§7.11's requirement that a panel record its source run, table, script and
filter is satisfied by none of it.

On the sealed numbers the derived caption is already more accurate than the
prose it replaces: `cell_composition_r2` fails in 3 of 3 regions, but
`cell_music` also fails in 1 of 3 and `cell_scmd` in 1 of 1 fitted, which the
single-arm sentence never said.

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

### Which prediction number panel b carries (corrected 2026-09-23)

AGENTS.md §7.3 names one primary v2 prediction endpoint: `r2_pred_oof`, the
**end-to-end** out-of-fold R² from Module 03, in which the locus screen and the
residualization are learned inside the outer training donors too. Module 02 also
emits an `r2_oof` from the nested CV inside its joint-feature elastic net; §4
separates the two standards, and that one is **model-level**.

`fig-all-20260922-c` plotted Module 02's `r2_oof` in panel b under the axis
label "Held-out R²". The two standards are invisible on that axis, and no panel
of any figure in that run named an `lsp-*` run, so Module 03 reached the
manuscript nowhere. Panel b now reads `r2_pred_oof` from the accepted
`lsp-AA-{region}-20260825` runs, joined on `vmr_id` after asserting that
Module 02 and Module 03 agree on `vmr_set_id`; Module 02's `r2_oof` stays in
panel a, relabelled "Model-level OOF R²". Module 03's runs consumed the
pre-rescore `lgv-AA-*-20260823` for their locus screen, which is why the
identity check is on the catalog rather than on the upstream Module 02 run ID --
`r2_pred_oof` is a genotype-to-phenotype quantity and carries no score in it.

The substitution raises the top-decile median in all three regions (caudate
0.870→0.881, DLPFC 0.775→0.794, hippocampus 0.762→0.785). That it is favourable
is not why it was made.

§7.3 also requires that negative `r2_pred_oof` be retained rather than replaced
by `cor2_oof`. Panel b floors nothing and drops nothing: the median and
quartiles are taken on the raw column, most low-decile loci are negative, and
the panel's source table now carries `n_r2_negative` and `n_r2_missing` per
decile so the retention is auditable from the table.

**A new figure run is required for this to reach the manuscript.**
`fig-all-20260922-c` is sealed and still carries the model-level statistic.

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
  Figure 2 panel a's ρ² estimator label, in the file destined for the
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
