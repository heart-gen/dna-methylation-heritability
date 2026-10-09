# 11_integrated_manuscript_outputs — one source of truth

Consumes only accepted immutable upstream runs and produces every manuscript number, table, and figure.

**Status (2026-10-08): accepted as `fig-all-20261008`.** The run cites 51
upstream runs: every module 01-10 at its current accepted run, including 02b
and 02c, Module 01's QC refresh, and Module 06's accepted non-gating positive
control. It replaces `fig-all-20261006-b`. That run predated three
acceptances: the donor-robust tier-2 SE in Module 08, which withdrew the one
tier-2 difference; the span-guard S-LDSC cells; and the external positive
control, which the S-LDSC supplement now shows as panel b.

Four runs were superseded on the way, all recorded in the run-cleanup ledger (kept locally, not distributed):

- `fig-all-20261006-b`: accepted 2026-10-06, superseded 2026-10-08. Its tier-2
  panel claimed a DLPFC-hippocampus difference in meQTL-burden slope that the
  donor-robust SE withdrew, and its S-LDSC caption said there was no positive
  control. Do not cite either.
- `fig-all-20260922-c`: the accepted run before `fig-all-20261006-b`, now stale against its upstreams.
- `fig-all-20261006`: never sealed. It stopped at Figure 2 for the
  `all_individuals` arm, which has no Module 03 run (see **Which prediction
  number panel b carries**).
- `fig-all-20261006-a`: sealed and passed every automated check. Review of the
  rendered figures then found an unlabelled "NA" facet in the Module 08
  sensitivity supplement and two partly hidden in-panel labels (see the
  2026-10-06 entries below).

`fig-all-20260922-c` was the first run a reader can reproduce from the run itself. It
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
the run-cleanup ledger, whose first cleanup tranche deleted them on 2026-09-22.
Their `manifest.tsv` and `output_checksums.tsv` are kept under
the local deleted-run provenance archive, so the superseded builds remain auditable and this
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
  The locked analysis plan lists the products this module owes and sets no pass/fail
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
| fig-all-20261008 | all | crossregion | see `tables/exclusions-and-denominators.tsv` (6 cells: arm × region) | 2026-10-08 | Kynon J. M. Benjamin | ACCEPTED_MANUSCRIPT_OUTPUTS | Built at `47f75d034` with `git_dirty = false` and a `code/` snapshot (39 files). 25 figures × PDF/SVG/PNG, 83 panel source tables, 13 tables; 212 files, 210 checksummed and verified, 0 writable (the 2 exclusions are the run's own `manifest.tsv` and `output_checksums.tsv`). 51 upstream runs cited: 44 accepted gating runs (every module 01-10 including 02b and 02c), 1 accepted non-gating run (`sldsc-AA-external-20261008`), 6 Module 01 QC-refresh (`vmrcatqc-*-20260826-a`), 0 unaccepted. Every PDF has panel source data, embedded fonts and a ToUnicode map; every figure ≤ 9.2 in tall. `manuscript-number-registry.tsv` (83 rows) **is** Supplementary Data 14; `analysis-to-claim-matrix.tsv` has 76, including 4 non-gating runs marked `gating = FALSE`. Against `fig-all-20261006-b`: the region/donor figure reads `rdg-AA-crossregion-20261008` (tier 2 now 0 claimed differences; ceilings from the current EA runs), the S-LDSC supplement reads the span-guard cells and gains panel b, the accepted positive control. Every other panel's plotted data is unchanged. **No gate script exists for this module**, so this decision certifies completeness, provenance and the claim constraints asserted at build time — not a computed pass. Supersedes `fig-all-20261006-b`. |

## Implemented figures

Upstream run IDs are **not listed here**. Every builder resolves them through
`00_shared/gates.R::require_accepted_upstream()` at build time, and each run
records exactly which it used in `tables/software-and-run-manifest.tsv` and in
every panel's `source_run_id`. A run-ID column in this README went stale on every
reacceptance, which is how Figure 2 once shipped on a retired run.

| Output | Content | Owning modules |
|---|---|---|
| `figure1_vmr_catalog[_all_individuals][_epic]` | **a** cohort **b** VMRs per chromosome **c** VMR width and CpGs per VMR **d** off-array coverage **e** genomic compartment **f** distance to nearest gene | 01 (+ the 01 QC refresh) |
| `figure2_local_genetic_control[_all_individuals]` | **a** estimator concordance vs locus geometry **b** held-out local SNP prediction (end-to-end OOF R²) across rank deciles **c** cross-region rank concordance **d** genic context across the rank | 02, 03 |
| `figure3_repeat_repressive_architecture` | **a** the BH family (quiescent, H3K9me3, LINE/L1) **b** complementary contrasts incl. the BrainScope ATAC union, and the H3K27me3 specificity control **c** the locked analysis sets incl. MuSiC (all regions) and scMD (where fitted) **d** the continuous gradient | 04 |
| `figure4_meqtl_burden_coupling` | **a** meQTL-positive CpG fraction across the rank **b** burden model with distal-null λ **c** coupling tests in Module 07's FDR family **d** coupled-VMR denominators | 05, 07 |
| `figure5_gwas_architecture_axis` | **a** every trait's axis estimate by GWAS category **b** schizophrenia against its own null distribution **c** psychiatric vs other traits **d** what a trait's depletion tracks | 09 (+ stages 17/18) |
| `figure_region_donor_generalization` + `_sensitivity` | Module 08 tiers; tier 2 claims only `difference_claimed_primary_claim_family`; the supplement shows every Module 04 analysis set, with the cell arm split into MuSiC (all regions) and scMD (caudate only) | 08 |
| `figureS_catalog_turnover[...]` | legacy-catalog turnover, **audit only** | 01 |
| `figureS_local_control_denominators[...]` | denominators and exclusion reasons | 02 |
| `figureS_local_control_audit_unbounded[...]` | unbounded joint estimate, **audit only** | 02 |
| `figureS_greml_benchmark` | GCTA-GREML recovery of absolute local h2 on **simulated** phenotypes: **a** mean estimate and **b** CI coverage on real AA cis-window genotypes, **c** the v1 AR(1) design (out of regime) | 02b |
| `figureS_cis_greml_sensitivity` | conventional cis-GREML on the **observed** phenotypes against the Module 02 score: **a** mean estimate by score decile, **b** ordering Spearman (score primary; HE and BSLMM descriptive), **c** convergence by decile. Estimator agreement on shared data, not replication; no per-VMR h2 | 02c |
| `figureS_partitioned_heritability` | S-LDSC τ z of the score conditional on VMR membership (two-annotation model). **A null** | 06 |
| `figureS_aging_axis` | **a** primary age gradient, signed-test-only where the magnitude gate withholds a proportion **b** gating sensitivities, verdict derived from the run | 09b |
| `figureS_environmental_axis` | stage B gradients; magnitude gate and the undetermined DLPFC FDR call rendered on the panel | 10 |
| `figureS_schizophrenia_application` | **a** the axis depletion **b** locus evidence tiers **c** prioritized loci | 09 |
| `table1_cohort` (`.tsv`, `.tex`) | donor demographics, both arms x three regions | 01 |
| `figureS_ancestry_pcs`, `figureS_sample_integrity` | genotype PCs over 1000 Genomes; cross-region swap screen | 01 |
| `manuscript-number-registry.tsv` | every citable number -> panel, run, table, column, filter. **This is Supplementary Data 14** | all of the above |
| `analysis-to-claim-matrix.tsv` (`.tex`) | claim -> module -> accepted run -> decision token | module READMEs |
| `exclusions-and-denominators.tsv` | donors, VMRs called, scored, eligible, and why excluded | 01, 02 |
| `supplementary-table-index.tsv` | the tracked `_m/combined/` deliverables | all modules |
| `software-and-run-manifest.tsv` | environment, git commit, upstream run IDs | this run |
| `results-summary.md`, `methods-summary.md` | manuscript-ready summaries; every number cites its registry key | this run |

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
The locked analysis plan forbids exposure results from defining the title, abstract,
primary groups or main causal interpretation.

The caveats are rendered on the panel and read from the accepted
`env-AA-*-20261003` tables rather than typed: the Fieller magnitude gate's
`relative_magnitude_reportable` (no family passes, so no percentage may be
stated), and the PI's ruling that DLPFC nicotine@schizophrenia's FDR call is
undetermined across bootstrap seeds (rendered "?"). The earlier caption's
"these p-values may be too small" was removed: T28 measured the ratio-scale
combined SE as calibrated.
A negative gradient is never evidence that exposure effects concentrate in
weakly controlled VMRs — that is the permanent variance-budget limitation, and it belongs in the Discussion.

### Figure 5 is trait-general, not schizophrenia-specific

The axis depletion is a property of trait-associated loci in general:
schizophrenia sits mid-distribution among the stage 17 traits in every region,
and psychiatric traits are indistinguishable as a category. The ranges are read
from stage 17 under the primary model at build time and printed into the panel
source data and `results-summary.md`; they are not restated here, because a
copy here once mixed the two models' ranges.
Figure 5 therefore shows the **distribution**, with schizophrenia marked in
place as one trait among the rest, which is what the locked analysis plan requires the text to say.

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
result in `_m/`, and this decision does not call for it.

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

The locked analysis plan names one primary v2 prediction endpoint: `r2_pred_oof`, the
**end-to-end** out-of-fold R² from Module 03, in which the locus screen and the
residualization are learned inside the outer training donors too. Module 02 also
emits an `r2_oof` from the nested CV inside its joint-feature elastic net; §4
separates the two standards, and that one is **model-level**.

`fig-all-20260922-c` plotted Module 02's `r2_oof` in panel b under the axis
label "Held-out R²". The two standards are invisible on that axis, and no panel
of any figure in that run named an `lsp-*` run, so Module 03 reached the
manuscript nowhere. Panel b now reads `r2_pred_oof` from the accepted
Module 03 runs, joined on `vmr_id` after asserting that
Module 02 and Module 03 agree on `vmr_set_id`; Module 02's `r2_oof` stays in
panel a, relabelled "Model-level OOF R²". The identity check is on the
catalog rather than on the upstream Module 02 run ID because `r2_pred_oof` is a
genotype-to-phenotype quantity and carries no score in it.

The substitution raises the top-decile median in all three regions (caudate
0.870→0.881, DLPFC 0.775→0.794, hippocampus 0.762→0.785). That it is favourable
is not why it was made.

§7.3 also requires that negative `r2_pred_oof` be retained rather than replaced
by `cor2_oof`. Panel b floors nothing and drops nothing: the median and
quartiles are taken on the raw column, most low-decile loci are negative, and
the panel's source table now carries `n_r2_negative` and `n_r2_missing` per
decile so the retention is auditable from the table.

`fig-all-20260922-c` is sealed and still carries the model-level statistic;
`fig-all-20261006-b` is the first accepted run with panel b on `r2_pred_oof`.

**The `all_individuals` arm has no panel b endpoint.** Module 03 was run for
`AA` and for the donor-group cells `all_individuals.AA` / `.EA`, never for the
pooled `all_individuals` arm. A cell cannot stand in for it: a cell's
`r2_pred_oof` is estimated in one donor group, while the pooled score ranks
every donor, so the join would put two donor sets in one panel. When an arm has
no accepted Module 03 run in any region, panel b is drawn as an explicit "not
run" panel, keeping the letters aligned with the AA figure, and ships no source
data. A partial set still stops at `require_accepted_upstream()`. Accepted as is
by the PI on 2026-10-06; running Module 03 on the pooled arm was not requested.

## The 2026-10-03 rebuild

What changed in the builders, and why each change was needed. In every case the
fix is the same: a reading the panel stated is now read from the upstream table,
or the panel stopped claiming something its upstream no longer licenses.

- **Figure 4.** The header and panel d described PSI coupling as "strong in
  caudate, thin in DLPFC, null in hippocampus" and quoted ABC as "n=312, 19
  coupled". Module 07's PSI identifier-join repair (`tsc-AA-*-20260925-b`)
  withdrew the first reading as an artefact -- splicing coupling holds in all
  three regions -- and both counts were stale. Panel c now renders exactly
  Module 07's FDR family (`in_fdr_family`) and records each excluded modality
  with Module 07's own `fdr_exclusion_reason`.
- **Module 08 figure, panel b.** It claimed and labelled every row with
  `difference_claimed` -- 39 under the repaired accounting, most of them
  `atac_*` rows outside the claim family -- and its label lookup lacked those
  outcomes, so it printed overlapping "NA · held-out R²". Tier 2 now claims
  `difference_claimed_primary_claim_family` (1: the meQTL-burden gradient), and
  the other rows are drawn hollow and unlabelled, as Module 08's README requires.
- **S-LDSC supplement.** It plotted `enrichment`, which for the signed
  continuous score annotation is a ratio over a signed sum, not a share, and
  which Module 06 marks `enrichment_interpretable = FALSE`. It now plots the
  score's τ z conditional on VMR membership, and stops if an uninterpretable
  enrichment would be drawn. Caption: "no detectable enrichment at this footprint".
- **Aging supplement.** Panel a now marks regions whose magnitude the Fieller
  gate withholds as "signed test only" (DLPFC). Panel b no longer prints "the
  age-responsive low-control VMRs are the arm-sensitive ones": that is the
  composition qualifier the 2026-10-01 reacceptance retired.
- **Environmental supplement.** The caption's "no percentage is identifiable"
  is now read from `relative_magnitude_reportable` (0 of 19), and its "these
  p-values may be too small" is gone -- T28 found the ratio-scale SE
  calibrated. DLPFC nicotine@schizophrenia renders "?" rather than "*", because
  the PI ruled its FDR call undetermined across bootstrap seeds. That ruling had
  no column upstream, so it is recorded in this module's
  `config/reporting-constraints.tsv`, snapshotted into each run, and
  `apply_reporting_constraint()` stops the build if it no longer matches
  exactly one row.
- **Figure 3** adds the BrainScope ATAC union to panel b and the scMD arm to
  panel c where fitted; the seven per-cell-type ATAC tracks go to source data
  only (a breakdown, not a cell-type identification).
- **Figure 5** derives its percentile, Wilcoxon and genomic-context ranges from
  stage 17 under the primary model. The old header mixed the two models' ranges
  (16th vs 17th percentile).
- **New:** `figureS_greml_benchmark` (Module 02b),
  `figureS_cis_greml_sensitivity` (Module 02c, added 2026-10-06 at the PI's
  request; its Results line sits in section 2 beside the score, since
  supplementary placement never decides whether a result is reported) and
  `12_methods_results_summaries.R`.

Three more fixes came out of the 2026-10-06 production builds:

- **Figure 2b, `all_individuals` arm.** The first build stopped there: the
  review draft had rendered AA only, so nothing had exercised the second arm
  against the Module 03 join. Panel b is now drawn as "not run" for an arm
  with no Module 03 run (see above).
- **Module 08 sensitivity supplement.** Module 08's rerun added the
  caudate-only `adjust_cell_composition_scmd` set, which the label map lacked,
  so `fig-all-20261006-a` drew a facet titled "NA". The arms are now named
  "MuSiC adjusted" and "scMD adjusted" as in Figure 3, an unlabelled set stops
  the build, and the caption names where scMD was not fitted.
- **Label placement.** The Module 08 panel b claim label now sits above-left
  of the curve with a leader line, since the claimed difference is mid-curve
  after the claim-family repair. The Figure 4b n/λ notes now start right of the
  zero line. Layout only.

The rebuild was accepted as `fig-all-20261006-b` on 2026-10-06.

`tests/test_caption_literals.py` (gitignored) fails on a result number typed
into a caption or row filter, and catches the old Figure 4 string.

## Table 1 and cohort QC (PI decision D3, 2026-08-26)

`04_table1_cohort.R` and `05_qc_sample_integrity.R` replace `sample_summary/`
and the parts of `qc_analysis/` that Module 01 does not cover.

The legacy Table 1 is not merely unmigrated, it reports the **wrong cohort**: it
derived its donor set from `vmr-analysis/all_individuals/{region}/_m/samples.txt`
(invalidated by V1) and honoured the retired sample blacklists. The replacement
reads the donor list from the accepted Module 01 run and hard-stops if the count
disagrees with the locked `design_n` in `config/cohorts.yml`.

### Closed: the readmitted donors (T4)

The 2026-09-22 note here reported 4 of the 8 donors readmitted by retiring the
legacy blacklists as flagged by the cross-region integrity screen, with a Fisher
p of 0.0064. That table was superseded: in `fig-all-20260922-c` the screen flags
9 of 120 donors, and 3 of the 8 readmitted donors (Br1927 is not flagged). The
unit of the test was never prespecified, so no p-value is quoted without its
unit. **The PI cleared all eight donors on 2026-09-25**; no donor is excluded,
`vmr_set_id` is unchanged, and the reasoning is recorded in the
`sample_blacklist` comment block of `config/cohorts.yml`.

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
`12_methods_results_summaries.R` (which cites the registry it built), then
`03_close_figure_run.R`. The old `step_2_table1_qc.sh` ran *after* the seal, so
no sealed run ever contained a `tables/` directory; it has been folded in and
removed.

Both submit wrappers accept `DRY_RUN=1` to print the plan without queueing.

`config/reporting-constraints.tsv` is this module's only configuration: PI
reporting rulings a panel must honour that no upstream table carries as a
column. The submit wrapper snapshots it into `code/module_config/` beside the
`_h/` and repo `config/` snapshots.

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

The legacy v1 trees were retired on 2026-10-06; their tracked files are
recoverable from the annotated tag `v1-legacy-final`.

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

This module follows the repository layout: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
