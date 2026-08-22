AGENTS.md — DNAm local-genetic-control revision

## 1. Scope and authority

This file governs the AJHG revision of the
`heart-gen/dna-methylation-heritability` repository. It applies from the
repository root unless a more specific `AGENTS.md` exists in a child directory.

The revision is **biology-forward**. Statistical calibration and local SNP
prediction are enabling analyses, not the manuscript's destination. All work
must serve the central biological question:

> How does local genetic variation organize interindividual DNA methylation
> variability across repeat-rich and repressive compartments of the adult
> human brain, and how does this architecture connect schizophrenia-risk
> variants to methylation and transcription?

Before changing code, read these files completely when available:

- `writing-notes/REVISION_GUIDE.md`
- `writing-notes/PIPELINE_AUDIT.md`
- `calibrated-simulation-analysis/README.md`
- `meqtl-validation/ANALYSIS_SUMMARY_meqtl_validation.md`
- `meqtl-validation/07_repeat_mappability_sensitivity/README.md`
- `meqtl-validation/07_repeat_mappability_sensitivity/CELLTYPE_LINE_L1_PLAN.md`
- `meqtl-validation/08_schizophrenia_risk_application/README.md`
- `meqtl-validation/08_schizophrenia_risk_application/_m/decision/PHASE7_DECISION.md`
- the README for the analysis being changed

Some revision notes may exist only on Quest. If a required note is unavailable,
record that fact and do not invent missing decisions.

## 2. Manuscript identity

### 2.1 Working biological thesis

The revision should test and, if supported, establish this sequence:

1. Brain VMRs show a continuous spectrum of local SNP-explained methylation
   variance.
2. Greater local genetic control is concentrated in H3K9me3-marked,
   quiescent, and repeat-rich genomic compartments.
3. LINE/L1 enrichment is robust in caudate and hippocampus but may be
   regionally heterogeneous rather than uniformly brain-wide.
4. VMRs with stronger local genetic control contain greater CpG cis-meQTL
   burden and are more often coupled to expression or splicing.
5. Schizophrenia-risk variants regulate methylation within a subset of these
   genetically anchored VMRs, including loci with transcriptional and external
   eQTL support.

### 2.2 Contribution to the field

The paper must not be framed as the first elastic-net DNAm predictor or the
first postmortem-brain meQTL map. Those precedents already exist.

The distinctive contribution is the combination of:

- WGBS coverage of methylation outside array-accessible CpGs;
- region-level variable methylation phenotypes rather than isolated CpGs;
- three disease-relevant brain regions;
- substantial representation of Black American donors;
- calibrated aggregate local SNP-explained variance;
- repeat- and heterochromatin-centered biological interpretation;
- CpG meQTL and transcription/splicing coupling;
- a prespecified schizophrenia-risk application.

Prediction is used to define which VMR phenotypes can be imputed into new
cohorts. It is secondary to the biological architecture.

### 2.3 Prohibited framing

Do not restore the previous genetically anchored versus exposure-associated
binary classification. Low local SNP variance or poor SNP prediction is not
evidence that a VMR is environmentally determined.

Do not claim that:

- local genetic control is uniformly high across VMRs;
- LINE/L1 enrichment is identical across all brain regions;
- repeat overlap demonstrates active retrotransposition;
- bulk-tissue results establish a cell type of origin;
- methylation mediates schizophrenia risk;
- an elastic-net SNP is a causal variant;
- donor-group differences establish ancestry-specific biology;
- WGBS distinguishes 5mC from 5hmC.

Exposure results may be retained only as descriptive or sensitivity analyses in
the supplement. They must not define the title, abstract, primary groups, or
main causal interpretation.

## 3. Legacy quantities and eventual removal

The following quantities and derived classifications are invalid for active v2
scientific use:

- legacy `r_squared_cv`;
- `h2_unscaled`;
- `r_squared_cv > 0.75`;
- high- versus low-prediction groups derived from that threshold;
- stacked, Venn, Sankey, enrichment, exposure, or disease results whose group
  definition depends on the withdrawn metric;
- `h2_en_calibrated` and any absolute locus-level PVE percentage or threshold
  derived from it (retired 2026-08-21; see 7.2). The raw
  `pve_cis_joint_calibrated` estimate is retained inside Module 02 runs for
  audit only and must not be reported as a percentage of variance explained.

Use a two-stage retirement process:

1. **During validation:** keep legacy code and outputs temporarily, clearly
   labeled `legacy_invalid_for_prediction_accuracy`, so old-versus-new
   comparisons and migration audits are possible.
2. **Before submission:** remove superseded legacy outputs and obsolete active
   scripts from the submission branch after every downstream consumer has a
   validated v2 replacement. Historical recovery must remain possible from a
   tagged Git commit, archived release, or existing Zenodo record.

Do not retain invalid generated outputs merely because they were previously
tracked. Do not delete a legacy directory until its row in
`MIGRATION_MANIFEST.tsv` is marked `validated_replacement` and
all consumers point to accepted v2 run IDs.

Active revision code must fail if `r_squared_cv`, `h2_unscaled`, or
`h2_en_calibrated` appears as an input to a v2 model, table, or figure. If a historical training-set statistic is
temporarily retained, name it `r_squared_insample` and never call it held-out,
cross-validated, or predictive accuracy.

All prior biological findings conditioned on the legacy metric are hypotheses
to retest, not numbers that can be carried into the revision.

## 4. Literature and novelty boundary

Agents must interpret results against these precedents:

- PrediXcan/PredictDB and FUSION establish cis-SNP molecular prediction,
  heritability screening, and held-out model evaluation.
- Released PredictDB code includes outer prediction folds with inner
  `cv.glmnet` tuning.
- FUSION reports fivefold out-of-sample model accuracy after full-data
  cis-heritability screening.
- The 2026 multi-ancestry monocyte WGBS atlas uses cis-heritability screening,
  elastic-net DNAm prediction, and MWAS.
- BrainSeq WGBS established widespread CpG and CpH meQTLs in DLPFC and
  hippocampus and connected brain meQTLs to schizophrenia risk.
- Postmortem neuronal and non-neuronal studies establish that L1 methylation is
  cell-contextual and potentially relevant to psychiatric disease.

Distinguish two evaluation standards:

- **Model-level out-of-fold:** held-out donors are not used to fit the final
  regression, but full-data screening or preprocessing may occur upstream.
- **End-to-end out-of-fold:** every phenotype-informed operation, including
  locus screening and residualization, is learned in outer-training donors.

Use end-to-end out-of-fold evaluation as the primary v2 prediction standard.
A literature-comparable full-data-screen sensitivity analysis is allowed if
clearly labeled. Do not claim that v2 is the first analysis to use OOF
prediction.

## 5. Active v2 repository organization

### 5.1 Revision modules and their order

New v2 work lives in numbered modules at the **repository root**:

```text
<repo root>/
├── AGENTS.md
├── MIGRATION_MANIFEST.tsv
├── config/                       # shared; extended with the v2 keys
│   ├── paths.yml
│   ├── cohorts.yml
│   ├── covariates.yml
│   ├── thresholds.yml
│   ├── prediction.yml
│   ├── repeat_annotations.yml
│   └── schizophrenia.yml
├── 00_shared/                    # library code; no _h/, since it emits no _m/
│   └── tests/
├── 01_vmr_catalog/
├── 02_local_genetic_variance/
├── 03_local_snp_prediction/
├── 04_repeat_repressive_architecture/
├── 05_cpg_meqtl_burden/
├── 06_transcription_splicing_coupling/
├── 07_region_donor_generalization/
├── 08_schizophrenia_risk_application/
└── 09_integrated_manuscript_outputs/
```

The numeric prefixes encode scientific dependency, in the same spirit as the
ordering used in `meqtl-validation`. Do not assign numbers based on perceived
importance or manuscript figure order.

Configuration is shared at `config/`, not duplicated per module. A module gets
its own `config/` only for files that are meaningless outside it (for example
`02_local_genetic_variance/config/acceptance-criteria.tsv`).

### 5.2 Required structure for every numbered analysis

Every analysis directory must use this contract:

```text
NN_analysis_name/
├── README.md
├── _h/                               # human-written executable code
│   ├── 00_prepare.*
│   ├── 01_analyze.*
│   ├── 02_summarize.*
│   ├── 03_plot.*
│   └── step_*.sh                     # SLURM launchers
├── _m/                               # machine-generated output only
│   ├── README.md
│   └── runs/{RUN_ID}/
├── config/                            # only when analysis-specific
└── tests/                             # gitignored smoke checks
```

Rules:

- `_h/` contains R, Python, shell, workflow, and submission scripts written by
  people or agents.
- The `_h/`/`_m/` split exists to separate code from generated results. A
  directory that produces **no** `_m/` therefore has no `_h/`: `00_shared/`
  holds its library files directly at its top level.
- `tests/` is gitignored everywhere. Tests are quick smoke checks run by hand
  before submitting a full SLURM array, not a tracked CI suite. They must still
  exist and pass locally before any production run.
- Run IDs are `{module}-{cohort}-{region}-{YYYYMMDD}`, with a letter suffix if
  more than one run lands on the same day. No `v{N}` component.
- `_m/` contains only generated tables, figures, logs, manifests, checkpoints,
  and run metadata.
- Do not manually edit a scientific result in `_m/`.
- A generated file must be reproducible from `_h/`, locked configuration, and
  declared inputs.
- Each `_m/README.md` documents what is tracked in Git, what remains on Quest,
  and how a run can be regenerated.
- Use immutable `_m/runs/{RUN_ID}/` directories. Never update a completed run
  in place.
- Do not store scripts inside `_m/` or generated results inside `_h/`.

### 5.3 Migration from current directories

Do not perform a single destructive repository-wide move. Build and validate
the v2 modules incrementally while retaining legacy material only for the
temporary comparison period described in Section 3.

`MIGRATION_MANIFEST.tsv` must contain:

- legacy path;
- v2 replacement path;
- scientific role;
- legacy commit SHA;
- replacement status;
- accepted v2 run ID;
- validation status;
- downstream consumers;
- retirement status.

Recommended mapping:

| Existing source | Active v2 destination |
|---|---|
| `vmr-analysis/` | `01_vmr_catalog/` |
| `calibrated-simulation-analysis/` | `02_local_genetic_variance/` |
| `local-snp-prediction/` | `03_local_snp_prediction/` |
| repeat and cell-composition modules in `meqtl-validation/` | `04_repeat_repressive_architecture/` |
| CpG mapping and VMR burden modules in `meqtl-validation/` | `05_cpg_meqtl_burden/` |
| expression/splicing modules in `meqtl-validation/` | `06_transcription_splicing_coupling/` |
| donor-group, cross-region, and downsampling modules | `07_region_donor_generalization/` |
| `meqtl-validation/08_schizophrenia_risk_application/` | `08_schizophrenia_risk_application/` |
| manuscript figures and consolidated tables | `09_integrated_manuscript_outputs/` |
| `environmental-analysis/` | remove before submission unless an approved supplemental sensitivity is migrated |

Do not duplicate region- or cohort-specific copies of the same code. Refactor
shared logic into parameterized functions in `00_shared/` or the owning
analysis `_h/`.

## 6. Dependency order and production gates

Analyses must run in this order:

1. `01_vmr_catalog`
2. `02_local_genetic_variance`
3. `03_local_snp_prediction`
4. `04_repeat_repressive_architecture`
5. `05_cpg_meqtl_burden`
6. `06_transcription_splicing_coupling`
7. `07_region_donor_generalization`
8. `08_schizophrenia_risk_application`
9. `09_integrated_manuscript_outputs`

No downstream production run may consume an upstream result until the upstream
README records a passing acceptance gate and immutable run ID.

VMR turnover never authorizes reuse of downstream numbers. Correcting VMR
definitions requires recomputing every analysis based on VMR membership,
boundaries, summaries, or classification.

## 7. Analysis contracts

### 7.1 `01_vmr_catalog`: corrected biological units

The primary catalog should use the PI-locked cohort and chromosome policy.
Default recommendation pending PI lock:

- corrected Black American VMR catalog as primary;
- all-individual VMR catalog as sensitivity;
- autosomes as the primary analysis;
- sex chromosomes reported separately or excluded with an explicit manifest.

Mandatory repairs:

- align methylation, phenotype, genotype-PC, and covariate rows by donor ID;
- fail on missing, duplicate, or reordered identifiers;
- remove hard-coded region filters and stale paths;
- use correct region-specific blacklists;
- ensure functions operate on their declared arguments;
- guard top-CpG indexing by observed CpG count;
- sort chromosomes and chunks numerically;
- verify identical FID/IID order before combining outputs;
- avoid unpinned `package:::` internals;
- use one parameterized implementation across regions and cohorts.

Required outputs:

- corrected VMR BED and phenotype matrix per region;
- CpG membership table;
- donor and covariate manifests;
- old-versus-new turnover table;
- array-coverage comparison;
- technical QC and exclusion table;
- immutable `vmr_set_id` propagated downstream.

### 7.2 `02_local_genetic_variance`: primary quantitative endpoint

**Superseded 2026-08-21.** The endpoint was `h2_en_calibrated`, a
simulation-calibrated elastic-net estimate of local SNP-explained methylation
variance. Independent zero-overlap validation (12,960/12,960 simulations)
returned `PASS_RELATIVE_GENETIC_CONTROL / FAIL_ABSOLUTE_LOCUS_PVE`: strong
ordering (Spearman 0.796) but 6/14 absolute-PVE gates failed (mean bias 0.0869,
RMSE 0.2127, null mean 0.0767). Estimator development has ended and
`h2_en_calibrated` is retired as a reported quantity.

The primary endpoint is now the within-cell relative score, described as:

> relative rank of local SNP contribution to methylation variance, within
> cohort by region

emitted by Module 02 as `local_snp_contribution_score` (empirical midrank
percentile), `local_snp_contribution_score_z` (primary model predictor), and
`local_snp_contribution_quartile` (secondary contrast).

Requirements:

- consume the frozen joint model by its pinned SHA-256 and apply it once;
- rank within each cohort-by-region cell, never pooled across cells;
- numerically interpret only within-domain loci, and report the domain coverage;
- retain the raw `pve_cis_joint_calibrated` estimate for audit only;
- carry `absolute_pve_interpretation_allowed = FALSE` on every emitted row;
- stop if the Stage 05 observed-score gate fails;
- rerun annotation sensitivity on corrected VMRs.

Prohibited: exact locus-level PVE percentages; absolute PVE thresholds
including 0.10; heritable/nonheritable classes; raw score-level comparison
across regions. See `02_local_genetic_variance/README.md` and
`config/local_genetic_control.yml` for the lock.

Changing the estimator, alpha grid, folds, repeats, lambda rule, screen, or raw
metric requires recalibration.

Treat local genetic control as continuous. Any quantile or extreme-group
presentation is secondary to the continuous model and must use a prespecified
threshold.

### 7.3 `03_local_snp_prediction`: secondary translational endpoint

Prediction answers whether a VMR can be imputed in a new cohort. It is not the
primary biological endpoint.

The primary estimate must be end-to-end out-of-fold. Within each outer split,
held-out donors must not influence:

- phenotype residualization, centering, or scaling;
- genotype imputation, scaling, MAF, missingness, or zero-variance filters;
- phenotype-informed SNP selection;
- locus-level cis-genetic signal screening;
- alpha, lambda, early stopping, or model-family selection.

All tuning occurs inside outer-training donors. A VMR that fails its
fold-internal screen receives the prespecified null prediction for held-out
donors; do not drop failed folds.

Required metrics:

- `r2_pred_oof = 1 - SSE/SST` from pooled OOF predictions;
- `cor2_oof` from pooled OOF predictions;
- RMSE and MAE;
- calibration intercept and slope;
- screening pass frequency;
- number of predictions per donor;
- fold and repeat diagnostics.

Negative `r2_pred_oof` values must be retained. Do not substitute `cor2_oof`
when `r2_pred_oof` is unfavorable.

Do not use the current ordinary-OLS Haseman-Elston p-value as the production
screen without donor-robust or permutation calibration. Lock the screen and
multiple-testing rule before viewing final v2 results.

A full-data-screen followed by nested OOF model evaluation may be included as a
literature-comparable sensitivity analysis. Label it explicitly.

### 7.4 `04_repeat_repressive_architecture`: primary biological analysis

This module tests whether increasing calibrated local genetic control is
associated with repeat-rich and repressive genomic compartments.

Primary outcomes:

- H3K9me3 overlap;
- quiescent-chromatin overlap;
- LINE/L1 overlap or overlap fraction.

Primary predictor:

- standardized continuous `local_snp_contribution_score_z` among eligible
  loci.

Secondary predictor:

- honest `r2_pred_oof`, used to determine whether the same compartments are
  enriched among imputable VMRs.

Required adjustment or matching variables include, where available:

- VMR length;
- CpG count and density;
- GC content;
- mean methylation and methylation variance;
- WGBS coverage;
- local tested-SNP opportunity and SNP proximity;
- mappability;
- segmental duplication and problematic-region overlap;
- broad genomic annotation;
- cell-composition-associated methylation properties.

Required sensitivities:

- high-mappability restriction;
- exclusion of SNP-proximal CpGs;
- exclusion of segmental duplications;
- RNA MuSiC cell-composition adjustment;
- caudate DNAm scMD adjustment when the integration gate passes;
- matched high- versus low-control VMR comparison as secondary evidence;
- public cell-type-resolved brain methylome overlap when feasible.

Interpretation gate:

- H3K9me3 and quiescent enrichment may be described as shared across regions
  only if direction and inference survive locked sensitivities in all three.
- LINE/L1 may be described as a multi-region result if it survives in at least
  two regions. If it survives only in caudate, present it as caudate-specific.
- Never hide a DLPFC reversal or null after the high-mappability restriction.

High-value extension, if technically reliable:

- separate LINE/L1 by RepeatMasker class, family, and subfamily;
- compare younger L1HS/L1PA with older L1M families;
- distinguish full-length from truncated elements;
- distinguish 5-prime promoter-containing from internal or 3-prime fragments;
- test joint LINE/L1 and H3K9me3 annotation;
- use a validated reference list for retrotransposition-competent L1s.

Do not infer activity, expression, or retrotransposition from overlap alone.

### 7.5 `05_cpg_meqtl_burden`: orthogonal genetic evidence

This module asks whether increasing calibrated VMR local genetic control is
associated with a greater fraction of constituent CpGs having conventional cis
meQTL support.

Requirements:

- use corrected CpG-to-VMR membership;
- preserve the locked donor and covariate models;
- report tested CpGs separately from prepared but untested CpGs;
- model continuous calibrated local variance as the primary predictor;
- use overdispersion-appropriate and donor-robust inference;
- adjust for tested CpG count, coverage, methylation properties, VMR geometry,
  SNP opportunity, mappability, repeats, and broad genomic annotation;
- retain matched extreme-group analysis only as secondary evidence;
- audit all concordance denominators;
- resolve genomic inflation before figure freeze.

Internal meQTL mapping is convergent evidence, not independent replication.
Positive-only public resources cannot provide an external gradient because
absence from the positive list is not a tested negative.

### 7.6 `06_transcription_splicing_coupling`: regulatory consequences

Test whether meQTL-supported or locally controlled VMRs are more likely to have
existing significant associations with:

- gene or transcript abundance;
- transcript usage or splicing.

Do not initiate an unbounded transcriptome-wide fishing analysis. Reuse the
prespecified expression and splicing analyses and clearly document their tested
universe.

Adjust for the number of tested features, VMR length, VMR-to-feature distance,
methylation variance, local SNP number, and applicable technical factors.

Allowed interpretation:

> Genetically regulated VMRs are more frequently transcriptionally coupled.

Forbidden interpretation:

> Methylation mediates the genetic effect on expression or splicing.

### 7.7 `07_region_donor_generalization`: boundaries of the biology

This module should establish what is shared and what is context-dependent
across:

- caudate, DLPFC, and hippocampus;
- Black American and white American donor analyses;
- matched donor subsets and sample-size-matched analyses.

Prioritize biological generalization of local variance, repeat enrichment,
meQTL burden, and effect direction. Cross-population predictor portability is
optional and secondary.

Do not attribute differences to ancestry-specific biology without eliminating
sample size, MAF, LD, SNP availability, assay, covariate, and brain-region
explanations. Use donor-group or population language approved by the PI.

### 7.8 `08_schizophrenia_risk_application`: required disease application

Phase 7 is part of the intended main-text biological story, conditional on its
results surviving corrected VMRs and the new local-genetic-control axis.

Preserve these design principles:

- define PGC schizophrenia loci independently of methylation results;
- keep risk-variant-CpG tests in their own FDR family;
- link risk variants to corrected VMRs through explicit CpG meQTL evidence;
- test association with calibrated local variance, not legacy predictability;
- evaluate transcription and splicing coupling;
- retain cross-region, shared-donor genotype-by-region, and caudate
  downsampling analyses;
- use external GTEx eQTL evidence as support, not proof of mediation;
- prioritize at most five illustrative loci using prespecified criteria.

Existing hero loci `rs8048039` and `rs13331198` remain candidates, not fixed
answers. Reprioritize after the corrected run.

Add an explicit integration analysis asking whether schizophrenia-linked VMRs
are enriched for:

- LINE/L1;
- H3K9me3;
- quiescent chromatin;
- high-mappability repeat intervals;
- expression or splicing coupling.

If this integration is positive, connect the disease application directly to
the repeat/repressive architecture. If it is null, present Phase 7 as a separate
proof of disease relevance.

Retain the main-text application only if the corrected run meets prespecified
Phase 7 criteria, including:

- at least one schizophrenia-risk locus with CpG meQTL support;
- enrichment along the new local-genetic-control axis in the primary region;
- at least one locus with transcriptional coupling;
- evidence that the caudate result is not solely a sample-size artifact;
- external shared genetic support for at least one prioritized locus.

Diagnosis associations are optional supportive evidence. A null best-VMR
diagnosis result does not invalidate genetic regulation. Secondary diagnosis
hits must be labeled exploratory.

Do not claim mediation, causality, or colocalization unless the corresponding
analysis has been run with adequate ancestry-matched LD and passes its own gate.

### 7.9 `09_integrated_manuscript_outputs`: one source of truth

This module consumes only accepted immutable upstream runs and creates:

- manuscript-number registry;
- main and supplementary tables;
- main and supplementary figures;
- figure source-data tables;
- analysis-to-claim matrix;
- exclusions and denominator table;
- software and run manifest;
- manuscript-ready Methods and Results summaries.

Do not manually assemble final figures from files copied across old directories.
Every figure panel must record its source run ID, table, script, and filter.

## 8. Prior results to retest, not assume

These results motivate the v2 biological hypotheses and must be recomputed
after the VMR and metric repairs:

| Existing result | v2 test |
|---|---|
| H3K9me3 enrichment in all three regions | continuous `local_snp_contribution_score_z` plus locked sensitivities |
| Quiescent-chromatin enrichment in all three regions | relative local-genetic-control score plus cell and mapping sensitivities |
| LINE/L1 enrichment in caudate and hippocampus; fragile DLPFC result | repeat-family model with high-mappability and cell-adjusted gates |
| Strong CpG meQTL-burden gradient | relative local-genetic-control score primary; OOF prediction secondary |
| Expression coupling in all regions and splicing coupling in caudate/DLPFC | corrected VMR membership and tested-universe model |
| Aggregate donor-group consistency | matched MAF/SNP opportunity and denominator audit |
| Phase 7: 31 caudate loci, 361 pairs, 38 VMRs, eight TX-coupled VMRs | complete rerun and reprioritization on corrected VMRs |
| Hero loci rs8048039 and rs13331198 | retain only if they pass corrected prioritization |

Do not preserve a numerical result merely because it is favorable or already
appears in a manuscript draft.

## 9. Reproducibility and run identity

Every production output must carry or link to:

- run ID and analysis name;
- Git commit;
- configuration checksum;
- input and annotation checksums;
- upstream run IDs;
- `vmr_set_id`;
- ordered donor-list checksum;
- region, donor group, and sample count;
- residualization formula;
- cis-window and SNP QC;
- fold assignment and random seeds when applicable;
- software environment;
- SLURM job identifiers;
- task reconciliation counts;
- output checksums.

Use deterministic seeds derived from run ID, region, VMR task, repeat, and fold.
Write outputs atomically. A combine step must reconcile expected, completed,
excluded, QC-failed, and computationally failed tasks. Production runs have
zero tolerance for unexplained computational failures.

Quest paths belong in configuration or environment variables, never inside
analysis functions.

## 10. Required tests

### 10.1 VMR and identity tests

- Shuffling input rows does not change results after ID-based alignment.
- Missing or duplicate donors fail loudly.
- Region-mismatched blacklists or paths fail.
- Observed sample counts match locked design counts.
- Chromosome ordering is numeric and deterministic.
- A small end-to-end run reproduces expected VMRs and checksums.

### 10.2 Variance and prediction tests

- Outer test donors never appear in fitting, tuning, screening, or
  phenotype-informed preprocessing.
- Imputation and QC are training-only.
- Fold assignment is deterministic.
- Null simulations yield approximately null held-out performance.
- Positive controls recover signal.
- Negative `r2_pred_oof` is retained.
- `r2_pred_oof` and `cor2_oof` are calculated distinctly and correctly.
- Active outputs do not contain legacy metrics.

### 10.3 Enrichment and disease tests

- Genomic coordinates and genome builds are consistent.
- Repeat overlaps are invariant to interval ordering.
- Mappability, segmental-duplication, and SNP-proximity exclusions use locked
  definitions.
- Enrichment nulls preserve VMR length, CpG density, and applicable genomic
  structure.
- Cell-composition sensitivity uses the declared sample set and PCs.
- CpG-to-VMR membership uses the accepted `vmr_set_id`.
- Schizophrenia loci are defined independently of methylation outcomes.
- FDR families are not combined after inspection.
- Cross-region and donor-group concordance denominators are explicit.

Run unit tests and a small smoke workflow before full SLURM arrays. A completed
SLURM job is not proof of scientific validity.

## 11. Manuscript and figure organization

The main Results should follow the biology rather than repository chronology:

1. A corrected WGBS VMR catalog across three brain regions.
2. Continuous local SNP control of VMR methylation variability.
3. Concentration of local genetic control in repeat-rich and repressive
   chromatin.
4. CpG meQTL burden and transcription/splicing coupling.
5. Shared and regionally heterogeneous architecture across donor groups and
   brain regions.
6. Schizophrenia-risk variants linked to genetically regulated VMRs.

Recommended main-figure logic:

- **Figure 1:** cohort, corrected VMR catalog, and WGBS off-array coverage;
- **Figure 2:** calibrated local SNP-explained variance and secondary held-out
  prediction;
- **Figure 3:** LINE/L1, H3K9me3, quiescent chromatin, mappability, and cell
  sensitivities;
- **Figure 4:** CpG meQTL burden and expression/splicing coupling;
- **Figure 5:** schizophrenia architecture and one or two corrected hero loci.

Cross-region, donor-group, matched, and predictor-resource details may be main
panels or supplements depending on space, but must remain available.

Use claim language such as:

- “local SNP-explained methylation variance”;
- “held-out local SNP prediction”;
- “repeat-rich and repressive genomic compartments”;
- “regionally heterogeneous LINE/L1 enrichment”;
- “transcriptionally coupled”;
- “schizophrenia-risk variants associate with methylation in genetically
  regulated VMRs.”

Always report denominators, exclusions, brain region, donor group, VMR set, and
the exact metric used.

## 12. PI decisions that must be locked

Record these in versioned configuration or a GitHub issue before production:

- primary VMR discovery cohort;
- role of the all-individual catalog;
- chromosome policy;
- VMR covariates and cell-composition strategy;
- cis-window and SNP QC;
- variance-calibration version;
- prediction screen and threshold;
- outer folds, repeats, inner folds, alpha grid, and lambda rule;
- model-release inclusion rule;
- RepeatMasker version and L1 subfamily definitions;
- primary genomic-enrichment model;
- Phase 7 success criteria and prioritized-locus rule;
- whether exposure analyses remain supplemental or are removed;
- prespecified interpretation threshold for VMR turnover.

Agents may recommend defaults but must not silently make these scientific
decisions.

## 13. Issues, branches, and handoff

- Use one GitHub issue per bounded change.
- Use one focused branch per issue; do not commit directly to `main`.
- Do not combine repository migration, estimator changes, and biological
  reinterpretation in one pull request.
- Preserve collaborator changes and do not rewrite unrelated files.

Every completed issue or pull request must state:

1. scientific question;
2. legacy and v2 paths affected;
3. files changed;
4. input and output run IDs;
5. tests executed and results;
6. sample counts and exclusions;
7. estimator, threshold, or annotation changes;
8. whether recalibration is required;
9. downstream analyses invalidated or unblocked;
10. biological claims supported, weakened, qualified, or withdrawn;
11. unresolved PI decisions or reviewer-fatal blockers.

## 14. Stop conditions

Stop production and request PI direction when:

- the primary cohort, chromosome policy, or VMR covariates are not locked;
- sample counts differ from the design without explanation;
- corrected VMR turnover exceeds the prespecified threshold;
- calibration acceptance fails or too many loci fall outside its domain;
- the prediction screen lacks validated type-I error;
- repeat enrichment disappears after locked mapping or cell-composition
  sensitivities;
- LINE/L1 conclusions differ materially across reasonable annotations;
- meQTL genomic inflation remains unresolved;
- concordance denominators are degenerate or undocumented;
- the corrected Phase 7 application fails its retention criteria;
- a proposed claim requires mediation, colocalization, cell-type specificity,
  or retrotransposition evidence that was not generated;
- computational failures remain in a production run.

Scientific validity and a coherent biological conclusion take precedence over
preserving favorable legacy numbers or a planned submission date.
