# 05_cpg_meqtl_burden — convergent genetic evidence

Asks whether a higher relative local SNP contribution score (`local_snp_contribution_score_z`, Module 02) is associated with a greater fraction of constituent CpGs having conventional cis-meQTL support.

**Status: the three accepted runs are superseded and a rerun is required
(2026-09-24).** They were mapped under a covariate design that is not the locked
one; see "The executed covariate model diverged from the lock" below. Their rows
stay in the **Accepted runs** table until replacements are accepted, because
nothing downstream has a v2 replacement yet, but no new downstream production run
should consume them.

See the **Accepted runs** table below for `cmb-AA-{caudate,dlpfc,hippocampus}-20260825`,
each `PASS_CPG_MEQTL_BURDEN_QC`. The caudate run carries a distal-null lambda of
1.139, recorded in its acceptance note. **Read the "Convergent evidence, not
independent replication" section below before citing any of it** -- the title of
this module says *convergent*, not *orthogonal*, and that is a substantive
claim about what the analysis can support, not a wording preference.

`_h/02_map_cpg_meqtl.py` (tensorqtl cis mapping) and `_h/01b_prepare_meqtl_inputs.py`
(BED / covariate / plink triple) were ported from the validated legacy
`meqtl-validation/01_cpg_meqtl_mapping/` implementation. The code was first
exercised on smoke run `cmb-smoke-AA-caudate-20260823`, chr22 only: 4,140 member
CpGs -> 4,083 tested, 57 excluded by the ENCODE blacklist
(`cpg_qc.exclude_blacklist`, which the previous draft declared but ignored),
accounting balanced. That run is marked `smoke_run = TRUE` and is not eligible
for the accepted-runs table.

Of the nine QC plots `config/meqtl_parameters.yml` requires, eight are built.
`covariate_model_comparison` is **deliberately skipped, not missed**: it needs a
second mapping run under an alternative covariate model, which has not been
commissioned. `_h/04_qc_plots.py` emits an explicit skip record rather than
quietly producing eight of nine.

Note that `02_map_cpg_meqtl.py` deliberately computes **no** q-values:
`fdr_family: per_brain_region` means FDR is applied once across all autosomes in
`02b_combine_meqtl.R`, not 22 times per chromosome.

### The executed covariate model diverged from the lock

`config/covariates.yml:primary_meqtl` has been locked since 2026-08-01 to

```text
M3a = agedeath + sex + primarydx + snpPC1-5 + methPC1-5
```

`_h/01b_prepare_meqtl_inputs.py` set `n_pc = 3` and added no methylation PC, so
the three accepted runs fitted `agedeath + sex + primarydx + snpPC1-3`. Nothing in
this module or in `00_shared/` read `config/covariates.yml`; the only reference in
the repository was `00_shared/runid.R`, which writes `config_covariates_sha256`
into the manifest. The lock was attested by every sealed run and enforced by
nothing.

**PI decision 2026-09-24: the lock is authoritative.** The evidence is the
decision pilot on branch `module05/covariate-model-pilot`
(`PILOT_COVARIATE_MODEL.md`): on chr10, caudate, M3a lowers the distal-null lambda
from 1.1651 to 1.1412, and every pi0-free discovery threshold favours it. The gain
is entirely the latent factors — an snpPC1-5-only arm gives 1.1701, marginally
*worse* than snpPC1-3 — so snpPC4-5 contribute nothing on their own.

What changed here:

- `_h/01a_estimate_latent_factors.py` (new) estimates methPC1-15 once per run,
  before the mapping array, and writes `results/latent-factor-provenance.tsv`;
- `_h/meqtl_covariates.py` (new) builds the design and
  `00_shared/covariate_lock.py` (new) expands the lock from config, so no stage
  types a term list;
- `_h/01b` asserts the matrix it built against the lock before writing it.

**The primary model is not free of cell composition.** The lock names a method for
methPC1-5 and no specification, so the recipe is resolved explicitly in code and
recorded per run (see `PROPOSED_CONFIG_CHANGE.md` at the repository root for the
config block that would pin it). The pilot found methPC1 is 72% explained by this
region's RNA MuSiC cell proportions (R² = 0.721; Oligo ρ = +0.766, p = 9.2e-31),
so M3a carries a substantial cell-composition adjustment into the primary scan
even though `config/covariates.yml:cell_composition` reads `sensitivity_only`. The
M6d sensitivity is `M3a + dnamCellPC1-3`, so its contrast is an increment over a
baseline that already carries that structure. This is a bulk-tissue correlation
between a methylation PC and an RNA-derived proportion estimate: collinearity, not
a cell type of origin (AGENTS.md §2.3). Every run writes it into
`results/interpretation-constraints.txt`.

### Genomic inflation is measured on distal cis pairs

`config/meqtl_parameters.yml:genomic_inflation` gates on lambda computed over
nominal cis pairs with |CpG-to-SNP distance| > 400 kb, not over all cis pairs.
Lambda over all cis pairs is not an inflation estimate: cis pairs are enriched
for true meQTLs, so a healthy scan has lambda > 1 by construction. The chr22
smoke shows this directly -- lambda falls monotonically with distance (1.234 all
pairs, 1.146 > 100 kb, 1.113 > 250 kb, 1.089 > 400 kb), the signature of real
signal rather than global inflation. The all-pair figure is still recorded as
`lambda_all_pairs` and the whole decay profile as
`qc/genomic-inflation-by-distance.tsv`. Because the gate depends on the nominal
pass, `04_qc_plots.py` runs **before** `03_vmr_burden.R` in the job graph.

### Reading Module 01's CpG matrices

`cpg/chr_N/cpg_meth.phen` is donors x every CpG on the chromosome -- chr1 is
5.2 GB and about two million columns. `fread` cannot open a file that wide at
all: it segfaults ("memory not mapped") while setting up per-column state,
before parsing a single row, and `select=` does not help because the crash
happens first. `01_prepare_cpg_set.R` therefore reads the header alone, works
out which columns are VMR members, extracts just those with `cut` (which
streams and does not care how wide a line is), and hands only the narrow file
to `fread`. chr1 completes in under two minutes at ~6 GB peak RSS -- almost all
of it the two-million-name header.

Because stage 01 runs on the submit host, a whole-genome production run spends
roughly 45 minutes there before the array is submitted. Worth moving into a
batch job if that becomes a nuisance.

### Known gaps

- (Resolved 2026-08-25.) `positive_control_public_meqtl_overlap` is implemented
  in `04_qc_plots.py`. `jaffe_dlpfc_450k_meqtl` and
  `schulz_hippocampus_array_meqtl`, harmonized to hg38, are registered in
  `inputs/supportfiles/_m/annotation_asset_manifest.tsv`; the plot discovers
  them by their notes field, so registering a further catalog needs no code
  change.

  Independence is deliberately NOT required here. A positive control asks only
  whether the scan recovers already-known meQTLs, so the Phase 3 exclusions do
  not apply and cohort overlap would if anything strengthen it. The two Phase 3
  rules that DO carry over are the assayed universe as denominator (the
  harmonized tables embed all 485,441 probes with an `external_assayed` flag)
  and an exact hg38 `(chrom, pos_1based)` join.

  The statistic is the recovery rate among externally supported CpGs against
  the rate among assayed-but-unsupported ones; the second column is what stops
  a scan that calls everything significant from scoring a perfect control. On
  the chr22 smoke (caudate, both references cross-tissue): Jaffe 0.615 vs 0.222
  (OR 5.58, p = 6.5e-07), Schulz 0.783 vs 0.423 (OR 4.91, p = 0.0015), on 172
  shared assayed CpGs.

  Caveats: only 4.2% of the scan's CpGs are on 450K, and **caudate has no
  tissue-matched public resource**, so its control is cross-tissue in both
  cases. `tissue_match` is recorded per row rather than used as a filter.

- `covariate_model_comparison` is still not implemented: it needs a second
  mapping run under an alternative covariate design, and which design to use is
  a PI decision. It remains recorded per run in
  `results/qc/qc-plots-not-produced.tsv`.
- (Resolved 2026-08-25.) Neither Bioconductor `qvalue` nor `py_qvalue` is in
  the `epigenomics` environment, but `py_qvalue` 0.1.0 **is** in `genomics`.
  `02b_combine_meqtl.R` now tries Bioconductor `qvalue`, then
  `_h/storey_qvalue.py` under `genomics`, then BH, recording which it used as
  `fdr_method_used`. No package was added to a shared environment.

  `storey_qvalue.py` must pass `lfdr_out=False`: py_qvalue's local-FDR branch
  does not return in any reasonable time even on 4,066 p-values. The helper
  also re-checks monotonicity and the pi0 bound before returning, so the R
  stage never has to trust an unvalidated 0.1.0 dependency.

  This was not cosmetic. On the chr22 smoke pi0 = 0.217, and Storey calls 2,098
  CpGs at FDR 0.05 against BH's 1,636 -- 28% more, in the numerator of the
  module's endpoint.

## Migrating from

`meqtl-validation/01_cpg_meqtl_mapping/` and `meqtl-validation/02_vmr_meqtl_burden/`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Requirements

Use the corrected CpG-to-VMR membership keyed by the accepted `vmr_set_id`.
Report tested CpGs separately from prepared-but-untested CpGs. Model the
continuous standardized local SNP contribution score as the primary predictor,
with overdispersion-appropriate and donor-robust inference. Top-versus-bottom
quartile analysis is secondary relative evidence only; quartile boundaries are
not biological cutoffs. Audit every concordance denominator, and resolve
genomic inflation before the figure freeze.

Internal meQTL mapping is convergent evidence, not independent replication.
Positive-only public resources cannot provide an external gradient, because
absence from a positive list is not a tested negative.

## Convergent evidence, not independent replication

The primary model regresses per-VMR meQTL burden on
`local_snp_contribution_score_z`. Both sides of that regression are derived
from the genotype-methylation covariance in the **same donors** (AA: n = 153
caudate, 118 DLPFC, 117 hippocampus) and from the same accepted Module 01 VMR
catalog. The predictor comes from Module 02's variance estimate on those
donors; the outcome is a cis-meQTL scan on those same donors' genotypes and
methylation. The two quantities are therefore not independent measurements of
local genetic control -- they are two summaries of one covariance structure.

A large, highly significant coefficient is expected by construction and is
**not** evidence of external validity. The production runs give
`local_snp_contribution_score_z` beta = 2.42 / 2.45 / 2.49 (caudate / DLPFC /
hippocampus) at z = 39.4 / 35.0 / 36.6. The defensible reading of that result
is internal consistency: three regions, mapped independently of one another,
agree on the effect size to within 3%. The magnitude of z is a property of the
design, not a finding, and must not be reported as one.

What this module therefore does and does not license:

- **Does**: confirm that the Module 02 score tracks a directly observable
  molecular quantity (the fraction of a VMR's CpGs carrying a significant cis
  meQTL) in the expected direction, consistently across regions.
- **Does**: support the relative, rank-based use of the score, which is the
  prespecified endpoint.
- **Does not**: replicate Module 02 in an independent sample.
- **Does not**: calibrate or validate any absolute PVE, effect size, or
  threshold.

Independent support comes only from the public-resource positive control
(`results/qc/positive_control_public_meqtl_overlap.tsv`), which is external to
these donors. That control is directional evidence with a known ceiling:
positive-only resources cannot supply a tested negative, so a CpG's absence
from a public list is not evidence against an meQTL there, and the recovery
contrast it reports is a floor rather than an estimate.

## Contract

This module follows: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.

## Accepted runs

| run_id                      | cohort | region      | vmr_set_id                         | accepted_on | accepted_by         | decision                 | notes                                                                      |
|-----------------------------|--------|-------------|------------------------------------|-------------|---------------------|--------------------------|----------------------------------------------------------------------------|
| cmb-AA-caudate-20260825     | AA     | caudate     | vmrset-AA-caudate-937a41979978     | 2026-08-28  | Kynon J.M. Benjamin | PASS_CPG_MEQTL_BURDEN_QC | Convergent evidence, not independent replication; distal-null lambda 1.139 |
| cmb-AA-dlpfc-20260825       | AA     | dlpfc       | vmrset-AA-dlpfc-856067dfe289       | 2026-08-28  | Kynon J.M. Benjamin | PASS_CPG_MEQTL_BURDEN_QC | Convergent evidence, not independent replication                           |
| cmb-AA-hippocampus-20260825 | AA     | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-08-28  | Kynon J.M. Benjamin | PASS_CPG_MEQTL_BURDEN_QC | Convergent evidence, not independent replication                           |
