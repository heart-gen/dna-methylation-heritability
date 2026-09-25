# 05_cpg_meqtl_burden — convergent genetic evidence

Asks whether a higher relative local SNP contribution score (`local_snp_contribution_score_z`, Module 02) is associated with a greater fraction of constituent CpGs having conventional cis-meQTL support.

**Status: the locked covariate model is implemented, the rerun is sealed
(2026-09-24), and the three runs are ACCEPTED (2026-09-25).**
`cmb-AA-{caudate,dlpfc,hippocampus}-20260924` were mapped under the locked `M3a`
design and each sealed `PASS_CPG_MEQTL_BURDEN_QC` on 8 of 8 criteria, including
the new `executed_covariate_design_matches_lock`. They are now the runs a
downstream production run must consume, which releases the AGENTS.md §6 gate for
`07_transcription_splicing_coupling`. The `-20260825` runs remain in
the **Accepted runs** table marked **superseded**, because AGENTS.md §3 keeps a
superseded row until every consumer points at the replacement; no new downstream
production run may consume them.

See the **Accepted runs** table below. The distal-null lambda is **1.137**
caudate, **1.142** DLPFC, **1.135** hippocampus, all inside the 0.9-1.15 gate.

**Lambda did not improve genome-wide, and must not be cited as evidence for
M3a.** The chr10 decision pilot measured 1.1651 -> 1.1412 in caudate; genome-wide
the same change is 1.139 -> 1.137 in caudate and 1.139 -> 1.135 in hippocampus,
while **DLPFC moved the wrong way**, 1.137 -> 1.142. The pilot's gain was real on
chr10 and did not carry. M3a is adopted because it is the PI-locked design
(F9, 2026-09-24), which is an argument about authority and not about lambda.

**Read the "Convergent evidence, not independent replication" section below
before citing any of it** -- the title of
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
the three superseded `-20260825` runs fitted
`agedeath + sex + primarydx + snpPC1-3`. Nothing in this module or in
`00_shared/` read `config/covariates.yml`; the only reference in
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
- `_h/01b` asserts the matrix it built against the lock before writing it, and
  `_h/02` asserts the matrix it hands to tensorqtl, in the process that fits it;
- `_h/04_check_burden.R` gains the criterion
  `executed_covariate_design_matches_lock`, backed by
  `00_shared/gates.R::meqtl_covariate_design_gate()`, which reads every
  `inputs/chr*.covariates.tsv` off disk. It compares designs rather than
  spellings (`age` → `agedeath`, `sex_M` → `sex`), fails on a column that maps to
  no locked term, and fails rather than passing vacuously when there is no
  covariate file to inspect. The older `continuous_predictor_is_primary`
  criterion is the shape of check this replaces: it compares a config value with
  itself, which is what let a design diverge silently through an acceptance.

**The primary model is not free of cell composition.** `latent_factor_policy`
named a method for methPC1-5 and no specification; the PI pinned the recipe as
`primary_meqtl.latent_factor_recipe` on 2026-09-24, and every run records which
source it resolved from (see `PROPOSED_CONFIG_CHANGE.md` in this directory for
what was decided and why). The chr10 pilot measured methPC1 as 72% explained by
this region's RNA MuSiC cell proportions (R² = 0.721; Oligo ρ = +0.766,
p = 9.2e-31), and each accepted run now carries the same finding in its own
`results/latent-factor-provenance.tsv:cell_composition_note`,
so M3a carries a substantial cell-composition adjustment into the primary scan
even though `config/covariates.yml:cell_composition` reads `sensitivity_only`. The
M6d sensitivity is `M3a + dnamCellPC1-3`, so its contrast is an increment over a
baseline that already carries that structure. This is a bulk-tissue correlation
between a methylation PC and an RNA-derived proportion estimate: collinearity, not
a cell type of origin (AGENTS.md §2.3). Every run writes it into
`results/interpretation-constraints.txt`.

### What the locked model changed in the reported gradient

The burden gradient is the number Modules 07, 08 and 11 consume, so the effect of
the rerun on it is recorded here rather than left to be rediscovered. Primary
model, `quasibinomial` on the continuous `local_snp_contribution_score_z`:

| region      | estimate 20260825 | estimate 20260924 | z 20260825 | z 20260924 | n_vmrs 20260825 | n_vmrs 20260924 |
|-------------|-------------------|-------------------|------------|------------|-----------------|-----------------|
| caudate     | 2.419             | **2.084**         | 39.45      | **42.70**  | 11,231          | 11,142          |
| dlpfc       | 2.449             | **2.094**         | 34.95      | **41.16**  | 9,216           | 9,134           |
| hippocampus | 2.495             | **2.367**         | 36.57      | **42.63**  | 9,140           | 9,053           |

The **estimate attenuates** in all three regions and the **z rises** in all
three: the gradient is smaller and better determined. The direction of the claim
is unchanged and so is the gate.

**Do not attribute the attenuation to the covariate model.** Two inputs changed
between these two dates -- the covariate design *and* the upstream Module 02
score run, see below -- so this is a two-factor comparison with one cell observed. The
table records what the accepted numbers are; it does not decompose why they moved,
and no arm was run that would. In particular nothing here licenses a statement
about how much of the gradient is cell composition: methPC1 is a methylation PC
correlated with an RNA-derived proportion estimate, the attenuation is not
decomposed, and a cell type of origin does not follow (AGENTS.md §2.3).

**The rerun also corrected the upstream score pointer, and that is a second,
independent defect.** The `-20260825` runs consumed `lgv-AA-{region}-20260823`,
which Module 02 superseded on 2026-09-17 for scoring under the pooled `p_eff`
floor 2.058 rather than AA's own 7.079 -- it over-admitted 90 caudate, 90 DLPFC and
99 hippocampus loci that AA's own support excludes. The `-20260924` runs consume
the accepted `lgv-AA-{region}-rescore-20260913`. That is the whole of the `n_vmrs`
change (89, 82 and 87 fewer VMRs), and it means the `-20260825` runs were already
superseded on their upstream pointer alone, before F9 is considered. The CpG
denominators are untouched: `prepared`/`tested`/`untested` are identical across the
two dates in all three regions, and `unaccounted = 0` in every case.

Significant CpGs rose with the locked model, 99,203 -> 105,553 in caudate at the
same FDR threshold and the same 189,998-CpG denominator, which is the power gain
the latent factors buy.

### The genotype QC donor set was the source pfile, not the analysis set

`_h/01b` wrote a `{chrom}.keep` file and never passed it to plink2, so `--maf`,
`--geno` and `--hwe` were evaluated over all 526 AA donors in the source pfile
rather than the region's analysis donors — 373 of them outside a 153-donor
caudate estimation set. `00_shared/locus_io.R` documents the opposite convention
for Modules 02 and 03, where the filters run *after* the group restriction
deliberately, so this was Module 05's divergence from project practice.

The keep file was also malformed — `{donor}\t{donor}`, the FID twice — while the
psam's IID is a chip barcode (`Br2585` → `3998646007_R01C01`). plink2 reads a
two-column `--keep` as FID/IID, so the file matched **zero** samples; had it been
passed as written, every array task would have failed. That is why the missing
flag went unnoticed.

Measured on chr10 (pilot, and reproduced by `tests/covariate_lock_smoke.py`): the
locked QC over the 153 donors leaves 414,363 variants against the executed run's
416,659, with 17,832 over-included and 15,536 wrongly excluded. 90.2% of the
over-included fall to tensorqtl's in-sample MAF filter, but the surviving 1,747
were tested and **all** of them exceed the locked `missingness_max: 0.05` in the
153 donors, because tensorqtl re-applies MAF and never missingness. Lambda barely
moves (1.1651 → 1.1619), so this is a denominator and QC-compliance defect rather
than a calibration one — not a reason to rerun by itself, and a fix any rerun must
carry.

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

**Current.** Accepted 2026-09-25 by the PI. These are the runs a downstream
production run must consume.

Acceptance was held until the cross-region sample-integrity screen was
adjudicated, because that was the only open item with the reach to invalidate
these runs: excluding a donor would have re-derived Module 01, changed
`vmr_set_id`, and taken every module with it. It was closed on 2026-09-25 with no
donor excluded, so the donor set behind these three runs is final --- see
`11_integrated_manuscript_outputs/T4_READMITTED_DONOR_ADJUDICATION.md` and the
affirmation in the `sample_blacklist` comment block of `config/cohorts.yml`, both
of which land on branch `qc/t4-readmitted-donor-adjudication`. Every other open
finding routes through a module this one does not consume, or is metadata that no
run reads back.

| run_id                      | cohort | region      | vmr_set_id                         | upstream_lgv_run_id                 | sealed_at            | accepted_on | accepted_by | decision                 | notes                                                                                    |
|-----------------------------|--------|-------------|------------------------------------|-------------------------------------|----------------------|-------------|-------------|--------------------------|------------------------------------------------------------------------------------------|
| cmb-AA-caudate-20260924     | AA     | caudate     | vmrset-AA-caudate-937a41979978     | lgv-AA-caudate-rescore-20260913     | 2026-09-24T11:20:05  | 2026-09-25  | Kynon J.M. Benjamin | PASS_CPG_MEQTL_BURDEN_QC | 8/8 criteria; locked M3a verified on 22 chromosome files; distal-null lambda 1.137; n_vmrs 11,142 |
| cmb-AA-dlpfc-20260924       | AA     | dlpfc       | vmrset-AA-dlpfc-856067dfe289       | lgv-AA-dlpfc-rescore-20260913       | 2026-09-24T11:33:49  | 2026-09-25  | Kynon J.M. Benjamin | PASS_CPG_MEQTL_BURDEN_QC | 8/8 criteria; locked M3a verified; distal-null lambda 1.142 -- the one region where lambda ROSE; n_vmrs 9,134 |
| cmb-AA-hippocampus-20260924 | AA     | hippocampus | vmrset-AA-hippocampus-2d907b892215 | lgv-AA-hippocampus-rescore-20260913 | 2026-09-24T11:54:46  | 2026-09-25  | Kynon J.M. Benjamin | PASS_CPG_MEQTL_BURDEN_QC | 8/8 criteria; locked M3a verified; distal-null lambda 1.135; n_vmrs 9,053                 |

All three are convergent evidence, not independent replication, and all three
carry the cell-composition constraint in
`results/interpretation-constraints.txt`. All three record `git_dirty = true`
against commit `8d12d3900`, so they are not byte-reproducible from the recorded
commit alone -- a provenance defect of the same class as F10, recorded rather than
hidden, and not a data defect: every config checksum is intact and the covariate
design is verified off disk by the gate.

**Superseded.** Retained per AGENTS.md §3 until every consumer points at the
replacement. **No new downstream production run may consume these.**

| run_id                      | cohort | region      | vmr_set_id                         | upstream_lgv_run_id           | accepted_on | superseded_by               | why superseded                                                                              |
|-----------------------------|--------|-------------|------------------------------------|-------------------------------|-------------|-----------------------------|---------------------------------------------------------------------------------------------|
| cmb-AA-caudate-20260825     | AA     | caudate     | vmrset-AA-caudate-937a41979978     | lgv-AA-caudate-20260823       | 2026-08-28  | cmb-AA-caudate-20260924     | Two defects: fitted snpPC1-3 with no methPC against the locked M3a (F9); and consumed a Module 02 run superseded 2026-09-17 (90 over-admitted loci) |
| cmb-AA-dlpfc-20260825       | AA     | dlpfc       | vmrset-AA-dlpfc-856067dfe289       | lgv-AA-dlpfc-20260823         | 2026-08-28  | cmb-AA-dlpfc-20260924       | Same two defects; 90 over-admitted loci upstream                                            |
| cmb-AA-hippocampus-20260825 | AA     | hippocampus | vmrset-AA-hippocampus-2d907b892215 | lgv-AA-hippocampus-20260823   | 2026-08-28  | cmb-AA-hippocampus-20260924 | Same two defects; 99 over-admitted loci upstream                                            |

**Not eligible.** Both chr22 smokes are `smoke_run = TRUE` and neither is
citable. `cmb-AA-caudate-20260924-smoke` **never sealed** -- it has no `decision`
and no `sealed_at`, because the chain failed; it is kept as evidence of the
failure and correctly has no decision row. `cmb-AA-caudate-20260924-smoke-b`
sealed `PASS_SMOKE_ONLY_NOT_ACCEPTABLE` at 2026-09-24T10:15:00 and is the passing
chr22 exercise that preceded the three production chains.
