# 05_cpg_meqtl_burden — orthogonal genetic evidence

Asks whether a higher relative local SNP contribution score (`local_snp_contribution_score_z`, Module 02) is associated with a greater fraction of constituent CpGs having conventional cis-meQTL support.

**Status: implemented, smoke-verified, not yet run in production.**
The previously missing `_h/02_map_cpg_meqtl.py` (tensorqtl cis mapping) and
`_h/01b_prepare_meqtl_inputs.py` (BED / covariate / plink triple) are written,
ported from the validated legacy `meqtl-validation/01_cpg_meqtl_mapping/`
implementation. Smoke run `cmb-smoke-AA-caudate-20260823` covers chr22 only:
4,140 member CpGs -> 4,083 tested, 57 excluded by the ENCODE blacklist
(`cpg_qc.exclude_blacklist`, which the previous draft declared but ignored),
accounting balanced. 153 donors, 6 covariates (age, snpPC1-3, sex_M,
diagnosis_Schizo), 109,117 chr22 variants.

The smoke ran the full chain -- 01, 01b, 02, 02b, 04_qc_plots, 03, 04_check,
05_finalize -- and sealed with `PASS_CPG_MEQTL_BURDEN_QC` in its smoke form
(`PASS_SMOKE_ONLY_NOT_ACCEPTABLE`), all seven gate criteria passing: 4,083
prepared = 4,066 tested + 17 untested, burden fraction inside [0, 1], lambda
reported and resolved. One chromosome in one cell is not a result, and the run
is marked `smoke_run = TRUE`.

Note that `02_map_cpg_meqtl.py` deliberately computes **no** q-values:
`fdr_family: per_brain_region` means FDR is applied once across all autosomes in
`02b_combine_meqtl.R`, not 22 times per chromosome.

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

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.

## Accepted runs

_(none)_

AGENTS.md §6 makes acceptance a human step: no row appears here until a
production run's gate stage passes and the PI records it. A smoke run is never
entered in this table.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
