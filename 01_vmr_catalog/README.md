# 01_vmr_catalog — corrected biological units

Defines the variably methylated regions (VMRs) every downstream v2 module is
built on, for two cohort arms across three brain regions, from one parameterized
codepath.

**Status: accepted.** All six arm x region production runs pass all five
acceptance criteria as of 2026-08-16; see *Acceptance gate* below for the run
IDs and `vmr_set_id`s downstream modules must cite (AGENTS.md 6).

## Why this module exists

The VMR sets in `vmr-analysis/` are invalid. From
`writing-notes/PIPELINE_AUDIT.md`:

> **V1 (critical).** `02b.res_var.R::filter_pheno()` reorders the methylation
> matrix with `match(valid_ids, brain_id)` but subsets the PC design matrix with
> `pc[V1 %in% valid_ids, ]`, which preserves `pc.csv` order. Response rows and
> design rows are different donors — caudate chr1 has 153/153 rows wrong, DLPFC
> 92/96, hippocampus 98/101.

The consequence is not cosmetic: the 99th-percentile SD cutoff that *defines* a
VMR shifted by +9.2%, and roughly 10% of seed CpGs changed (Jaccard 0.82). Every
VMR set, and therefore every result conditioned on VMR membership or boundaries,
has to be recomputed (AGENTS.md §6: "VMR turnover never authorizes reuse of
downstream numbers").

## Cohort arms

| Arm | Donors | Role |
|---|---|---|
| `AA` | Black American | **primary** |
| `all_individuals` | Black + white American | sensitivity |

Both arms run for `caudate`, `dlpfc`, and `hippocampus`. Definitions live in
`config/cohorts.yml`; there is one implementation, not six copies.

## Pipeline

| Step | Script | Scope | What it does |
|---|---|---|---|
| 0 | `00_new_run.R` | per arm × region | Mints the run ID, writes the provenance manifest |
| 1 | `00_prepare.R` | per chromosome | Donor selection, C→T masking, coverage QC, CpG matrix, covariates |
| 2 | `01_analyze.R` | per chromosome | Methylation PCs, residualization, per-CpG residual variance |
| 3 | `02_summarize.R` | per arm × region | Calls VMRs, per-VMR methylation, mints `vmr_set_id` |
| 3 | `03_plot.R` | per arm × region | Catalog diagnostics |
| 3 | `04_turnover.R` | per arm × region | Old-vs-new turnover, array coverage, QC/exclusion tables |
| 4 | `step_4.sh` | per VMR | Cis genotype window extraction |
| 5 | `05_close_run.R` | per arm × region | Checksums outputs, seals the run read-only |

Launchers: `step_1.sh` (autosomes), `step_1x.sh` (X/Y → `excluded/`),
`step_2.sh`, `step_3.sh`, `step_4.sh`, `step_5.sh`, and
`submit_vmr_catalog_workflow.sh` to chain them.

`step_5.sh` must be **last** — it makes the run directory read-only, so anything
that still needs to write (notably `step_4.sh`) has to finish first.

```bash
cd 01_vmr_catalog/_m && mkdir -p logs
../_h/submit_vmr_catalog_workflow.sh AA caudate
```

## Repairs applied

| ID | Defect | Fix | Where |
|---|---|---|---|
| V1 | PC design rows misaligned against methylation rows | `align_by_id()` reorders both sides and asserts | `00_shared/identity.R`, `01_analyze.R` |
| V2 | `region == "caudate"` hard-coded in all three BA_only copies | Region filter comes from `--region` | `00_prepare.R` |
| V3 | `here()` roots pointing at the renamed `heritability/` tree | All paths via `config/paths.yml` | `00_shared/config.R` |
| V4 | Unmasked X/Y VMRs merged into the primary catalog | Autosomes primary; X/Y to `excluded/` with a manifest | `config/thresholds.yml`, `00_prepare.R` |
| V5 | `remove_ct_snps()` read a global, ignoring its argument | Operates on its argument | `00_prepare.R` |
| V6 | `v_top[1:10^6, ]` fabricates NA rows on small chromosomes | `min(cap, nrow(v))` | `01_analyze.R` |
| V7 | Chunks combined in lexicographic order, donor order unchecked | `sort_chunks()` + `verify_fid_iid()` | `00_shared/chrom.R` |
| V8 | Headerless AA `.psam` read with `header=TRUE`, dropping Br2585 | `read_psam()` handles both forms | `00_shared/identity.R` |
| V9 | Aborted on telomere-proximal VMRs, dropping ~1% | Window clamps to 1 | `step_4.sh` |
| V10 | chr1's length used as the bound for every chromosome | `chrom_size()` per chromosome | `00_shared/slurm.sh` |
| V11 | 500 kb vs 1 Mb cis window between arms | One window from `config/thresholds.yml` | `step_4.sh` |
| V12 | `module load plink` | plink2 from `/projects/p32505/opt/bin/plink2` | `00_shared/slurm.sh` |
| — | Hippocampus read the **DLPFC** sample blacklist | Blacklists retired; per-region resolution retained, no fallback | `config/cohorts.yml`, `00_shared/config.R` |

### Sample blacklists are retired

The legacy `samples_blacklist.txt` files were not a QC exclusion. They existed
only to reconcile a stale AA phenotype file that was missing those donors, so
the all-individual arm could be matched to the Black American-only arm. v2 reads
one phenotype table for both arms, so nothing needs reconciling and all eight
donors (7 DLPFC, 1 hippocampus — all AA, all with usable WGBS and genotypes) are
**included**.

Consequence for the acceptance gate: the AA analysis sets should be **153 /
118 / 117** for caudate / DLPFC / hippocampus, not the 111 and 116 recorded in
`vmr-analysis/{dlpfc,hippocampus}/_m/samples.txt`, which are post-blacklist.

## Outputs

Per accepted run, under `_m/runs/{RUN_ID}/`:

- `vmr/vmr.bed`, `vmr/vmr_catalog.tsv` — the catalog, carrying `vmr_set_id`
- `vmr/cpg_vmr_membership.tsv` — CpG → VMR membership (consumed by 05)
- `vmr/phenotypes/{chr}_{start}_{end}_meth.phen` — per-VMR methylation
- `vmr/donors_plink.txt` — `--keep` list, in catalog donor order
- `vmr/sd_cutoffs.tsv` — per-chromosome SD cutoff and counts
- `cpg/chr_{N}/`, `covs/chr_{N}/`, `pca/chr_{N}/` — per-chromosome intermediates
- `qc/vmr_turnover.tsv`, `qc/array_coverage.tsv`, `qc/technical_qc.tsv`,
  `qc/exclusions.tsv`
- `excluded/` — sex chromosomes, with the reason recorded
- `manifest.tsv`, `task_reconciliation.tsv`, `output_checksums.tsv`

## Acceptance gate

A run is production only when all of these hold, recorded in the table below:

1. Observed donor counts match `config/cohorts.yml` `design_n` for all six
   arm × region cells. (`design_n` is `null` until the first run establishes it
   and the PI locks it.)
2. `task_reconciliation.tsv` reports zero unaccounted and zero failed tasks.
3. Two invocations of the same smoke configuration produce identical output
   checksums.
4. Turnover against the legacy catalog is reported and reviewed. Substantial
   turnover is *expected* — it is the V1 repair — but a catastrophic difference
   means a new bug.
5. `config/*.yml` carry `pi_locked: true` for every key the run consumes.

| Run ID | Arm | Region | n | VMRs | vmr_set_id | Gate | Notes |
|---|---|---|---|---|---|---|---|
| `vmrcat-AA-caudate-20260815-a` | AA | caudate | 153 | 260 | `vmrset-AA-caudate-3968cebd6d9c` | **smoke, not production** | chr22 only, `--allow-unlocked`; data deleted 2026-08-16, manifest kept in `_m/runs/retired/` |
| `vmrcat-AA-caudate-20260816` | AA | caudate | **153** | 11,530 | `vmrset-AA-caudate-937a41979978` | **all five pass** | all 22 autosomes |
| `vmrcat-AA-dlpfc-20260816` | AA | dlpfc | **118** | 9,572 | `vmrset-AA-dlpfc-856067dfe289` | **all five pass** | first full-scale run, all 22 autosomes |
| `vmrcat-AA-hippocampus-20260816` | AA | hippocampus | **117** | 9,497 | `vmrset-AA-hippocampus-2d907b892215` | **all five pass** | all 22 autosomes |
| `vmrcat-all_individuals-caudate-20260816` | all_individuals | caudate | **282** | 11,463 | `vmrset-all_individuals-caudate-cb5519d7d2ad` | **all five pass** | sensitivity arm |
| `vmrcat-all_individuals-dlpfc-20260816` | all_individuals | dlpfc | **173** | 9,374 | `vmrset-all_individuals-dlpfc-e88f46904afb` | **all five pass** | sensitivity arm; turnover outlier, see forensic audit below |
| `vmrcat-all_individuals-hippocampus-20260816` | all_individuals | hippocampus | **177** | 9,365 | `vmrset-all_individuals-hippocampus-809f8de0db2d` | **all five pass** | sensitivity arm |

All six runs seal with `smoke_run FALSE` and reconcile at 22 expected / 22
completed / 0 excluded / 0 qc_failed / 0 failed / 0 unaccounted / 0 unexpected,
which clears criterion 2 for every cell. `design_n` is locked for all six cells
in `config/cohorts.yml` as of 2026-08-16, so criterion 1 is met: every observed
count matches. Criterion 4 is met for five cells directly and for
`all_individuals` dlpfc through the forensic audit below. Criterion 3 is met by the paired smoke runs recorded below.

The `config_cohorts_sha256` staleness noted below for `vmrcat-AA-dlpfc-20260816`
applies to all six runs for the same reason, and now to the `all_individuals`
lock as well.

### Criterion 3: checksum stability

Two independent invocations of the same smoke configuration, submitted
separately on 2026-08-16: `vmrcat-AA-caudate-20260816-a` and
`vmrcat-AA-caudate-20260816-b`, both `VMR_CHROMS=21,22`, `ALLOW_UNLOCKED=1`,
`smoke_run TRUE`. Neither is production; they exist only to answer this gate.

Both produced the **same `vmr_set_id`, `vmrset-AA-caudate-ae953670d7eb`**, 367
VMRs, 7,278 CpGs in VMRs, 153 donors, the same `donor_checksum`, and 410 output
files each.

Of the 410 files, **408 are byte-identical**. The two that differ are
`pca/chr_21/pca.pdf` and `pca/chr_22/pca.pdf`, which differ only in the
`/CreationDate` and `/ModDate` strings the PDF device embeds at write time --
identical file sizes, identical content otherwise. That is a property of the
graphics device, not of the computation.

The `donor_checksum` also matches the retired chr22-only smoke run from
2026-08-15 and the production `vmrcat-AA-caudate-20260816`, so donor resolution
is now reproducible across three independent invocations, two chromosome
scopes, and two days.

**Criterion 3 passes.** The catalog itself -- VMR boundaries, membership,
phenotypes, covariates, and the `vmr_set_id` derived from them -- is
bit-reproducible.

### Cis-variant density and the legacy 1 Mb window

Per-VMR variant counts in the 500 kb cis window, from each run's
`plink_format/extraction_log/`:

| Cell | VMRs | 0 variants | 1–99 | <100 total | median variants |
|---|---|---|---|---|---|
| AA caudate | 11,530 | 141 (1.22%) | 25 (0.22%) | **1.44%** | 5,858 |
| AA dlpfc | 9,572 | 176 (1.84%) | 19 (0.20%) | **2.04%** | 5,764 |
| AA hippocampus | 9,497 | 178 (1.87%) | 15 (0.16%) | **2.03%** | 5,799 |
| all_individuals caudate | 11,463 | 153 (1.33%) | 57 (0.50%) | **1.83%** | 2,643 |
| all_individuals dlpfc | 9,374 | 181 (1.93%) | 45 (0.48%) | **2.41%** | 2,607 |
| all_individuals hippocampus | 9,365 | 178 (1.90%) | 43 (0.46%) | **2.36%** | 2,611 |

**500 kb is not variant-limiting.** The median VMR carries thousands of cis
variants, and the excluded fraction is 1.4–2.4% in every cell. Widening to 1 Mb
would not recover most of it: the great majority of excluded loci have *zero*
variants, not few, and they sit in pericentromeric heterochromatin and
acrocentric short arms (chr1, chr9, chr20, chr21, chr15 account for ~90% of
them) where the imputation panel has no coverage at any window size. Doubling
the window there doubles an empty interval.

**The arm asymmetry is a genotype-file property, not a window problem.** The
`1–99 variant` count is ~2× higher in `all_individuals` than in the matching AA
region. The cause is variant density in the source files:

    inputs/genotypes/TOPMed_LIBD.AA                    15,865,072 variants / 526 samples
    inputs/genotypes/all_individuals/TOPMed_LIBD        7,109,096 variants / 1,938 samples

The combined file carries 3.7× the samples but 2.2× *fewer* variants, and the
median cis window is correspondingly thinner (2,643 vs 5,858 in caudate). The
two files were evidently filtered under different MAF/INFO rules. This is the
most likely reason the legacy all-individuals arm used a 1 Mb window while
`BA_only` used 500 kb (V11): the wider window compensated for a sparser panel.
It made the two arms incomparable, which is why v2 fixes the window at 500 kb
for both and lets the density difference show up honestly as a variant count.

The density gap is worth carrying into `02_local_genetic_variance`: an
AA-versus-all_individuals difference in local genetic variance is partly a
difference in how many variants each arm's elastic net can see. It is not a
module 01 gate.

### `vmrcat-AA-dlpfc-20260816`

The first full-scale run of this module. Against the acceptance gate:

1. **Donor count — 118, as expected.** Not the legacy 96, and not the 111 in
   `vmr-analysis/dlpfc/_m/samples.txt`. This is the observable signature of the
   V2 repair plus the retired blacklist, and it is now *demonstrated* rather than
   reconstructed from the inputs. `design_n` in `config/cohorts.yml` is still
   `null` pending the PI lock, so this criterion is not formally met.
2. **Reconciliation clean.** 22 expected, 22 completed, 0 excluded, 0 qc_failed,
   0 failed, 0 unaccounted, 0 unexpected. `step_4`: 9,572/9,572 completed.
3. **Checksum stability across two invocations — NOT DONE.**
4. **Turnover reported.** 9,572 v2 vs 10,229 legacy autosomal VMRs (plus 143
   legacy sex-chromosome VMRs dropped by V4); 5,131 overlapping, 4,441 novel,
   5,066 lost; Jaccard 0.315; **46.4% of v2 VMRs novel**, within the
   `max_vmr_turnover: 0.95` gate. This is well above the 31% seen in the chr22
   caudate smoke run, which is the expected direction: DLPFC is the region where
   V1 misaligned 92 of 96 donors, so it had the most to correct.
5. Configs locked; `smoke_run FALSE`. `git_dirty true` (untracked `AGENTS.md`).

**`config_cohorts_sha256` in the sealed manifest is stale by design.** The
manifest records `e4797e0b…`, the hash of `config/cohorts.yml` *before*
`design_n` was locked; the file now hashes to `7dec8a52…`. The locked values are
a record of what these runs observed, so writing them necessarily changes the
file after the runs that produced them. Nothing about the computation differs —
`design_n: null` makes `assert_expected_n()` warn, a set value makes it stop —
and re-deriving these catalogs requires the *pre-lock* config, which is what the
manifest points at. Runs submitted after 2026-08-16 record the locked hash.

**Cis genotype extraction (`step_4`, first execution ever).** 9,396 VMRs
extracted, median 5,764 variants; 106 windows clamped at a chromosome start and
133 at an end — V9 and V10 both firing on real data, where the legacy code would
have aborted on those 106.

**176 VMRs (1.8%) have no cis variants** and are recorded as
`no_cis_variants`: chr9 54, chr1 49, chr21 30, chr20 27, chr15 7, chr22 7,
chr5 1, chr14 1. These fall in pericentromeric heterochromatin and acrocentric
short arms, where imputation panels carry no variants. A further 19 VMRs have
fewer than 100 cis variants and are excluded by
`config/thresholds.yml: cis.min_cis_variants`.

Both exclusions track genomic features rather than occurring at random, so the
analyzed set is **not** a uniform sample of the catalog. `02_local_genetic_variance`
records them as QC failures rather than dropping them silently, and the methods
report the counts per region.

### Forensic audit: `vmrcat-all_individuals-dlpfc-20260816`

Forensic audit performed 2026-08-16. This run produced the largest
old-versus-new discrepancy of the six full
arm-by-region catalogs: 9,374 v2 VMRs versus 9,859 legacy autosomal VMRs, with
2,607 v2 intervals overlapping a legacy interval, 6,767 novel v2 intervals,
7,229 lost legacy intervals, base-pair Jaccard 0.131, and 72.2% of v2 VMRs
novel. The result is an unusually severe manifestation of V1 in the legacy
DLPFC comparator, not evidence that the all-individuals sensitivity arm or the
v2 DLPFC catalog is biologically unstable.

The legacy all-individuals DLPFC residualization reordered the methylation
matrix to phenotype-table donor order but retained `pc.csv` row order for the
methylation-PC design. Only 17 of 166 row positions matched by accident; 149
donors were paired with another donor's PC values. The AA DLPFC branch carried
the same defect through a different donor ordering and a separately estimated
PC matrix, so the two legacy branches did not receive the same permutation even
though they shared most donors.

Four diagnostics localize the discrepancy to this legacy residualization:

1. **Realignment recovers v2.** In 5,000-CpG blocks from chromosomes 1, 2, 6,
   11, 16, and 22, reordering the frozen legacy PC rows by donor ID before
   residualization yielded Pearson correlations of 0.9978--0.9996 with the v2
   residual SDs. The as-coded legacy SDs correlated only 0.8774--0.9486 with
   v2. This sensitivity retained the legacy 166-donor set, so the seven donors
   restored in v2 cannot explain the turnover.
2. **The effect is chromosome-wide.** Per-chromosome VMR base-pair Jaccard
   ranged from 0.066 to 0.221 across all 22 autosomes. No single failed or
   mismatched chromosome drives the aggregate result.
3. **DLPFC was unusually sensitive to permuting the design.** Across
   chromosomes, the first five legacy methylation PCs explained a median 33.0%
   of variance in all-individuals DLPFC, compared with 29.1% in AA DLPFC, 26.8%
   in all-individuals hippocampus, and 13.4% in all-individuals caudate. The
   legacy DLPFC 99th-percentile residual-SD cutoff was therefore substantially
   inflated; for example, it was 0.1068 versus 0.0647 in v2 on chromosome 1 and
   0.1018 versus 0.0654 on chromosome 22. Because this cutoff selects the seed
   CpGs passed to `regionFinder3()`, the distortion propagates directly into
   VMR membership and boundaries.
4. **Correct alignment restores the expected cross-arm relationship.** The
   autosomal AA-versus-all-individuals base-pair Jaccard in DLPFC increased from
   0.312 between the two corrupted legacy catalogs to 0.684 between the two v2
   catalogs. The corresponding within-version AA-versus-all comparisons were
   0.691 in legacy and 0.665 in v2 for caudate, and 0.679 in legacy and 0.697 in
   v2 for hippocampus. DLPFC is therefore the only region with anomalous
   cross-arm disagreement before the repair and ordinary cross-arm agreement
   after it.

The turnover calculation itself was independently reproduced using the
arm-specific legacy path, autosomes only, and the same interval-union
definition as `04_turnover.R`. The audit also rules out the retired sample
blacklist, restored donors, sex-chromosome policy, HDF5SE source, and a single
corrupted chromosome as primary explanations. The full frozen legacy catalog
was not regenerated under corrected ordering; the six-chromosome-block
realignment is a forensic sensitivity, not an accepted replacement run.

**Acceptance interpretation.** The 72.2% novel fraction should remain reported
as the observed legacy-to-v2 turnover, but it must not be interpreted as
DLPFC-specific biology or as a general property of the all-individuals arm. It
is explained by branch-specific legacy donor permutation interacting with
strong DLPFC methylation-PC structure. This explanation clears turnover as
evidence of a new v2 implementation defect.

## What has been verified

A chr22-only caudate AA smoke run exercised every step end-to-end:

- **489,722 CpGs × 153 donors** survive C→T masking and coverage QC
  (154 in the phenotype table → 153 after intersecting BSobj and `.psam`).
- **260 VMRs** called; 260 phenotype files written; membership table covers
  4,140 CpGs; reconciliation clean (1 completed, 21 QC-failed for the
  chromosomes this smoke run did not build, 0 unaccounted, 0 failed).
- Run sealed: 289 outputs checksummed, directory read-only.

**V1 shuffle invariance, on real data.** Permuting the methylation rows and the
PC rows independently and re-running the residualization changed the per-CpG
residual SD by at most `5.6e-17` — floating-point noise. The alignment holds.

**What V1 costs, on real data.** Permuting only the design rows — the legacy
failure mode — on all 489,722 chr22 CpGs:

| | Correctly aligned | Design misaligned |
|---|---|---|
| 99th-pct residual SD cutoff | 0.074733 | 0.083297 (**+11.5%**) |
| Seed CpGs above cutoff | 4,898 | — |
| Jaccard of seed-CpG sets | — | 0.575 |
| Seed CpGs that change | — | **26.9%** |

Direction and magnitude agree with the audit's +9.2% and Jaccard 0.82 measured
on the legacy files. (The audit's figures come from the *actual* legacy
permutation, which partially overlaps the correct order; a fully random
permutation is a worse case, hence the larger numbers here.)

**V4 corroborated.** The legacy caudate catalog carries **431** sex-chromosome
VMRs, matching the audit's reported 3× excess exactly. v2 holds them out of the
primary catalog.

**Turnover, chr22 only.** 260 v2 vs 244 legacy VMRs; 180 overlap, 80 novel,
68 legacy lost; Jaccard 0.51; 31% of v2 VMRs novel. Substantial turnover is the
expected consequence of the V1 repair, not a surprise.

### Verification not yet done

- **Full-scale runs.** AA DLPFC is done. Caudate and hippocampus AA, and all
  three `all_individuals` cells, have not been run. Hippocampus shares DLPFC's
  HDF5SE BSobj source (V14) and has never been run at any scale.
- **`design_n`.** AA DLPFC observed **118**, matching the reconstruction — the
  PI still has to lock it in `config/cohorts.yml`. Caudate (**153**) and
  hippocampus (**117**) remain reconstructed from the inputs, not observed.
- **Checksum stability across two invocations** (acceptance criterion 3).
- **Array-coverage panel.** `inputs/supportfiles/_m/array_cpg_manifest_hg38.bed.gz`
  does not exist, so `04_turnover.R` skips the off-array comparison and records
  that it did. AGENTS.md §11 Figure 1 needs it before the figure freeze.

## Smoke tests

Gitignored, hand-run before submitting an array:

```bash
Rscript -e 'source("00_shared/load.R"); testthat::test_dir("01_vmr_catalog/tests")'
```

## `_m/` contents

`_m/runs/` is gitignored — it holds large generated matrices and stays on Quest.
What is tracked is this README, the acceptance table above, and the run IDs it
cites. A run is regenerated by checking out the recorded git commit, restoring
the recorded config checksums, and rerunning `submit_vmr_catalog_workflow.sh`.
Completed runs are made read-only; never edit a result in `_m/`.
