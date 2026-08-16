# 01_vmr_catalog — corrected biological units

Defines the variably methylated regions (VMRs) every downstream v2 module is
built on, for two cohort arms across three brain regions, from one parameterized
codepath.

**Status: implemented, no accepted production run yet.** See *Acceptance gate*
below. Nothing downstream may consume this module until a run ID is recorded
here with a passing gate (AGENTS.md §6).

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

| Run ID | Arm | Region | n | vmr_set_id | Gate | Notes |
|---|---|---|---|---|---|---|
| `vmrcat-AA-caudate-20260815-a` | AA | caudate | 153 | `vmrset-AA-caudate-3968cebd6d9c` | **smoke, not production** | chr22 only, `--allow-unlocked` |

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

- **Full-scale runs.** Only chr22 caudate AA. No other chromosome, region, or
  arm has been run.
- **`design_n`.** DLPFC and hippocampus AA must be checked against the expected
  **118** and **117** rather than the legacy 96 and 101 — recovering those ~20
  donors per region is the observable signature of the V2 fix plus the retired
  blacklist, and it has **not** been demonstrated by a run yet. The counts are
  reconstructed from the inputs (phenotype → BSobj → psam), not observed.
- **Checksum stability across two invocations** (acceptance criterion 3).
- **`step_4.sh`** genotype extraction has never been executed.
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
