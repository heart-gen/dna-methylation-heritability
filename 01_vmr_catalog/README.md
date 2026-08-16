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

Launchers: `step_1.sh` (autosomes), `step_1x.sh` (X/Y → `excluded/`),
`step_2.sh`, `step_3.sh`, `step_4.sh`, and
`submit_vmr_catalog_workflow.sh` to chain them.

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
| — | Hippocampus read the **DLPFC** sample blacklist | Per-region resolution, no fallback | `00_shared/config.R` |

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
| _(none accepted)_ | | | | | | |

### Verification not yet done

- **Full-scale runs.** Only chr22 caudate AA has been exercised end-to-end.
- **`design_n`.** chr22 caudate AA gives n=153 (phenotype ceiling 154). DLPFC and
  hippocampus AA must be checked against the audit's expected 111 and 116 rather
  than the legacy 96 and 101 — that recovery is the observable signature of the
  V2 fix.
- **Array-coverage panel.** `inputs/supportfiles/_m/array_cpg_manifest_hg38.bed.gz`
  does not exist yet, so `04_turnover.R` skips the off-array comparison and says
  so. AGENTS.md §11 Figure 1 needs it before the figure freeze.

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
