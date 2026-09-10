# 01b_estimation_cells — donor-group estimation cells on a pooled catalog

This module materializes an **estimation cell**: a donor subset of an accepted,
sealed Module 01 run, packaged so that Modules 02 and 03 can be fit inside it
without changing anything else about the analysis.

## Why it exists

`config/cohorts.yml` defines two discovery **arms**, `AA` and
`all_individuals`, and until 2026-09-10 every stage ran inside one arm end to
end: VMR discovery, cis-genotype extraction, Module 02 and Module 03 all read
the same arm's pgen and the same arm's donors.

That is the wrong split for the donor-group question. The biologically correct
design — the one v1 used, and the one the PI locked on 2026-09-06 (AGENTS.md
§7.7) — separates the two things the arm was conflating:

- **discovery happens once, pooled.** VMRs are called on all individuals, so
  both donor groups are evaluated on one fixed locus set.
- **SNP modeling happens within group.** Modules 02 and 03 are the
  SNP-consuming stages; they run separately in Black American and in
  non-Hispanic white American donors, so MAF, LD and effect estimation are
  within-group.

AGENTS.md §7.7 also forbids the alternative: *"Do not contrast `AA` against
`all_individuals`: those are nested, and a set-versus-superset comparison is not
a donor-group contrast."* The contrast is `all_individuals.AA` against
`all_individuals.EA`, and only that pair.

This is the v1 design rebuilt on the v2 run contract. v1 wrote `samples-AA.txt`
/ `samples-EA.txt` at discovery time
(`vmr-analysis/all_individuals/{region}/_h/01.get_cpg_stats.R:112-133`), ran
plink2 twice per VMR window with `--keep`
(`.../_h/step_5.sh:95-122`), modeled each group separately
(`local-snp-prediction/all_individuals/{region}/_h/01.elastic-net.R`,
`population`-parameterized) and recombined by an inner join on locus
coordinates with `_AA` / `_EA` suffixes
(`.../venn_diagram/_h/02.compare_cohort.R:53-62`).

## Why it is a separate module

Runs are immutable (AGENTS.md §5.2) and the pooled catalogs
`vmrcat-all_individuals-{region}-20260816` are accepted and **sealed**.
Per-group genotype extraction cannot be added to them in place. So this module
reads a sealed Module 01 run and writes a *new* immutable run that satisfies the
same directory contract `00_shared/locus_io.R` already consumes — `vmr/`,
`covs/`, `plink_format/` — which is what keeps the change out of Modules 02
and 03 almost entirely.

**This module never re-derives VMRs, methylation residuals or phenotypes.** It
re-partitions donors over a fixed locus set. Stage 04 refuses to seal a run
whose `vmr.bed` or `vmr_catalog.tsv` differs from its source by so much as a
checksum, because a divergence there would put the two donor groups on
different loci and void the design.

## The cell token

`{catalog_cohort}.{estimation_group}` — `all_individuals.AA`,
`all_individuals.EA`. A **bare arm token** (`AA`, `all_individuals`) parses as a
cell whose estimation group equals its cohort, which is what every
pre-2026-09-10 run is. That is what lets the token go wherever the cohort token
went — run IDs, README acceptance tables, `00_shared/gates.R`, every downstream
manifest field — with no schema change and no accepted run changing meaning.

`EA` is deliberately **not** an arm. Making it one would let Module 01 discover
VMRs in EA donors only, which is precisely the design being replaced.

Resolve a token with `parse_cell()` / `cell_def()` in `00_shared/config.R`;
`cohort_def()` still takes a discovery arm and refuses a cell rather than
silently returning `NULL`.

## Stages

| Stage | Script | Contract |
|---:|---|---|
| 00 | `_h/00_new_run.R` | Verify the upstream Module 01 run is accepted, non-smoke, and on the cell's **discovery arm**; derive the group donor list from the pooled list; prove the cells on that catalog are disjoint and cover every pooled donor; open the run and copy `vmr.bed` + `vmr_catalog.tsv` verbatim. |
| 01 | `_h/step_1_group_pca.sh` → `_h/01_write_genotype_pcs.R` | LD-pruned plink2 PCA on the pooled pgen restricted to this cell's donors; normalize to `covs/genotype_pcs.tsv` (`FID IID snpPC1..k`). |
| 02 | `_h/step_2_extract_cis.sh` (array) | Extract each VMR's cis window for this cell's donors, writing `plink_format/chr_N/TOPMed_LIBD-{group}.{start}_{end}.bed`. Carries the V9/V10/V11/V12 repairs from `01_vmr_catalog/_h/step_4.sh` unchanged. |
| 03 | `_h/03_subset_covariates.R` | Restrict the catalog's covar/qcovar and per-VMR phenotypes to the cell's donors, in the cell's donor order. Copies values; recomputes nothing. |
| 04 | `_h/04_close_run.R` | Reconcile every VMR as extracted or `no_cis_variants`, cross-check the extraction log against the filesystem, verify the locus set and the PC coverage, checksum and seal. |

Every `step_*.sh` reads `RUN_ID` from the environment and takes no positional
arguments, so any single stage can be resubmitted by hand against an existing
run without editing the workflow script.

## Within-group genotype PCs

Locus-level covariates in `00_shared/locus_io.R` have always been
`age + sex + diagnosis`, with population structure removed upstream by
residualizing methylation on **pooled** `snpPC1-3` during discovery. That is
correct for pooled discovery, but it does not control structure *within* an
AA-only or EA-only estimation set: restricted to one group, the pooled leading
PCs are close to constant.

PI decision 2026-09-10 (`config/covariates.yml`,
`estimation_cells.genotype_pcs`): estimation cells recompute genotype PCs in
their own donors and carry them as locus-level covariates. `load_observed_locus()`
appends them when `covs/genotype_pcs.tsv` is present. **When that file is absent
the covariate matrix is byte-identical to the pre-2026-09-10 behaviour**, which
is what lets the six sealed accepted Module 02 runs reproduce exactly.

## Recombination

Downstream loaders stay strictly single-cell —
`00_shared/gates.R::load_local_genetic_control()` errors on mixed cohorts by
design and that guard is not loosened. Recombination is two explicit stages:

- `02_local_genetic_variance/_h/14_combine_donor_group_cells.R` →
  `_m/combined/local-genetic-control-donor-group-{region}.tsv`, per-locus with
  `_AA` / `_EA` suffixes. **No pooled rank, no cross-cell score difference.**
  The score is a within-cell midrank percentile and
  `config/local_genetic_control.yml` locks `rank_scope: cohort_by_region`; a
  pooled rank would be exactly the quantity AGENTS.md §7.6 prohibits. Ordering
  agreement (Spearman, quartile concordance) is what the table supports.
- `03_local_snp_prediction/_h/07_stack_donor_group_predictions.R` →
  `_m/combined/oof-predictions-per-donor-donor-group-{region}.tsv`, the two
  cells' per-donor OOF predictions concatenated with a `donor_group` column.
  Stacking predictions is well defined (same phenotype, same loci, disjoint
  donors); stacking *accuracy* is not, so `r2_pred_oof` stays per cell.

## Interpretation constraints

The two cells differ in sample size (EA is roughly half of AA in each region),
MAF spectrum, LD structure and SNP availability. AGENTS.md §7.7: *"Do not
attribute differences to ancestry-specific biology without eliminating sample
size, MAF, LD, SNP availability, assay, covariate, and brain-region
explanations."* AGENTS.md §7.6 separately forbids raw score-level comparison
across cells, so compare **ordering**, never levels.

## Running it

```bash
cd 01b_estimation_cells/_m && mkdir -p logs
COHORT=all_individuals GROUP=EA REGION=dlpfc \
  VMR_RUN_ID=vmrcat-all_individuals-dlpfc-20260816 \
  ../_h/submit_estimation_cell.sh
```

`DRY_RUN=TRUE` prints the plan without submitting; `ALLOW_UNLOCKED=TRUE` stamps
the run `smoke_run=TRUE` so it can never be accepted. Stage 00 runs on the
submit host — it mints the run ID every later stage needs, and the VMR count it
copies is what sizes the stage-02 array.

A completed SLURM chain is **not** acceptance. Stage 04 must seal the run and a
human must record it below (AGENTS.md §6).

## Expected donor counts

Reconstructed from the locked counts in `config/cohorts.yml`; the two cells
partition the pooled analysis set exactly.

| Region | pooled | `all_individuals.AA` | `all_individuals.EA` |
|---|---:|---:|---:|
| caudate | 282 | 153 | 129 |
| dlpfc | 173 | 118 | 55 |
| hippocampus | 177 | 117 | 60 |

`design_n` is `null` for every cell until the first accepted run reports
observed counts and the PI signs off. While null,
`00_shared/identity.R::assert_expected_n()` records and warns; once set, a
mismatch is a hard stop.

## Acceptance gate

A run is acceptable when **all five** hold:

1. every VMR in the catalog is accounted for as extracted or `no_cis_variants`,
   with zero unaccounted and zero computational failures;
2. `vmr.bed` and `vmr_catalog.tsv` checksum-match the source Module 01 run;
3. the cell's donors are a strict subset of the pooled donors, and the cells on
   that catalog are disjoint and cover every pooled donor;
4. the observed donor count matches the reconstruction in `config/cohorts.yml`;
5. `covs/genotype_pcs.tsv` covers exactly the cell's donors, in order, with no
   missing or zero-variance component.

### Gate results

Produced by `_h/05_report_acceptance.R`, which re-checks all five against the
sealed artifacts rather than trusting Stage 04's in-flight assertions.

| run_id | cell | region | donors | VMRs | PCs | vmr_set_id | gate | notes |
|---|---|---|---|---:|---:|---|---|---|
| `estcell-all_individuals.AA-caudate-20260910` | all_individuals.AA | caudate | **153** | 11,463 | 3 | `vmrset-all_individuals-caudate-cb5519d7d2ad` | **all five pass** | 153 no cis variant; loci from `vmrcat-all_individuals-caudate-20260816` |
| `estcell-all_individuals.AA-dlpfc-20260910` | all_individuals.AA | dlpfc | **118** | 9,374 | 3 | `vmrset-all_individuals-dlpfc-e88f46904afb` | **all five pass** | 181 no cis variant; loci from `vmrcat-all_individuals-dlpfc-20260816` |
| `estcell-all_individuals.AA-hippocampus-20260910` | all_individuals.AA | hippocampus | **117** | 9,365 | 3 | `vmrset-all_individuals-hippocampus-809f8de0db2d` | **all five pass** | 178 no cis variant; loci from `vmrcat-all_individuals-hippocampus-20260816` |
| `estcell-all_individuals.EA-caudate-20260910` | all_individuals.EA | caudate | **129** | 11,463 | 3 | `vmrset-all_individuals-caudate-cb5519d7d2ad` | **all five pass** | 153 no cis variant; loci from `vmrcat-all_individuals-caudate-20260816` |
| `estcell-all_individuals.EA-dlpfc-20260910` | all_individuals.EA | dlpfc | **55** | 9,374 | 3 | `vmrset-all_individuals-dlpfc-e88f46904afb` | **all five pass** | 181 no cis variant; loci from `vmrcat-all_individuals-dlpfc-20260816` |
| `estcell-all_individuals.EA-hippocampus-20260910` | all_individuals.EA | hippocampus | **60** | 9,365 | 3 | `vmrset-all_individuals-hippocampus-809f8de0db2d` | **all five pass** | 178 no cis variant; loci from `vmrcat-all_individuals-hippocampus-20260816` |

## Accepted runs

Machine-readable, in the schema `00_shared/gates.R::read_accepted_runs()`
parses. Module 02 additionally requires the row to say "all five pass".

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| estcell-all_individuals.AA-caudate-20260910 | all_individuals.AA | caudate | vmrset-all_individuals-caudate-cb5519d7d2ad | 2026-09-10 |  | ALL_FIVE_CRITERIA_PASS | 153 donors, 11463 VMRs (153 no cis variant), 3 within-group PCs; loci from vmrcat-all_individuals-caudate-20260816 |
| estcell-all_individuals.AA-dlpfc-20260910 | all_individuals.AA | dlpfc | vmrset-all_individuals-dlpfc-e88f46904afb | 2026-09-10 |  | ALL_FIVE_CRITERIA_PASS | 118 donors, 9374 VMRs (181 no cis variant), 3 within-group PCs; loci from vmrcat-all_individuals-dlpfc-20260816 |
| estcell-all_individuals.AA-hippocampus-20260910 | all_individuals.AA | hippocampus | vmrset-all_individuals-hippocampus-809f8de0db2d | 2026-09-10 |  | ALL_FIVE_CRITERIA_PASS | 117 donors, 9365 VMRs (178 no cis variant), 3 within-group PCs; loci from vmrcat-all_individuals-hippocampus-20260816 |
| estcell-all_individuals.EA-caudate-20260910 | all_individuals.EA | caudate | vmrset-all_individuals-caudate-cb5519d7d2ad | 2026-09-10 |  | ALL_FIVE_CRITERIA_PASS | 129 donors, 11463 VMRs (153 no cis variant), 3 within-group PCs; loci from vmrcat-all_individuals-caudate-20260816 |
| estcell-all_individuals.EA-dlpfc-20260910 | all_individuals.EA | dlpfc | vmrset-all_individuals-dlpfc-e88f46904afb | 2026-09-10 |  | ALL_FIVE_CRITERIA_PASS | 55 donors, 9374 VMRs (181 no cis variant), 3 within-group PCs; loci from vmrcat-all_individuals-dlpfc-20260816 |
| estcell-all_individuals.EA-hippocampus-20260910 | all_individuals.EA | hippocampus | vmrset-all_individuals-hippocampus-809f8de0db2d | 2026-09-10 |  | ALL_FIVE_CRITERIA_PASS | 60 donors, 9365 VMRs (178 no cis variant), 3 within-group PCs; loci from vmrcat-all_individuals-hippocampus-20260816 |

All six rows were recorded on 2026-09-10 from the evidence
`_h/05_report_acceptance.R` produced, which re-checks the five criteria against
the SEALED artifacts rather than trusting Stage 04's in-flight assertions.

`accepted_by` is blank pending PI review. `read_accepted_runs()` does not
require it, so Modules 02 and 03 will consume these runs; the column records who
looked at the run, which is a human act (AGENTS.md 6) and not something the
evidence script can supply.

**Passing this gate does not mean a cell can produce a defensible Module 02
score.** These five criteria establish that the donor re-partition is faithful
and complete -- nothing more. In particular the EA cells' donor counts (129, 55,
60) are absent from `allowed_n` in
`02_local_genetic_variance/config/joint-pve-characterized-support.tsv`, so every
EA locus would be labelled outside the frozen model's characterized domain and
Stage 05 would fail the run. An EA cell needs its own observed-regime grid
(`02_local_genetic_variance/_h/submit_observed_regime_cell.sh`) before its n may
enter that file. The AA cells' counts (153, 118, 117) are already present.
