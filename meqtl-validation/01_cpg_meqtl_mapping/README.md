# Phase 1: Internal CpG-level cis-meQTL mapping

## Question

Do CpGs within predefined VMRs show conventional cis-meQTL associations in each brain region?

## Design (prespecified)

- Unit: CpG methylation (non-residualized `cpg_meth.phen` from Alexis repo)
- Discovery: Black American donors
- Cis window: ±500 kb
- FDR: separate family per brain region
- Engine: TensorQTL (`/projects/p32505/opt/envs/genomics` conda env)
- Covariates: **locked M3a** — age, sex, diagnosis, snpPC1–5, methPC1–5 (`config/covariates.yml`; see `PHASE1_LOCK_DECISION.md`)
- Genotypes: filtered population-specific pfile (`meqtl_AA` or `meqtl_EA`) with `FID=IID=BrNum`
- EA outputs are isolated under `prepared/EA/` and `tensorqtl/EA/` (AA primary M3a untouched)
- Submit with `POPULATION=EA` for stratified mapping, e.g. `sbatch --export=ALL,POPULATION=EA ../_h/step_1.sh`
- EA-M3a caudate sensitivity (EA-estimated methPC1–5): `step_ea_m3a_caudate.sh` → `step_ea_m3a_tensorqtl.sh` → `02_vmr_meqtl_burden/_h/step_ea_m3a_burden.sh` (writes `tensorqtl/EA/M3a/` and `02_.../_m/EA_M3a/`; does not overwrite EA M0)

## Run (Quest)

Submit from `_m/`. Steps are intentionally separate for CPU vs GPU and easier troubleshooting.

```bash
cd meqtl-validation/01_cpg_meqtl_mapping/_m
mkdir -p logs

sbatch ../_h/step_1.sh   # CPU: preflight + covariates
sbatch ../_h/step_2.sh   # CPU: CpG phenotype BEDs + merge
sbatch ../_h/step_3.sh   # CPU: genotype filtering (plink2)
sbatch ../_h/step_4.sh   # GPU (gengpu): TensorQTL cis mapping
sbatch ../_h/step_5.sh   # CPU: QC summaries
sbatch ../_h/step_6.sh   # CPU: real CpG coverage from stats.rda
sbatch ../_h/step_7.sh  # CPU: QQ + λ_GC by qval stratum

# Phase 1 covariate / latent sensitivity (does not overwrite primary tensorqtl/):
sbatch --export=ALL,REGION=caudate --array=0 ../_h/step_8_latent.sh
sbatch --dependency=afterok:<step8_jobid> --export=ALL,REGION=caudate ../_h/step_9_covar_models.sh
sbatch --dependency=afterok:<step9_jobid> --export=ALL,REGION=caudate ../_h/step_10_compare_models.sh
```

After coverage array finishes:

```bash
for r in caudate dlpfc hippocampus; do
  python3 ../_h/07_aggregate_vmr_coverage.py --region "$r"
done
```

Single region example:

```bash
sbatch --export=ALL,REGION=caudate --array=0 ../_h/step_1.sh
```

CPU fallback for TensorQTL (if GPU queue is unavailable):

```bash
sbatch --partition=short --export=ALL,DEVICE=cpu ../_h/step_4.sh
```

## Scripts

| Script | Purpose | Resources |
|---|---|---|
| `00_preflight.py` | Sample inclusion lists | CPU (`step_1`) |
| `01_prepare_covariates.py` | Covariate matrix | CPU (`step_1`) |
| `02a_prepare_cpg_bed.py` | Per-chromosome VMR CpG BED | CPU (`step_2`) |
| `02b_merge_cpg_beds.py` | Merge autosomal BEDs + bgzip/tabix | CPU (`step_2`) |
| `03_prepare_genotypes.sh` | plink2 filter | CPU (`step_3`) |
| `04_tensorqtl_map.py` | cis TensorQTL | GPU (`step_4`) |
| `05_qc_summarize.py` | λ_GC, histograms, lead SNP table | CPU (`step_5`) |
| `06_extract_cpg_coverage.R` | Real depth from `stats.rda` BSobj | CPU (`step_6`) |
| `07_aggregate_vmr_coverage.py` | VMR-level mean coverage | CPU (post step_6) |
| `08_calibration_plots.py` | QQ plots + λ_GC by qval stratum | CPU (`step_7`) |
| `09_estimate_latent_factors.py` | PCA residual latents (PEER-style) | CPU (`step_8_latent`) |
| `10_prepare_covariate_models.py` | Build M0–M5 covariate matrices | CPU (`step_8_latent`) |
| `11_compare_covariate_models.py` | λ_NS / n_sig / external overlap | CPU (`step_10_compare`) |
| `12_lock_decision.py` | Write `PHASE1_LOCK_DECISION.md` | CPU (`step_10_compare`) |
| `step_9_covar_models.sh` | Sensitivity TensorQTL per model | GPU |

## Outputs

- `{region}/_m/preflight/` — inclusion lists
- `{region}/_m/prepared/` — BED phenotypes (bgzip+tabix), covariates (samples × covariates), CpG–VMR map
- `{region}/_m/genotypes/` — filtered pfiles (`FID=IID=BrNum`) + `genotype_qc_summary.tsv`
- `{region}/_m/tensorqtl/` — cis QTL results
- `{region}/_m/tensorqtl/qc/` — QC summaries
- `{region}/_m/tensorqtl/qc/calibration/` — QQ plots + λ_GC by qval stratum

## Notes

- Conda env: `/projects/p32505/opt/envs/genomics` (tensorqtl + plink2 + bgzip/tabix)
- Source `TOPMed_LIBD.AA.psam` lacks a plink2 header; step_3 stages a headered copy under `_m/genotype_source/`
- Imputation INFO/R2 is absent from the staged pvar because the distributed panels
  were already filtered upstream. Step 3 verifies the resource-construction logs:
  AA was built from `LIBD_postmortem.AA.info_0.8...`; AA+EA was built from
  `Merged_Info0.8/...`. An optional `IMPUTATION_R2_IDS` list can impose an
  additional filter, but it is not required. The unverified escape hatch remains
  debugging-only and is recorded in the QC summary.
- TensorQTL loads genotypes in chromosome chunks; variant contigs are chr-prefixed to match phenotype BEDs
- Submit scripts from `_m/` so relative `../_h/` paths resolve (SLURM spool copies break `$0`-based paths)
- Phen matrices are already depth-filtered in Alexis prep (`cov>=5` in ≥80% samples); `step_6` attaches per-CpG `mean_coverage` for Phase 2 matching
- Primary CpG significance uses Storey **qval** (per-region FDR); `step_7` reports λ_GC overall and within qval strata (does not change calls)
- **Phase 1 locked to M3a** for AA primary — see `PHASE1_LOCK_DECISION.md`. EA stratified uses M0 covariates under `prepared/EA/` / `tensorqtl/EA/`.
