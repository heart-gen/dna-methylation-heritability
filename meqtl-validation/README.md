# meQTL validation

CpG-level cis-meQTL validation of VMR local genetic predictability.
Uses the repo `_h` (scripts) / `_m` (outputs) convention. Submit from each module `_m/`.

**Analysis summary:** [`ANALYSIS_SUMMARY_meqtl_validation.md`](ANALYSIS_SUMMARY_meqtl_validation.md). Phases 1–7 and the cross-region, donor-group, downsampling, and cell-type sensitivity analyses are complete.

| Module | Phase | Status | Submission |
|---|---|---|---|
| [`00_data_audit/`](00_data_audit/) | 0 | Done | `sbatch ../_h/step_1.sh` |
| [`01_cpg_meqtl_mapping/`](01_cpg_meqtl_mapping/) | 1 | **Locked M3a** | primary promoted; see `PHASE1_LOCK_DECISION.md` |
| [`02_vmr_meqtl_burden/`](02_vmr_meqtl_burden/) | 2 | **Pass on M3a** | `step_1` aggregate → `step_2` models |
| [`03_external_meqtl_validation/`](03_external_meqtl_validation/) | 3 | **Pass** (Jaffe/Schulz) | `step_1` init → `step_2` external test |
| [`04_cross_region_sharing/`](04_cross_region_sharing/) | 4a | **Pass on M3a** + Exp3 Pass (lead retention + G×region) | `step_1` → `step_3_downsample` → `step_4_gxregion` → `step_2` |
| [`05_donor_group_comparison/`](05_donor_group_comparison/) | 4b | **Pass** (EA meQTL+burden; Exp2 depth complete) | `POPULATION=EA` Phase1/2 then `step_1` → `step_2` |
| [`06_transcription_splicing_integration/`](06_transcription_splicing_integration/) | 5 | **Pass on M3a** | `step_1` expression → `step_2` PSI → `step_3` |
| [`07_repeat_mappability_sensitivity/`](07_repeat_mappability_sensitivity/) | 6 | **Pass on M3a** | `step_2_tech_joins` then `step_1`; plan: `CELLTYPE_LINE_L1_PLAN.md` |
| [`08_schizophrenia_risk_application/`](08_schizophrenia_risk_application/) | 7 | **Complete** (proof-of-application) | Analyses 1–7 + Level3 + locus panels (`_m/locus_panels/figures/`) |
| [`09_libd_eqtl_mapping/`](09_libd_eqtl_mapping/) | support | **Debug** | BrainSeq/LIBD AA gene eQTL; see `EQTL_DEBUG_TODO.md` |
| [`10_downsampling_caudate/`](10_downsampling_caudate/) | Exp3 formal | **Pass (complete)** | official TensorQTL perm-FDR remap (30× N=111); claim snapshot written |
| [`11_celltype_compartment_sensitivity/`](11_celltype_compartment_sensitivity/) | 6 cell-type | **Pass** | MuSiC/DNAm cellPC_r2 adj; keep LINE/L1 main-figure row |

Shared helpers: [`_lib/`](_lib/).
Config: [`config/`](../config/).
Data dictionary: [`inputs/data_dictionary/`](../inputs/data_dictionary/).

## Example (Phase 1)

```bash
cd meqtl-validation/01_cpg_meqtl_mapping/_m
mkdir -p logs
sbatch ../_h/step_1.sh   # CPU: preflight + covariates
sbatch ../_h/step_2.sh   # CPU: CpG BEDs
sbatch ../_h/step_3.sh   # CPU: genotypes
sbatch ../_h/step_4.sh   # GPU: TensorQTL
sbatch ../_h/step_5.sh   # CPU: QC
```
