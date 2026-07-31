# meQTL validation

CpG-level cis-meQTL validation of VMR local genetic predictability.
Uses the repo `_h` (scripts) / `_m` (outputs) convention. Submit from each module `_m/`.

**Overall readiness:** [`../AJHG_SUBMISSION_READINESS.md`](../AJHG_SUBMISSION_READINESS.md) (updated 2026-07-30). Decision: **CONDITIONAL-GO**.

| Module | Phase | Status | Submission |
|---|---|---|---|
| [`00_data_audit/`](00_data_audit/) | 0 | Done | `sbatch ../_h/step_1.sh` |
| [`01_cpg_meqtl_mapping/`](01_cpg_meqtl_mapping/) | 1 | Usable; lock pending | `step_1–5` + `step_8–10` covariate sensitivity |
| [`02_vmr_meqtl_burden/`](02_vmr_meqtl_burden/) | 2 | **Pass** | `step_1` aggregate → `step_2` models |
| [`03_external_meqtl_validation/`](03_external_meqtl_validation/) | 3 | **Pass** (Jaffe/Schulz) | `step_1` init → `step_2` external test |
| [`04_cross_region_sharing/`](04_cross_region_sharing/) | 4a | **Pass** (downsample remap pending) | `step_1` → `step_2` summary |
| [`05_donor_group_comparison/`](05_donor_group_comparison/) | 4b | Partial (EA paths ready; meQTL not run) | `sbatch ../_h/step_1.sh` |
| [`06_transcription_splicing_integration/`](06_transcription_splicing_integration/) | 5 | **Pass** | `step_1` expression → `step_2` PSI → `step_3` |
| [`07_repeat_mappability_sensitivity/`](07_repeat_mappability_sensitivity/) | 6 | **Pass** | `sbatch ../_h/step_1.sh` |
| [`08_schizophrenia_risk_application/`](08_schizophrenia_risk_application/) | 7 | Deferred | stub `step_1.sh` |

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
