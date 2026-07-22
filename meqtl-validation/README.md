# meQTL validation

CpG-level cis-meQTL validation of VMR local genetic predictability.
Uses the repo `_h` (scripts) / `_m` (outputs) convention. Submit from each module `_m/`.

| Module | Phase | Submission |
|---|---|---|
| [`00_data_audit/`](00_data_audit/) | 0 | `sbatch ../_h/step_1.sh` |
| [`01_cpg_meqtl_mapping/`](01_cpg_meqtl_mapping/) | 1 | `step_1` prep → `step_2` TensorQTL → `step_3` QC |
| [`02_vmr_meqtl_burden/`](02_vmr_meqtl_burden/) | 2 | `step_1` aggregate → `step_2` models |
| [`03_external_meqtl_validation/`](03_external_meqtl_validation/) | 3 | `step_1` init → `step_2` external test |
| [`04_cross_region_sharing/`](04_cross_region_sharing/) | 4 | `sbatch ../_h/step_1.sh` |
| [`05_donor_group_comparison/`](05_donor_group_comparison/) | 4 | `sbatch ../_h/step_1.sh` |
| [`06_transcription_splicing_integration/`](06_transcription_splicing_integration/) | 5 | `step_1` expression → `step_2` PSI |
| [`07_repeat_mappability_sensitivity/`](07_repeat_mappability_sensitivity/) | 6 | `sbatch ../_h/step_1.sh` |
| [`08_schizophrenia_risk_application/`](08_schizophrenia_risk_application/) | 7 | Deferred stub `step_1.sh` |

Shared helpers: [`_lib/`](_lib/).
Config: [`config/`](../config/).
Data dictionary: [`inputs/data_dictionary/`](../inputs/data_dictionary/).

## Example (Phase 1)

```bash
cd meqtl-validation/01_cpg_meqtl_mapping/_m
mkdir -p logs
sbatch ../_h/step_1.sh          # prep all regions (array)
sbatch ../_h/step_2.sh          # TensorQTL (after prep)
sbatch ../_h/step_3.sh          # QC summaries
```
