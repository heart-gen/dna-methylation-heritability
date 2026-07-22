# meQTL validation (strengthening)

CpG-level cis-meQTL validation of VMR local genetic predictability lives under
`meqtl-validation/`, using the repo `_h` (scripts) / `_m` (outputs) convention.

| Module | Phase | Status |
|---|---|---|
| [`data_audit/`](data_audit/) | 0 | Complete (earlier pilot audit) |
| [`cpg_meqtl_mapping/`](cpg_meqtl_mapping/) | 1 | Pipeline ready; run TensorQTL on cluster |
| [`vmr_meqtl_burden/`](vmr_meqtl_burden/) | 2 | Scripts ready; needs Phase 1 outputs |
| [`external_meqtl_validation/`](external_meqtl_validation/) | 3 | Catalog + workspace; needs downloads |
| [`cross_region_sharing/`](cross_region_sharing/) | 4 | Script ready |
| [`donor_group_comparison/`](donor_group_comparison/) | 4 | Script ready |
| [`transcription_splicing_integration/`](transcription_splicing_integration/) | 5 | Script ready |
| [`repeat_mappability_sensitivity/`](repeat_mappability_sensitivity/) | 6 | Robustness table scaffold |
| [`schizophrenia_risk_application/`](schizophrenia_risk_application/) | 7 | **Deferred** |

Shared helpers: [`_lib/`](_lib/).
Config: [`config/`](../config/).
Data dictionary: [`inputs/data_dictionary/`](../inputs/data_dictionary/).
