# meQTL validation

CpG-level cis-meQTL validation of VMR local genetic predictability.
Uses the repo `_h` (scripts) / `_m` (outputs) convention. Submit from each module `_m/`.

> **Repair-v2 status (updated 2026-08-15):** the repair-v2 production rerun
> **has been executed** (2026-08-08 → 2026-08-09) across all modules in the order
> given in [`REPAIR_V2.md`](REPAIR_V2.md). All `_m/` tables below are schema-v2
> outputs. The earlier "outputs stale — must be regenerated" banner is obsolete.
>
> **Two claims failed on regeneration.** Phase 3 external validation (module `03`) was
> an input-labeling bug and is **fixed and passing** as of 2026-08-15. The
> caudate-not-sample-size claim (module `10`) genuinely fails and should be retired;
> modules `04`, `05`, and `08` still report the pre-repair caudate verdict. See
> [`../AJHG_SUBMISSION_READINESS.md`](../AJHG_SUBMISSION_READINESS.md) §4 and §7 before
> using those tables.

**Historical analysis summary:** [`ANALYSIS_SUMMARY_meqtl_validation.md`](ANALYSIS_SUMMARY_meqtl_validation.md)
(pre-repair; retained for provenance only — do not quote its numbers).

## Module status

| Module | Phase | Rerun | Result | Notes |
|---|---|---|---|---|
| [`00_data_audit/`](00_data_audit/) | 0 | n/a | Done | |
| [`01_cpg_meqtl_mapping/`](01_cpg_meqtl_mapping/) | 1 | 2026-08-08 | **Complete** | M3a; AA + EA. n_sig: caudate 80,281 / DLPFC 50,288 / hip 54,898 |
| [`02_vmr_meqtl_burden/`](02_vmr_meqtl_burden/) | 2 | 2026-08-08 | **Pass** | schema v2; AA adj-prespec coef 3.12 / 1.87 / 2.15; EA 3.43 / 1.12 / 1.38 |
| [`03_external_meqtl_validation/`](03_external_meqtl_validation/) | 3 | 2026-08-15 | **Pass** | `Criterion 5: 2/2`. Jaffe/DLPFC coef +0.68 (p=6.5e-13); Schulz/hippocampus +0.48 (p=1.3e-22); universe from 450K manifest |
| [`04_cross_region_sharing/`](04_cross_region_sharing/) | 4a | 2026-08-09 | **Pass — audit required** | Direction concordance 1.000×3 and OR 91/1237/22 look degenerate; also still reports caudate≠N `True` |
| [`05_donor_group_comparison/`](05_donor_group_comparison/) | 4b | 2026-08-09 | **Pass** | ρ 0.51 / 0.12 / 0.30; MAF/LD-matched now 2/3 (hip Δ 0.079→0.314) |
| [`06_transcription_splicing_integration/`](06_transcription_splicing_integration/) | 5 | 2026-08-08 | **Pass** | expression 3/3; PSI 2/3 |
| [`07_repeat_mappability_sensitivity/`](07_repeat_mappability_sensitivity/) | 6 | 2026-08-08 | **Pass (3/4 claims 3/3)** | LINE/L1 2/3; DLPFC LINE/L1 null under cellPC (p=0.60) and reversed under high-mappability (0.864) |
| [`08_schizophrenia_risk_application/`](08_schizophrenia_risk_application/) | 7 | 2026-08-08 | **Decision stale** | `PHASE7_DECISION.md` predates module `10` summarize and cites the failed caudate≠N criterion; matched-permutation arm p=0.224 |
| [`09_libd_eqtl_mapping/`](09_libd_eqtl_mapping/) | support | — | **Debug** | BrainSeq/LIBD AA gene eQTL; QC-failed; see `EQTL_DEBUG_TODO.md` |
| [`10_downsampling_caudate/`](10_downsampling_caudate/) | Exp3 formal | 2026-08-09 | **FAIL** | Common-universe rate ratio 0.943 / 0.928; `criterion_not_solely_sample_size=False` |
| [`11_celltype_compartment_sensitivity/`](11_celltype_compartment_sensitivity/) | 6 cell-type | 2026-08-08 | **Pass** | `keep_main_figure_with_cell_composition_row`; caudate + ≥1 region only |
| [`12_supplementary_data/`](12_supplementary_data/) | packaging | 2026-08-09 | **Repackage needed** | 6 archives packaged; they encode the failing `03` and `10` results; run log is 0 bytes |

Shared helpers: [`_lib/`](_lib/).
Config: [`config/`](../config/).
Data dictionary: [`inputs/data_dictionary/`](../inputs/data_dictionary/).

## Open work

1. ~~Phase 3 fix~~ — **done 2026-08-15.** The Jaffe GEO "allPairs" file is a
   significant-results table, not a tested universe; the assayed universe is now built
   from the 450K manifest. Criterion 5 = 2/2. See
   [`03_external_meqtl_validation/PHASE3_DIAGNOSIS.md`](03_external_meqtl_validation/PHASE3_DIAGNOSIS.md).
2. **Retire the caudate≠N claim** in modules `04`, `05` (`phase4_claim_summary.tsv`,
   written by the same summarizer) and `08` (`PHASE7_DECISION.md`); they read module
   `04`'s BH-FDR lead-retention snapshot, not module `10`'s official permutation-FDR
   result.
3. **Audit module `04` concordance denominators** before any figure uses them.
4. **Commit the repair** — ~45 files / ~1,800 insertions are uncommitted.
5. **Repackage supplementary data** after 1–3.

## Verification

```bash
python3 -m unittest discover -s meqtl-validation/tests -p 'test_*.py'
```

Currently 6 tests, all passing. Coverage is light relative to the ten repaired
invariants in [`REPAIR_V2.md`](REPAIR_V2.md); passing tests are not evidence that a
claim holds.

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

Full dependency-ordered rerun sequence: [`REPAIR_V2.md`](REPAIR_V2.md) §"Recommended
production rerun order".
