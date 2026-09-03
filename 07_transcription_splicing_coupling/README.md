# 07_transcription_splicing_coupling — regulatory consequences

Tests whether meQTL-supported or locally controlled VMRs are more likely to have existing significant associations with gene/transcript abundance or transcript usage/splicing.

**Status: implemented, production run 2026-09-02, awaiting PI acceptance.**
Gated on `05_cpg_meqtl_burden` acceptance (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID"), which is met
(`cmb-AA-*-20260825`). Runs `tsc-AA-{caudate,dlpfc,hippocampus}-20260902` all
returned `PASS_TX_COUPLING_QC`. No run is accepted yet — see **Accepted runs**.

## Migrating from

`meqtl-validation/06_transcription_splicing_integration/` and `meqtl-validation/09_libd_eqtl_mapping/`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Scope discipline

Do **not** initiate an unbounded transcriptome-wide fishing analysis. Reuse the
prespecified expression and splicing analyses and document their tested universe
explicitly. Adjust for the number of tested features, VMR length, VMR-to-feature
distance, methylation variance, local SNP number, and applicable technical
factors.

Allowed: "Genetically regulated VMRs are more frequently transcriptionally
coupled."
Forbidden: "Methylation mediates the genetic effect on expression or splicing."

## Pipeline

| Stage | Script | Purpose |
|---|---|---|
| 00 | `_h/00_new_run.R` | Gate on accepted 01, 02 and 05; require one shared `vmr_set_id`. |
| 01 | `_h/01_build_feature_links.R` | Rebuild VMR→feature links on accepted VMRs; write the tested universe. |
| 02 | `_h/02_run_local_associations.R` | Fit `feature ~ VMR methylation + covariates` per pair. |
| 03 | `_h/03_test_coupling.R` | The three coupling tests per modality. |
| 04 | `_h/04_apply_gates.R` | Acceptance gate and interpretation constraints. |
| 05 | `_h/05_plot.py` | Figures. |
| 06 | `_h/06_finalize_run.R` | Seal the run (sealing is not acceptance). |

Submit one cell with `_h/submit_transcription_splicing.sh <cohort> <region>`;
`DRY_RUN=1` prints the job graph, `SMOKE_N=1` permits unaccepted upstreams.

## Why the legacy tables could not be reused

The legacy coupling script consumed `architecture_model_input.tsv` from
`local-snp-prediction/.../regulatory_context/_m/`. Those tables are keyed to
**pre-repair VMRs** and carry `h2_category`, `r_squared_cv` and `h2_unscaled` as
columns — the first two banned by AGENTS.md §3, the third by Module 02's
terminal decision. AGENTS.md §6 also forbids carrying downstream numbers across
VMR turnover, and a VMR's methylation summary is a function of its boundary. So
the links are rebuilt and every pair is refitted; only the *method* is reused.

## Pair-level model

`feature ~ VMR_methylation + Age + Sex + RIN + MoD + mito_mapping_rate +
percent_assigned + cell proportions (asin-sqrt)`.

There are ~2.7M pairs. Because every pair within a modality shares the same
donors and covariate matrix, the covariates are residualised out once by QR and
each pair reduces to a dot product. This is algebraically identical to fitting
the full model per pair — verified against `lm()`, matching t and p to 4+
significant figures — not an approximation.

### PSI missingness restricts the tested universe

PSI events are frequently unquantified in a subset of donors: measured on
caudate, the median event is NA in 62% of donors and only ~17% are quantified in
every donor. An event is tested only where it is quantified in every donor of
the analysis set (`normalisation.psi.max_na_fraction`), so every pair is fitted
on the same complete design. **This biases the retained PSI set toward
constitutively quantified events and must be stated in Methods.** The declared
and realised universes are both written out
(`results/tested-universe.tsv`, `results/{modality}-realised-universe.tsv`).

## Scope boundary: the internal LIBD eQTL map is not used

`meqtl-validation/09_libd_eqtl_mapping/` is **out** of this module's acceptance
gate. Its genome-wide QC repair is open (~1–2 eGenes at FDR 0.05; see that
directory's `EQTL_DEBUG_TODO.md`, tasks A1–D1 unchecked). The coupling analysis
does not depend on it — per AGENTS.md §7.6 it reuses the prespecified local
association screen rather than running a transcriptome-wide discovery. Enabling
`internal_libd_eqtl_support_arm` while that repair is open fails the gate.

## Acceptance gate

1. accepted 01, 02 and 05 runs for the cell, sharing one `vmr_set_id`;
2. every enabled modality ran;
3. tested universe above the configured minimum VMR and pair counts;
4. at least one test produced a finite estimate;
5. no banned column reached a model frame;
6. the LIBD eQTL arm is off;
7. Stage 04 decision `PASS_TX_COUPLING_QC`;
8. immutable Stage 06 checksums and a manual README acceptance record.

A null coupling result is a reportable finding, not a gate failure.

## Accepted runs

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| _(none)_ | | | | | | | |

AGENTS.md §6 makes acceptance a human step: no row appears here until a
production run's gate stage passes and the PI records it. A completed SLURM job
is not acceptance.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
