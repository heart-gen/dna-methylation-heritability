# 07_transcription_splicing_coupling — regulatory consequences

Tests whether meQTL-supported or locally controlled VMRs are more likely to have existing significant associations with gene/transcript abundance or transcript usage/splicing.

**Status: expression coupling accepted (AA, 2026-09-08). The PSI (splicing)
results of all three accepted runs are WITHDRAWN as of 2026-09-23 and require a
rerun — see "The PSI identifier join was broken" below. Expression and ABC are
unaffected and were verified so.**
Gated on `05_cpg_meqtl_burden` acceptance ("No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID"), which is met
(`cmb-AA-*-20260825`). Runs `tsc-AA-{caudate,dlpfc,hippocampus}-20260902` all
returned `PASS_TX_COUPLING_QC`. See **Accepted runs**.

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
columns — the first two banned, the third by Module 02's
terminal decision. Also forbid carrying downstream numbers across
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

## The PSI identifier join was broken (found and fixed 2026-09-23)

Two defects, both confirmed independently against the delivered files. The
evidence and the fix are documented at the top of `_h/psi_features.R`.

**1. `psi_uid` is a row position, not an identifier.** It is literally
`p{row index}` in whatever order a region's `psi-annotation.tsv` is written, and
the three regions ship the same 690,907 events in *different* orders:

| pair | psi_uids naming the same event |
|---|---|
| caudate vs dlpfc | 93,580 / 690,907 |
| caudate vs hippocampus | **0** / 690,907 |
| dlpfc vs hippocampus | 309,180 / 690,907 |

`config/transcription_splicing.yml`'s `annotation.psi` names one path,
`inputs/counts/psi-annotation.tsv`, which is a tracked git **symlink into the
caudate delivery**. Stage 01 therefore built every region's PSI links from
caudate's row order, and stage 02 joined them with `rownames(rse) %in% ...`.
Because the identifier exists in every region, the join always succeeded: no
error, no warning, no reduced feature count. So **caudate's PSI analysis is
intact**, DLPFC analysed a different event on 86.5% of its links, and
hippocampus on 100% of them. The sealed significant-pair rates are the
fingerprint: 419/353,728 (caudate), 36/177,537 (DLPFC), 9/243,226 (hippocampus).

**2. The delivered `rse-psi.*.RData` objects are internally misaligned.**
`rowData(rse)`'s columns are the region's annotation in file order
(`rowData$psi_uid` is p0, p1, …, and `rowRanges` agrees with it), while the rows
themselves carry a permutation of the same id set. A chrY event cannot be
quantified in a female donor, so which half the assay values follow was settled
from the data in all three regions: labelling by `rownames(rse)` puts the chrY
rows at 0.90–0.92 NA in females against 0.67 in males (difference +0.221 to
+0.229), labelling by `rowData` gives +0.004 to +0.049. **`rownames(rse)` is
authoritative; `rowData` is read only as a keyed lookup table, never positionally
against the rows.**

**The fix.** A PSI feature is keyed on `event_info` + `gene_id`, which is unique
(690,907/690,907; `event_info` alone is not — 308 events are annotated to two
genes). The annotation is read from beside *this region's* assay, with no
fall-back to the caudate symlink. Every link table records its `region` and
`feature_namespace`, stage 02 refuses one that is not its own, and stage 02
verifies the annotation against the assay's own metadata before fitting
anything — which also catches the case where both stages resolve the same wrong
path. An identifier that does not resolve is an error; a legacy positional link
table is refused by shape. Smoke check:
`tests/test_psi_identifier_join.R`.

**Expression and ABC are unaffected, and this was verified rather than assumed.**
`rownames(rse_gene)` equals `rowData(rse_gene)$gene_id` equals
`gene-annotation.tsv`'s `gene_id`, in identical order, in all three regions; the
gene annotation file is byte-identical across the three deliveries (one MD5); and
the identifier is a versioned Ensembl gene ID, not a row position.

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
does not depend on it — it reuses the prespecified local
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

**The PSI columns of these three runs are withdrawn (2026-09-23).** Their
expression and ABC results stand; their splicing results must be recomputed with
the repaired identifier join, and the PSI sentence below ("strong in caudate,
thin in DLPFC, null in hippocampus") is exactly the artefact the defect
produces. Caudate's PSI numbers are expected to survive largely unchanged,
because caudate is the region whose annotation was read; DLPFC and hippocampus
carry no information about splicing as run.

Permitted claim: genetically regulated VMRs are more frequently transcriptionally
coupled. Forbidden: methylation mediates the genetic effect on expression or
splicing. Nearest-gene expression supports that sentence in all three AA cells.
PSI is strong in caudate, thin in DLPFC (24 coupled VMRs), and null in
hippocampus. ABC links are underpowered and not a claim. LIBD eQTL arm is off
and is not part of this acceptance. PSI completeness filtering must be stated
in Methods.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| tsc-AA-caudate-20260902 | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-09-08 | Kynon J.M. Benjamin | PASS_TX_COUPLING_QC | 6/9 tests FDR-significant (2 local-control, 4 meQTL); PSI 227 coupled VMRs |
| tsc-AA-dlpfc-20260902 | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-09-08 | Kynon J.M. Benjamin | PASS_TX_COUPLING_QC | 5/9 tests; nearest-gene expression all three predictors; PSI 24 coupled VMRs |
| tsc-AA-hippocampus-20260902 | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-09-08 | Kynon J.M. Benjamin | PASS_TX_COUPLING_QC | 5/9 tests; nearest-gene expression all three predictors; PSI null (7 coupled VMRs) |

## Contract

This module follows: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
