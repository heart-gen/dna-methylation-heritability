# 07_transcription_splicing_coupling — regulatory consequences

Tests whether meQTL-supported or locally controlled VMRs are more likely to have existing significant associations with gene/transcript abundance or transcript usage/splicing.

**Status: accepted (AA, 2026-09-25), expression and splicing both.** The
accepted runs are `tsc-AA-{caudate,dlpfc,hippocampus}-20260925-b`. They carry the
repaired PSI identifier join and the locked per-modality power floor, and they
supersede the `-20260902` runs, whose PSI results were withdrawn 2026-09-23 — see
"The PSI identifier join was broken" below. **Splicing coupling is now a
three-region result, not a caudate-specific one**, which changes the
corresponding row of the locked analysis plan. `expression_abc` is fitted and surfaced but
sits **outside the coupling-test FDR family** on power grounds and is not a
claim; see **The ABC exclusion**.
Gated on `05_cpg_meqtl_burden` acceptance ("No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID"), which is met by
`cmb-AA-*-20260924` — the locked-covariate runs accepted 2026-09-25, which are
what these runs consume. The `-20260902` runs consumed `cmb-AA-*-20260825`, now
superseded. All three sealed `PASS_TX_COUPLING_QC`. See **Accepted runs**.

## Migrating from

`meqtl-validation/06_transcription_splicing_integration/` and `meqtl-validation/09_libd_eqtl_mapping/`.

The legacy v1 trees were retired on 2026-10-06; their tracked files are
recoverable from the annotated tag `v1-legacy-final`.

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

Permitted claim: genetically regulated VMRs are more frequently transcriptionally
coupled. Forbidden: methylation mediates the genetic effect on expression or
splicing.

Nearest-gene expression supports the permitted claim in all three AA cells, and
**since the PSI repair so does splicing** -- the pre-repair reading "strong in
caudate, thin in DLPFC, null in hippocampus" was the identifier-join defect, not
biology. Splicing coupling is now a three-region result rather than a
caudate-specific one, which changes the corresponding row of the locked analysis plan.
`expression_abc` is excluded from the coupling-test FDR family on power grounds
and is **not a claim**; see **The ABC exclusion** below. The LIBD eQTL arm is off
and is not part of this acceptance. PSI completeness filtering must be stated in
Methods.

Coupled-VMR counts below are over **all tested VMRs in the modality**. The
superseded rows' counts (227 / 24 / 7) were over the narrower association model
frame; both are defensible denominators, and the earlier table did not say which
it used. State the denominator wherever these counts appear.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| tsc-AA-caudate-20260925-b | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-09-25 | Kynon J.M. Benjamin | PASS_TX_COUPLING_QC | 6/6 eligible tests FDR-significant (2 local-control, 4 meQTL); FDR family 6, `expression_abc` excluded; 11,528 VMRs / 2,679,486 pairs; PSI 253 coupled VMRs |
| tsc-AA-dlpfc-20260925-b | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-09-25 | Kynon J.M. Benjamin | PASS_TX_COUPLING_QC | 6/6 eligible tests FDR-significant (2 local-control, 4 meQTL); FDR family 6, `expression_abc` excluded; 9,570 VMRs / 2,093,522 pairs; PSI 87 coupled VMRs (34 pre-repair) |
| tsc-AA-hippocampus-20260925-b | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-09-25 | Kynon J.M. Benjamin | PASS_TX_COUPLING_QC | 6/6 eligible tests FDR-significant (2 local-control, 4 meQTL); FDR family 6, `expression_abc` excluded; 9,495 VMRs / 2,064,462 pairs; PSI 134 coupled VMRs (8 pre-repair) |

Provenance is uniform across the three cells: `vmrcat-AA-{region}-20260816` ->
`lgv-AA-{region}-rescore-20260913` -> `cmb-AA-{region}-20260924` -> this run,
sealed 2026-09-25T15:21 at commit `e7bbcdaa7`, `smoke_run = FALSE`.

### The ABC exclusion, and why it is a restriction rather than a result

`config/transcription_splicing.yml:gates:min_vmrs_tested` is **500** and is
`pi_locked`. `expression_abc` links 305 (caudate), 250 (DLPFC) and 243
(hippocampus) VMRs, so it has been below that floor in every cell since the
module was written. It nonetheless entered the FDR family until 2026-09-25,
because `04_apply_gates.R` evaluated the floor against `max()` across modalities
and `expression_nearest_gene` links ~10,000 VMRs. A floor satisfied by the
strongest cell is not a floor.

The consequence was not cosmetic. `any_meqtl_support` crossed with the coupled
outcome is a 2x2 whose off-cell was **empty** in all three regions -- 15/305,
12/250, 8/243 coupled VMRs, every one of them with meQTL support -- so the
logistic coefficient had no maximum likelihood estimate and ran to the boundary
(estimate near 17.6, SE near 0.4, p underflowing). Sorted by p, that
unidentified coefficient was the module's single most significant "finding".

Both defects are now closed: the floor is applied per modality, and excluded
tests carry `q = NA`, `in_fdr_family = FALSE` and an
`fdr_exclusion_reason` naming the VMR count against the locked floor. The
estimates are **retained and reported**, not deleted -- the same
restriction-not-removal mechanism Module 04 uses for technically confounded
caudate. `results/coupling-power-analysis.tsv` (stage `_h/08`) carries the
justification as a minimum detectable odds ratio: ABC 4.51 / 5.52 / 7.91 against
1.25-2.04 for the powered modalities.

**Two DLPFC PSI results became significant only because the family shrank from 9
to 6**, at q = 4.25e-02 against a 0.05 threshold (`meqtl_proportion` and
`any_meqtl_support`, from q = 5.46e-02). They must be written as marginal. This
is disclosed rather than discovered later: shrinking a BH family raises every
surviving q-value's neighbours, and the locked analysis plan forbids recombining
FDR families after inspection, so the change was made on the locked power floor
and not on any view of the results.

### Two config proposals, closed 2026-09-27 without changing the lock

Both were recommendations, neither was required for correctness, and the PI
declined both. The substance is recorded here so it is not rediscovered as a
defect.

**1. `annotation.psi` is an unused key, and stays one.**
`config/transcription_splicing.yml` names a single `annotation.psi` path for a
table that is per region, and that path is a tracked git symlink into the caudate
delivery -- the origin of the identifier-join defect described above. The code no
longer reads it: `_h/psi_features.R::psi_annotation_path()` derives each region's
annotation from the per-region entry the config does carry,
`assay_files.psi.{region}`, on the convention that the annotation describing an
assay is delivered beside that assay. **There is deliberately no fall-back to
`annotation.psi`**, and stage 02 then verifies the resolved table against the
assay's own metadata and stops if they disagree, so a wrong resolution is fatal
rather than silent. `annotation.gene` is left alone because it is genuinely
single -- the gene annotation is byte-identical across the three deliveries. The
residual risk is that `annotation.psi` still *reads* as though it were the
annotation in use; this paragraph, not a config edit, is what stops the next
reader wiring it back in.

**2. The coupling-test FDR family stays declared in code.** `association.fdr_family`
(`modality_within_cell`) governs the **pair-level** FDR -- which VMR-to-feature
links are significant, and so what `any_sig_fdr` means. The family for the nine
coupling tests themselves is `_h/03_test_coupling.R:233`,
`in_fdr_family := power_eligible & is.finite(p)`, borrowing
`association$fdr_method` because there is no `coupling` key to borrow from. It is
stamped onto every run (`in_fdr_family`, `fdr_family_size`,
`modalities_excluded_from_fdr_family`), so it is recoverable from a run's outputs
without reading the code.

**One consequence binds a downstream module.** Because the family is not a config
key, a consumer cannot learn it from `config/`; it must read the `in_fdr_family`
column out of this module's tables. `08_region_donor_generalization` does not --
its harvest spec omits `outcome_role`, so `expression_abc` enters 08's claim
family although this module excluded it on power grounds. That is a Module 08
defect, recorded in that module's README, and it is the reason this proposal was
worth writing down rather than deleting.

### Superseded

`tsc-AA-{caudate,dlpfc,hippocampus}-20260902` (accepted 2026-09-08). Withdrawn on
two independent grounds: the PSI identifier-join defect (2026-09-23), and the
vacuous power floor described above. Their nearest-gene expression conclusions
are reproduced by the current runs; their PSI columns and their ABC q-values must
not be cited. Notes retained for audit: caudate 6/9 tests, PSI 227 coupled VMRs
(model-frame denominator); DLPFC 5/9, PSI 24; hippocampus 5/9, PSI null at 7.

`tsc-AA-{caudate,dlpfc,hippocampus}-20260925` -- the first rerun of the day.
It repaired PSI but still admitted `expression_abc` to the FDR family, so its
q-values are over a 9-test family and its ABC rows report a separated
coefficient as a finding. Superseded by `-20260925-b`.

`tsc-AA-{caudate,dlpfc,hippocampus}-20260925-a` -- **abandoned, never sealed.**
Killed mid-chain by a cluster-wide root cancellation at ~13:30-13:56 on
2026-09-25 that also took 3,019 Module 03 array tasks across 17+ nodes, with no
logs and `ExitCode 0:0`. No scientific content; do not cite.

## Contract

This module follows: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
