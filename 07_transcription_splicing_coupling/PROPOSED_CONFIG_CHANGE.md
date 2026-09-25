# Proposed change to a pi_locked config — `config/transcription_splicing.yml`

Branch `module07/psi-identifier-join`, 2026-09-23. **Not applied.** AGENTS.md §12:
agents may recommend, never silently decide.

## Status: the code fix does NOT depend on this

The repair on this branch is correct under the config exactly as it stands, and
the module should be rerun without waiting for a PI decision here. This file
records a structural recommendation, not a blocker.

## What is wrong with the locked keys

```yaml
annotation:
  gene: inputs/counts/gene-annotation.tsv
  psi: inputs/counts/psi-annotation.tsv
```

`annotation.psi` names **one** path for a table that is **per region**, and that
path is a tracked git symlink into the caudate delivery:

```
inputs/counts/psi-annotation.tsv -> /projects/b1213/resources/processed-data/\
r-variables/caudate/_m/psi-annotation.tsv
```

The identifier that table supplies, `psi_uid`, is a row position (`p{row index}`),
and the three regions order the same 690,907 events differently, so the caudate
table renames every event for the other two regions. It joins successfully
anyway, which is why the defect was silent for three accepted runs. Detail and
evidence: `07_transcription_splicing_coupling/README.md` and
`_h/psi_features.R`.

`annotation.gene` is a single path for a table that genuinely *is* single — the
gene annotation is byte-identical across the three deliveries (one MD5) — so that
key is correct as written and is left alone.

## What the code does instead, and why that is not the whole answer

`_h/psi_features.R::psi_annotation_path()` derives the region's annotation from
the per-region entry the config **does** carry, `assay_files.psi.{region}`: the
annotation that describes an assay is delivered in the same directory as that
assay. There is deliberately no fall-back to `annotation.psi`. Stage 02 then
verifies the resolved table against the assay's own metadata and stops if it does
not match, so a wrong resolution is fatal rather than silent.

That is sound, and it is also a *convention* — "the annotation sits beside the
assay" — rather than something the configuration states. `annotation.psi` is now
an unused key that still reads as though it were the PSI annotation the module
uses, which is the kind of divergence that invites the next person to wire it
back in.

## Proposed diff

```diff
--- a/config/transcription_splicing.yml
+++ b/config/transcription_splicing.yml
@@
 annotation:
   gene: inputs/counts/gene-annotation.tsv
-  psi: inputs/counts/psi-annotation.tsv
+  # PER REGION, and not optional. psi_uid is a ROW POSITION (p0, p1, ... in file
+  # order) and the three deliveries order the same 690,907 events differently, so
+  # one shared path renames every event for two of the three regions without
+  # failing to join. A single path here was the cause of the 2026-09-23 PSI
+  # identifier defect; see 07_transcription_splicing_coupling/README.md.
+  psi:
+    caudate: inputs/counts/psi-annotation.caudate.tsv
+    dlpfc: inputs/counts/psi-annotation.dlpfc.tsv
+    hippocampus: inputs/counts/psi-annotation.hippocampus.tsv
```

Applying it needs three new per-region symlinks under `inputs/counts/`, replacing
the one that is there:

```
psi-annotation.caudate.tsv     -> .../r-variables/caudate/_m/psi-annotation.tsv
psi-annotation.dlpfc.tsv       -> .../r-variables/dlpfc/_m/psi-annotation.tsv
psi-annotation.hippocampus.tsv -> .../r-variables/hippocampus/_m/psi-annotation.tsv
```

**Nothing under `inputs/` has been touched on this branch.** That directory is
shared, and the existing `inputs/counts/psi-annotation.tsv` symlink still has one
other consumer: the withdrawn v1 tree
`local-snp-prediction/{all_individuals,BA_only}/tissue_comparison/regulatory_context/_h/00.regulatory_context_utils.R:652`
and `.../04.feature_proximity.R:74`, which carry the same defect. Repointing or
removing the symlink would change what those legacy scripts read, so it is a PI
call, not a side effect of this fix.

## Recommendation

Either apply the diff above, or leave `annotation.psi` in place and record in the
config that it is **not** the module's PSI annotation source. The second option
costs nothing and removes the trap; the first is the honest structure. What should
not happen is `annotation.psi` staying as it reads today, because it looks usable
and is not.

---

# Second proposal — declare the coupling-test FDR family

Branch `module07/abc-power-floor-fdr-family`, 2026-09-25. **Not applied.**

## Status: the code fix does NOT depend on this either

The power-floor fix on this branch enforces keys that are **already locked**:
`gates.min_vmrs_tested: 500` and `gates.min_pairs_tested: 1000`. No key was
added, removed or reinterpreted. `expression_abc` links 243-305 VMRs across the
three AA regions, so it has always been below the locked floor, and
`03_test_coupling.R` has always warned about it. What changed is that the floor
is now evaluated **per modality**, which is how the key is written, instead of
against `max()` across modalities — where `expression_nearest_gene`'s ~10,000
VMRs satisfied it for everything else and the floor bound nothing.

So the exclusion needs no PI config act. It applies a PI decision already on the
books.

## What is genuinely missing from config

**The coupling-test FDR family is nowhere declared.** `association.fdr_family`
is declared, as `modality_within_cell`, but that key governs the **pair-level**
FDR — which VMR-to-feature links count as significant, and so what
`any_sig_fdr` means. The family for the nine coupling tests themselves
(3 predictors x 3 modalities) is implicit in one line of code:

```r
res[, q := p.adjust(p, method = ts$association$fdr_method)]
```

It even borrows `association$fdr_method`, because there is no `coupling` key to
borrow from. A reader of `config/` cannot tell what the coupling-test family is,
and `00_new_run.R` cannot refuse a run whose code and declared family disagree.

Proposed addition under the existing `coupling:` block:

```yaml
coupling:
  # The FDR family for the coupling tests themselves: predictors x modalities
  # within one cell, restricted to modalities that clear gates.min_vmrs_tested
  # and gates.min_pairs_tested. A modality below either floor is reported with a
  # nominal p and no q-value, so it neither spends multiplicity nor contributes
  # evidence -- the pattern 06_partitioned_heritability uses for its membership
  # tau. Regions are never pooled.
  fdr_method: BH
  fdr_threshold: 0.05
  fdr_family: powered_modalities_x_predictors_within_cell
  fdr_family_excludes_underpowered_modalities: true
```

`fdr_method` and `fdr_threshold` are restated here rather than inherited, so the
coupling family can be changed without silently changing the pair-level family.
Their values are the ones the code already uses, so restating them changes no
result.

## What the PI should know before applying it

The restriction **changes reported significance**, in both directions, and the
direction is not uniformly conservative:

- Four results that were counted as significant are withdrawn. Three are
  `expression_abc / any_meqtl_support` in caudate, DLPFC and hippocampus — all
  three completely separated (empty off-cell, boundary estimate ~17.6, SE ~0.4,
  p underflowing to 0 in caudate), so none was ever an estimate. The fourth is
  hippocampus `expression_abc / local_genetic_control` at q = 0.022, which looks
  ordinary but sits in a modality whose minimum detectable odds ratio is 3.16.
- Two results **become** significant, both DLPFC PSI: `meqtl_proportion` and
  `any_meqtl_support`, each q 0.0546 -> 0.0425. Shrinking the family from nine
  tests to six makes BH less conservative for marginal tests. Both are marginal
  and should be written as such, not as findings.

That second bullet is why this belongs in front of the PI rather than in a commit
message. AGENTS.md 10.3 says FDR families are not combined after inspection, and
shrinking one after inspection raises the same concern. The defence is that the
exclusion criterion is a **locked floor applied to a modality's SIZE** — a
property of the ABC enhancer-gene link set, fixed before any test was run and
independent of every p-value — and not a criterion chosen because of what the
tests returned. The timing is nonetheless post-inspection, and the README should
say so plainly rather than present six tests as if nine had never been computed.

`results/coupling-power-analysis.tsv` (new stage 08) is the supporting evidence:
minimum detectable odds ratio at 80% power, alpha 0.05, per modality per
predictor, computed from the design's margins rather than its results.
