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
