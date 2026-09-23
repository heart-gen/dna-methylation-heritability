# Proposed config change — NON-BLOCKING, comment only

Branch: `module04/music-cell-adjustment` (F5 + F6 in
`04_repeat_repressive_architecture`), 2026-09-23.

## Nothing is blocked

The MuSiC/scMD fix is implemented entirely in `_h/` and `00_shared/` and is
correct under `config/repeat_annotations.yml` exactly as it stands. That config
is `pi_locked` and was not touched. It already declared, at
`sensitivities` (lines 330-334):

```yaml
  adjust_cell_composition: true
  cell_composition_rna_music: true
  cell_composition_dnam_scmd: true    # caudate only, when the integration gate passes
```

so the locked configuration asked for a MuSiC-derived adjustment in every region
and a gated scMD adjustment before the code did either. The two modality keys are
now realized in `_h/01_build_features.R` (which modality each column is built
from) and consumed in `_h/02_test_association.R` (which arm gates where).
`config/cell_deconvolution.yml`, which `covariate_sources.cell_composition_pcs`
names as the source, already carries `paths.rna_proportions_template`,
`paths.dnam_proportions_template` and the `validation` thresholds the gate reads,
so no new key was needed anywhere.

## The one thing a PI may wish to write down

The modality keys carry no explanatory comment, and the code was miscited against
them for a month: `_h/02_test_association.R` claimed the `low_cell_composition`
subset arm realized `cell_composition_rna_music`, while the column it subset on
was scMD-derived. The code comment is now correct; the config could say the same
thing so the next reader does not have to infer it from an R script.

Proposed diff (comments only — no key, value or behaviour changes):

```diff
   # Refits every model with cell_composition_r2 ADDED to the adjustment set.
   # An adjustment-variant arm, not a subset arm.
   adjust_cell_composition: true
+  # The two keys below name the DECONVOLUTION MODALITY, not an arm. Realized in
+  # _h/01_build_features.R: cell_composition_r2 is built from RNA MuSiC in every
+  # region, and both composition arms (low_cell_composition,
+  # adjust_cell_composition) therefore gate in all three. The scMD quantity is
+  # built as cell_composition_r2_scmd and drives its own arm,
+  # adjust_cell_composition_scmd, only where the integration gate passes.
   cell_composition_rna_music: true
   cell_composition_dnam_scmd: true    # caudate only, when the integration gate passes
```

If this is accepted, it changes no fit, no gate and no run: the config checksum
in a future run's manifest would differ from the 2026-09-06 runs' for a
comment-only reason, which is the only cost.

## One decision the PI may want to make explicit (also not blocking)

`cell_composition_dnam_scmd` carries no `gating:` key, and every sensitivity in
this block is gating unless it says otherwise -- `exclude_snp_proximal_cpgs` has
`descriptive_only: true` and `matched_measurability` has `gating: false`. The
implementation follows that reading: where the integration gate passes, the
`adjust_cell_composition_scmd` arm joins the survival conjunction.

The consequence is asymmetric by construction and worth seeing before the rerun:
**caudate's conjunction gains one arm that DLPFC's and hippocampus's do not
have**, so a caudate outcome must now also survive the scMD-adjusted refit. On
the accepted run caudate carries the quiescent survival (LINE/L1 is already set
aside there as technically confounded), so if the scMD refit attenuates it, a
3/3 quiescent claim could become 2/3 for a reason that exists in one region only.
That is the locked config's rule working as written, not a defect, and the claims
table now names the asymmetry in `sensitivity_arms_not_fitted`. But if the PI
intends the scMD arm to be reported and never to break a claim, the one-line
change is:

```diff
-  cell_composition_dnam_scmd: true    # caudate only, when the integration gate passes
+  cell_composition_dnam_scmd:
+    # caudate only, when the integration gate passes
+    gating: false
```

which would need a matching two-line change in `_h/02_test_association.R` to
write that arm to its own file, exactly as `exclude_snp_proximal` is handled.
No such change has been made: the code implements the config as it stands.
