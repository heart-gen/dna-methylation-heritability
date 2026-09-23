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
