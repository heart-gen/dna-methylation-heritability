# Proposed change to `config/partitioned_heritability.yml` (pi_locked)

**Raised by:** Module 06 two-annotation model, `module06/two-annotation-model`
**Date:** 2026-09-23
**Status:** proposed, not applied. `config/partitioned_heritability.yml` carries
`pi_locked: true`, so per AGENTS.md §12 an agent may recommend this and must not
make it.

## Is the change required?

**No — the code is correct and complete without it.** The two-annotation model is
implemented and tested against the config exactly as it stands today. No key was
added, removed or reinterpreted, and no default was invented to stand in for a
missing key.

The annotation set is a **code-level contract** in
`06_partitioned_heritability/_h/annotations.py`: the two names, their column
order, their roles, which one carries the primary hypothesis, and which one's
enrichment is interpretable. Every stage reads it from that one file, and the
file is snapshotted into `runs/{RUN_ID}/code/_h` with the rest of `_h/`, so each
run records the annotation set it used. Stage 07 additionally stamps
`annotations_in_model`, `fdr_family_annotation` and
`tau_conditional_on_vmr_membership` onto `partitioned-h2-decision.tsv` from the
observed metrics, so the estimand is recoverable from a run's outputs alone
without consulting either the code or the config.

What the config change would buy is that the estimand becomes **declared** rather
than **implemented**: a PI reading `config/` would see which annotations are in
the model and which one the hypothesis rides on, and `00_new_run.R` could refuse
a run whose code contract and declared contract disagree. That is the reason
`writing-notes/ISSUE_06_two_annotation_model.md` asks for it under **Work**:

> `config/partitioned_heritability.yml`: declare both annotations explicitly and
> record which one carries the primary hypothesis.

It is a provenance improvement, not a correctness fix.

## Exact proposed diff

```diff
--- a/config/partitioned_heritability.yml
+++ b/config/partitioned_heritability.yml
@@ -15,9 +15,34 @@
 # ---------------------------------------------------------------- annotation
 annotation:
   # The only admissible predictor. Module 02's terminal decision is
   # PASS_RELATIVE_GENETIC_CONTROL_FAIL_ABSOLUTE_LOCUS_PVE, so the rank score is
   # usable and the PVE magnitude is not.
   score_column: local_snp_contribution_score_z
   continuous: true
+  # The S-LDSC model carries TWO annotations of ours, in this order. The order
+  # is positional in the .annot.gz, the .l2.ldscore.gz score columns and the
+  # .l2.M_5_50 entries; downstream identification is by name.
+  #
+  # One column cannot ask this module's question. local_snp_contribution_score_z
+  # is a standardized within-cell midrank percentile, so zero is the value of a
+  # median-ranked VMR -- the same value a SNP in no VMR would receive. Without a
+  # membership term there is nothing to absorb a VMR-versus-genome difference,
+  # and tau blurs the within-VMR gradient with a membership effect. With
+  # vmr_tested in the model, tau on the score is conditional on membership.
+  # See writing-notes/ISSUE_06_two_annotation_model.md.
+  annotations:
+    - name: VMR_TESTED
+      role: membership
+      form: binary
+      # Derived from the same overlap test as the score, from no threshold on
+      # it. This is NOT the banned v1 partition: it indicates membership in the
+      # tested universe, not a class of local genetic control.
+      description: SNP overlaps a VMR in the tested universe
+      enrichment_interpretable: true
+    - name: LOCAL_SNP_CONTRIBUTION_Z
+      role: score
+      form: continuous
+      description: relative rank of local SNP contribution, within cohort by region
+      # LDSC's enrichment denominator is the annotation's total value, which for
+      # a signed score is a signed sum rather than a share.
+      enrichment_interpretable: false
+  primary_hypothesis_annotation: LOCAL_SNP_CONTRIBUTION_Z
+  # PRESPECIFIED before any two-annotation result was seen: the frozen trait
+  # family stays on the score annotation only, so the family size does not
+  # change and no already-computed q-value is revised. The membership tau is
+  # reported descriptively, with a nominal p and no q.
+  fdr_family_annotation: LOCAL_SNP_CONTRIBUTION_Z
   # Guardrails asserted in code, not merely documented.
   forbid_thresholding: true
   forbid_grouping: true
```

## What would change in the code if this is accepted

Additive and small; nothing already written would be rewritten.

- `_h/annotations.py`: add a `check_against_config(cfg)` that compares
  `ANNOT_COLUMNS`, the roles, `PRIMARY_ANNOT` and `ENRICHMENT_INTERPRETABLE`
  against `cfg["annotation"]`, and aborts on any disagreement. Stages 03, 05a and
  06 already load the config, so each calls it where it currently calls
  `assert_score_source_column()`.
- `_h/00_new_run.R`: extend the existing pre-flight checks to assert
  `annotation.primary_hypothesis_annotation` and
  `annotation.fdr_family_annotation` are both `LOCAL_SNP_CONTRIBUTION_Z`, beside
  the `score_column` check already there, so a drifted config stops the run
  before a run directory exists.
- `_h/07_fdr_and_gates.R`: read the FDR-family annotation's role from config
  rather than from the hard-coded `"score"`, keeping the current value as the
  only accepted one.

None of this is needed for the fix to be correct; it is what turns the code
contract into a declared one.

## Note for whoever applies it

`config_partitioned_heritability_sha256` is recorded in every run manifest, so
editing this file changes the checksum and the three 2026-09-08 runs will no
longer match the live config. That is already true of them for a larger reason —
their estimand is withdrawn and they need rerunning (see
`06_partitioned_heritability/README.md`, **Accepted runs**) — so the natural
order is: apply the config change, then rerun all three cells, rather than
rerunning first.
