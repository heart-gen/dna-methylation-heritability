# Joint-PVE validation-independence recovery lock

**Locked:** 2026-08-21, before generating or running the replacement
validation simulations.

The validation manifest in `lgv-joint-pve-validate-20260820` is not eligible
for terminal acceptance. Although it shares no scenario seeds with the frozen
joint-estimator development run, 384 of its simulation seeds occur in the
earlier, inspected estimator-settings screen
`estimator-screen-20260818/config/screen-scenarios.tsv`. Consequently, the
associated decision run `lgv-joint-pve-decision-20260821` is retained only as
a forensic audit and is not a scientifically independent validation.

This recovery is not another estimator-development iteration. The following
objects remain frozen exactly as they existed before the contaminated
validation was evaluated:

- model family, features, transformations, coefficients, and model checksum;
- expanded simulation design and number of validation scenarios;
- all 14 acceptance gates and thresholds;
- fitting/calibration partitions, weights, and interval construction;
- BSLMM, elastic-net, HE, and effective-rank settings.

The only scientific-design change is replacement of the validation seed
offset from `500000000` to `700000000`. Before submission, the complete
12,960-row replacement manifest must have zero seed overlap with every prior
seed-bearing simulation manifest under the active and legacy Module 02 run
trees. The audit must be written into the immutable validation run.

The replacement validation run is `lgv-joint-pve-validate-20260821a`; the
derived decision run is `lgv-joint-pve-decision-20260821a`. The same binary
terminal rule applies without alteration:

- all hard gates pass: use the frozen estimator as simulation-calibrated local
  cis-SNP PVE;
- any hard gate fails: end absolute-PVE estimator development and pivot the
  manuscript to a relative/local-genetic-control axis.

No result from `lgv-joint-pve-validate-20260820` may be used to tune, replace,
or reinterpret the frozen estimator or gates.
