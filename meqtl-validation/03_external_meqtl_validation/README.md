# Phase 3: External public brain meQTL validation

## Question

Are CpGs / VMRs with higher local genetic predictability more likely to have
independent published brain cis-meQTL evidence?

## Catalog

See `inputs/data_dictionary/_m/public_meqtl_resources.tsv` and
`inputs/data_dictionary/_h/README.md`.

| Resource | Role | Status |
|---|---|---|
| `jaffe_dlpfc_450k_meqtl` | Primary external, DLPFC-matched (450K) | Universe rebuilt from array manifest (2026-08-15) |
| `schulz_hippocampus_array_meqtl` | Primary external, hippocampus-matched (450K) | Restored to primary; universe from array manifest |
| `brainseq_wgbs_meqtl` | **Not external** — cohort overlaps this discovery data (BrainSeq/LIBD) | Synapse deferred; even if obtained, do not count as independent Phase 3 validation |
| `brainseq_wgbs_meqtl_scz_subset` | Exploratory only (Nature SCZ-risk tables) | WGBS, no genome-wide universe → `not_estimable` |

**Note:** BrainSeq full catalogs are same-/overlapping-cohort with this project's WGBS+genotype donors. Prefer Jaffe and Schulz for independent external support. See `_m/harmonized/brainseq_wgbs_meqtl.PENDING_SYNAPSE.txt`.


## Run

```bash
cd meqtl-validation/03_external_meqtl_validation/_m
mkdir -p logs
sbatch ../_h/step_1.sh          # init / download checklist
sbatch ../_h/step_0.sh          # download + hg38 harmonize + VMR overlap
# After Phase 2 burden tables exist:
sbatch ../_h/step_2.sh          # enrichment models + summary
# Or locally:
#   conda activate /projects/p32505/opt/envs/genomics
#   python3 ../_h/04_run_and_summarize.py
```

Primary model: VMR any-external-support ~ continuous predictability (+ n CpGs / tech covariates).  
Preferred tissue pairings are `primary`; other region×resource tests are `secondary_cross_region`.

## Assayed universe (repair-v2 + 2026-08-15 fix)

Only CpGs documented as assayed by a resource enter its denominator. **The assayed
universe is the Illumina 450K manifest**, not the published results table.

This corrects two opposite errors:

| Version | Denominator | Jaffe/DLPFC VMR support rate | Problem |
|---|---|---:|---|
| Pre-repair (≤2026-07-30) | all WGBS VMR CpGs | 0.24 | Counted never-assayed WGBS CpGs as negatives |
| Repair-v2 (2026-08-08) | results-table rows only | 1.000 | Only positives; outcome constant → perfect separation |
| **Current (2026-08-15)** | **450K manifest ∩ VMR CpGs** | **0.625** | — |

Jaffe and Schulz are both 450K studies (100.0% and 99.9% of their significant probes
are on the manifest), so both are eligible for the tissue-matched primary test.
`build_array_resource()` refuses any resource whose significant probes fall below 95%
manifest coverage, so a non-array catalog cannot be given manufactured negatives.
The lifted manifest is cached at `_m/support/450k_universe_hg38.tsv.gz`
(485,441/485,512 probes lift to hg38).

BrainSeq is WGBS with no genome-wide universe and remains `not_estimable` — and is
same-cohort, so it would be exploratory regardless.

Full root-cause analysis: [`PHASE3_DIAGNOSIS.md`](PHASE3_DIAGNOSIS.md).

**Caveat for Methods:** only ~4% of WGBS VMR CpGs are on the 450K array (DLPFC
6,658/154,325). The array under-samples distal intergenic sequence, which is where
the high-predictability VMR class lives. This test is a conservative, array-accessible
probe of the central axis.

Note: Jaffe × caudate and Schulz × non-hippocampus are `secondary_cross_region` —
interpret tissue-matched primary tests for the claim; do not pool platforms.

Outputs:

- `{region}/external_support_model_{resource}.tsv`
- `{region}/external_matched_{resource}.tsv`
- `_m/external_support_models_all.tsv`
- `_m/external_support_primary_summary.tsv`
- `_m/phase3_criterion5_verdict.tsv`
