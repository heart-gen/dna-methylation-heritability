# Phase 1 lock TODO — covariate / latent-factor sensitivity

**Status:** **CLOSED** — locked to **M3a** (2026-08-01)  
**Date opened:** 2026-07-30  
**Closed:** 2026-08-01  

## Decision

Primary model = **M3a** (`agedeath + sex + primarydx + snpPC1–5 + methPC1–5`) for caudate, DLPFC, and hippocampus.

Rationale: directional λ_NS improvement vs M0 and increased FDR-significant CpG discoveries in all three regions, with retained external meQTL enrichment. See `PHASE1_LOCK_DECISION.md`.

## Completed

1. Latent estimation (methPC k=5/10/15) for all regions
2. Sensitivity TensorQTL for M1/M2/M3a/M3b/M3c
3. Cross-region comparison summaries
4. Lock decision recorded and primary outputs promoted from M3a
5. `config/covariates.yml` updated (`lock_status: locked_M3a`)

## Remaining follow-up (downstream, not Phase 1 lock)

1. ~~Recompute Phase 2 VMR meQTL-burden on promoted M3a calls~~ **done 2026-08-01**
2. Spot-check Phase 3–6 consumers of primary lead tables
3. Optional manuscript panel: M0 vs M3a λ_NS / n_sig

## Pointers

- Decision: `PHASE1_LOCK_DECISION.md`
- Plan: `PHASE1_COVARIATE_SENSITIVITY_PLAN.md`
- Choice TSV: `_m/phase1_lock_choice.tsv`
- Promote log: `_m/phase1_m3a_promote.log`
- Config: `config/covariates.yml`
