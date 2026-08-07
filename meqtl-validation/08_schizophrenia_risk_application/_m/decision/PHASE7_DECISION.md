# Phase 7 retain / supplement decision

**Decision:** `retain_main_text_proof_of_application`

## Criteria (§12.12)

| Criterion | Pass |
|---|---|
| ≥1 SCZ-risk locus with CpG meQTL | True (n_loci=31) |
| Predictability enrichment (caudate) | True |
| ≥1 locus with meth–TX coupling | True |
| Caudate not solely sample size | True |
| Level-3 shared genetic support | True (n_pass=3; via GTEx) |
| Diagnosis layer (optional) | tested=True; FDR-positive=True |

## Rationale

Meets §12.12/§12.13 architecture criteria for a focused caudate application. Diagnosis FDR hits among secondary VMRs only (n_any=1, n_best=0); treat as supportive exploratory case–control signal, not the primary claim. Present as regulatory application, not mediation. Level-3 relies on GTEx; internal LIBD eQTL map remains QC-failed for genome-wide discovery.

## Language

**Allowed:** Schizophrenia-risk variants associate with CpG methylation in genetically anchored caudate VMRs, including loci with transcriptional coupling and external eQTL support; not evidence that methylation mediates schizophrenia risk.

**Forbidden:** methylation mediates SCZ risk; causal variants from elastic-net; ancestry-specific SCZ methylation; formal colocalization (not run).

## Supporting snapshots

- Architecture: `caudate/architecture_decision_snapshot.tsv`
- TX: `caudate/tx_decision_snapshot.tsv`
- Downsample: `caudate_downsample/downsample_claim_snapshot.tsv`
- G×region: `gxregion/gxregion_claim_snapshot.tsv`
- Level 3: `level3/level3_claim_snapshot.tsv`
- Diagnosis: `diagnosis/diagnosis_claim_snapshot.tsv`
- Machine-readable: `decision/phase7_retain_decision.tsv`

## Deferred

- Level 4 coloc (requires strong dual regional QTL signal + ancestry-matched LD)
- Formal mediation / MR
- LIBD genome-wide eQTL QC repair (`09_libd_eqtl_mapping/EQTL_DEBUG_TODO.md`)
