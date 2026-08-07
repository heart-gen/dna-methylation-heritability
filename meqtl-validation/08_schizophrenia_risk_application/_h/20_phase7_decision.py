#!/usr/bin/env python3
"""Consolidate Phase 7 claim snapshots into retain / supplement / omit decision.

Implements the prespecified decision rules for the schizophrenia application.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
MODULE_M = PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m"


def _read_one(path: Path) -> dict:
    if not path.exists():
        return {"_missing": True, "path": str(path)}
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return {"_missing": True, "path": str(path)}
    return df.iloc[0].to_dict()


def _truthy(val) -> bool:
    if isinstance(val, bool):
        return val
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return False
    return str(val).strip().lower() in {"true", "1", "yes", "pass", "t"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=str(MODULE_M / "decision"))
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    arch_c = _read_one(MODULE_M / "caudate/architecture_decision_snapshot.tsv")
    tx_c = _read_one(MODULE_M / "caudate/tx_decision_snapshot.tsv")
    down = _read_one(MODULE_M / "caudate_downsample/downsample_claim_snapshot.tsv")
    gx = _read_one(MODULE_M / "gxregion/gxregion_claim_snapshot.tsv")
    l3 = _read_one(MODULE_M / "level3/level3_claim_snapshot.tsv")
    dx = _read_one(MODULE_M / "diagnosis/diagnosis_claim_snapshot.tsv")
    xreg = _read_one(MODULE_M / "cross_region/cross_region_claim_snapshot.tsv")
    phase7 = MODULE_M / "phase7_cross_region_summary.tsv"
    n_caudate_loci = 0
    if phase7.exists():
        p7 = pd.read_csv(phase7, sep="\t")
        row = p7[p7["region"].astype(str).str.lower() == "caudate"]
        if len(row):
            n_caudate_loci = int(row.iloc[0].get("n_loci", row.iloc[0].get("n_sig_loci", 0)))

    if not arch_c.get("_missing"):
        n_caudate_loci = int(float(arch_c.get("n_sig_risk_loci", n_caudate_loci) or 0))
    elif phase7.exists():
        p7 = pd.read_csv(phase7, sep="\t")
        crow = p7[p7["region"].astype(str).str.lower() == "caudate"]
        if len(crow):
            for col in ["n_sig_loci", "n_loci"]:
                if col in crow.columns:
                    n_caudate_loci = int(float(crow.iloc[0][col]))
                    break

    crit = {
        "c1_ge1_risk_locus_meqtl": _truthy(arch_c.get("criterion_ge1_locus_meqtl", False))
        or n_caudate_loci >= 1,
        "c2_predictability_enrichment_caudate": _truthy(
            arch_c.get("criterion_enrich_higher_predictability", False)
        ),
        "c3_tx_coupling_ge1": _truthy(tx_c.get("criterion_ge1_locus_with_tx", False))
        or int(float(tx_c.get("n_sig_vmr_with_tx_coupling", 0) or 0)) >= 1,
        "c4_caudate_not_solely_N": _truthy(down.get("criterion_not_solely_sample_size", False)),
        "c5_level3_shared_genetic_support": _truthy(l3.get("pass", False)),
        "c6_diagnosis_optional": (not dx.get("_missing", False)),
    }
    primary_pass = all(
        [
            crit["c1_ge1_risk_locus_meqtl"],
            crit["c2_predictability_enrichment_caudate"],
            crit["c3_tx_coupling_ge1"],
            crit["c4_caudate_not_solely_N"],
            crit["c5_level3_shared_genetic_support"],
        ]
    )

    dx_positive = False
    if not dx.get("_missing"):
        dx_positive = _truthy(dx.get("any_diagnosis_association_fdr", False)) or _truthy(
            dx.get("pass_optional_disease_layer", False)
        )

    gx_frac = float(gx.get("frac_sig_joint_interaction", 0) or 0) if not gx.get("_missing") else float("nan")

    n_best_dx = int(float(dx.get("n_best_vmr_primary_fdr_sig", 0) or 0)) if not dx.get("_missing") else 0
    n_any_dx = int(float(dx.get("n_any_vmr_primary_fdr_sig", 0) or 0)) if not dx.get("_missing") else 0
    n_tx_dx = int(float(dx.get("n_transcript_primary_fdr_sig", 0) or 0)) if not dx.get("_missing") else 0
    dx_strong = (n_best_dx + n_tx_dx) > 0

    if primary_pass and dx_strong:
        decision = "retain_main_text"
        rationale = (
            "Meets §12.12 architecture criteria (meQTL, predictability enrichment, TX coupling, "
            "N-matched caudate excess, Level-3 shared genetic support) and shows diagnosis association "
            "at a best VMR or linked transcript feature."
        )
    elif primary_pass:
        decision = "retain_main_text_proof_of_application"
        dx_note = (
            f" Diagnosis FDR hits among secondary VMRs only (n_any={n_any_dx}, n_best={n_best_dx}); "
            "treat as supportive exploratory case–control signal, not the primary claim."
            if dx_positive and not dx_strong
            else " Diagnosis layer null or weak at best VMR/transcript features (allowed)."
        )
        rationale = (
            "Meets §12.12/§12.13 architecture criteria for a focused caudate application."
            + dx_note
            + " Present as regulatory application, not mediation. "
            "Level-3 relies on GTEx; internal LIBD eQTL map remains QC-failed for genome-wide discovery."
        )
    elif crit["c1_ge1_risk_locus_meqtl"] and crit["c2_predictability_enrichment_caudate"]:
        decision = "supplement_proof_of_application"
        rationale = (
            "Partial success: risk-meQTL and predictability enrichment present, but missing TX, "
            "N-matched regional, or Level-3 support for full main-text retention."
        )
    else:
        decision = "omit_or_minimal_supplement"
        rationale = "Stop criteria approached: insufficient multi-omic SCZ application support."

    rows = [
        {
            "decision": decision,
            "primary_success_criteria_pass": primary_pass,
            "n_caudate_sig_loci": n_caudate_loci,
            **{k: v for k, v in crit.items()},
            "diagnosis_fdr_positive": dx_positive,
            "gxregion_frac_joint_interaction": gx_frac,
            "level3_n_pass": l3.get("n_loci_level3_pass", ""),
            "level3_via": "GTEx" if int(float(l3.get("n_loci_gtex_caudate_gene_variant", 0) or 0)) > 0 else "",
            "downsample_median_n_sig_pairs": down.get("downsample_median_n_sig_pairs", ""),
            "rationale": rationale,
            "allowed_language": (
                "Schizophrenia-risk variants associate with CpG methylation in genetically anchored "
                "caudate VMRs, including loci with transcriptional coupling and external eQTL support; "
                "not evidence that methylation mediates schizophrenia risk."
            ),
            "forbidden_language": (
                "methylation mediates SCZ risk; causal variants from elastic-net; ancestry-specific "
                "SCZ methylation; formal colocalization (not run)."
            ),
        }
    ]
    write_tsv(outdir / "phase7_retain_decision.tsv", rows)

    # Markdown memo
    md = outdir / "PHASE7_DECISION.md"
    r = rows[0]
    md.write_text(
        f"""# Phase 7 retain / supplement decision

**Decision:** `{r['decision']}`

## Criteria (§12.12)

| Criterion | Pass |
|---|---|
| ≥1 SCZ-risk locus with CpG meQTL | {r['c1_ge1_risk_locus_meqtl']} (n_loci={r['n_caudate_sig_loci']}) |
| Predictability enrichment (caudate) | {r['c2_predictability_enrichment_caudate']} |
| ≥1 locus with meth–TX coupling | {r['c3_tx_coupling_ge1']} |
| Caudate not solely sample size | {r['c4_caudate_not_solely_N']} |
| Level-3 shared genetic support | {r['c5_level3_shared_genetic_support']} (n_pass={r['level3_n_pass']}; via {r['level3_via'] or 'n/a'}) |
| Diagnosis layer (optional) | tested={r['c6_diagnosis_optional']}; FDR-positive={r['diagnosis_fdr_positive']} |

## Rationale

{r['rationale']}

## Language

**Allowed:** {r['allowed_language']}

**Forbidden:** {r['forbidden_language']}

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
"""
    )
    print(f"Decision: {decision}")
    print(f"Wrote {outdir / 'phase7_retain_decision.tsv'} and {md}")


if __name__ == "__main__":
    main()
