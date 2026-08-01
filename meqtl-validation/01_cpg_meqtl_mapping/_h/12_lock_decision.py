#!/usr/bin/env python3
"""Apply Phase 1 lock go/no-go rules and write PHASE1_LOCK_DECISION.md."""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path

import pandas as pd

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
PHASE1 = PROJECT / "meqtl-validation" / "01_cpg_meqtl_mapping"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions", nargs="+", default=["caudate", "dlpfc", "hippocampus"])
    p.add_argument("--pilot-region", default="caudate")
    return p.parse_args()


def load_comparison(region: str) -> pd.DataFrame:
    path = PHASE1 / region / "_m" / "covariate_sensitivity" / "comparison_summary.tsv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def candidate_from_pilot(tab: pd.DataFrame) -> tuple[str, str]:
    """Return (model_id, reason) for caudate pilot."""
    ok = tab[tab["status"] == "ok"].copy()
    if ok.empty or "M0" not in set(ok["model_id"]):
        return "M0", "missing pilot comparison; default M0"
    m0 = ok.loc[ok["model_id"] == "M0"].iloc[0]
    lam0 = float(m0["lambda_gc_nonsignificant"])
    or0 = float(m0["enrichment_or"])
    n0 = int(m0["n_significant_fdr"])

    best = None
    best_score = -1e9
    notes = []
    for _, r in ok.iterrows():
        mid = r["model_id"]
        if mid == "M0":
            continue
        lam = float(r["lambda_gc_nonsignificant"])
        or_ = float(r["enrichment_or"])
        n_sig = int(r["n_significant_fdr"])
        improve = (lam0 - lam) >= 0.15 or lam <= 1.20
        retain_ext = (or_ / or0) >= 0.80 if or0 > 0 else False
        not_collapse = n_sig >= 0.5 * n0
        if improve and retain_ext and not_collapse:
            # prefer larger λ improvement, then higher external OR
            score = (lam0 - lam) + 0.01 * (or_ / or0)
            notes.append(f"{mid}: Δλ_NS={lam0-lam:.3f}, OR_ratio={or_/or0:.3f}, n_sig={n_sig}")
            if score > best_score:
                best_score = score
                best = mid
        else:
            notes.append(
                f"{mid}: fail improve={improve} retain_ext={retain_ext} "
                f"not_collapse={not_collapse} (λ_NS={lam:.3f})"
            )
    if best is None:
        return "M0", "no model met λ_NS improvement + external retention; lock M0. " + "; ".join(notes)
    return best, f"selected {best} on pilot metrics. " + "; ".join(notes)


def main() -> None:
    args = parse_args()
    pilot = load_comparison(args.pilot_region)
    chosen, reason = candidate_from_pilot(pilot) if not pilot.empty else ("M0", "no pilot table")

    lines = [
        "# Phase 1 lock decision",
        "",
        f"**Date:** {date.today().isoformat()}",
        f"**Pilot region:** {args.pilot_region}",
        f"**Chosen primary covariate model:** `{chosen}`",
        "",
        "## Rationale",
        "",
        reason,
        "",
        "## Go/no-go rules applied",
        "",
        "- Lock new model if λ_NS improves by ≥0.15 or reaches ≤1.20, external enrichment OR ≥80% of M0, and n_sig ≥50% of M0.",
        "- Else lock M0 and document residual inflation as likely true-signal burden.",
        "",
        "## Per-region comparison snapshots",
        "",
    ]
    for region in args.regions:
        tab = load_comparison(region)
        lines.append(f"### {region}")
        lines.append("")
        if tab.empty:
            lines.append("_comparison_summary.tsv missing — TensorQTL sensitivity not finished._")
            lines.append("")
            continue
        keep = tab[tab.get("status", "ok") == "ok"] if "status" in tab.columns else tab
        cols = [c for c in [
            "model_id", "n_significant_fdr", "lambda_gc_nonsignificant",
            "enrichment_or", "enrichment_or_vs_M0",
        ] if c in keep.columns]
        lines.append("```")
        lines.append(keep[cols].to_string(index=False))
        lines.append("```")
        lines.append("")

    lines += [
        "## Next steps",
        "",
        f"1. If `{chosen}` ≠ M0: promote `covariates_{chosen}.txt` to primary `prepared/covariates.txt` after confirming DLPFC/hippocampus, then re-run primary TensorQTL only if replacing M0 outputs.",
        "2. Recompute Phase 2 burden if significance sets change materially.",
        "3. Update `config/covariates.yml` primary formula to match the locked model.",
        "4. Phase 2/3 external architecture tests already support the predictability axis; re-check only if calls change.",
        "",
        "## Status",
        "",
        (
            "**PROVISIONAL** — awaiting completed sensitivity TensorQTL for non-M0 models."
            if pilot.empty or (pilot["status"] == "missing_cis_qtl").all()
            else f"**DECISION RECORDED:** lock `{chosen}` pending cross-region confirmation."
        ),
        "",
    ]

    out = PHASE1 / "PHASE1_LOCK_DECISION.md"
    out.write_text("\n".join(lines))
    # also write machine-readable choice
    choice_path = PHASE1 / "_m" / "phase1_lock_choice.tsv"
    choice_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{
        "pilot_region": args.pilot_region,
        "chosen_model": chosen,
        "date": date.today().isoformat(),
        "reason": reason[:500],
    }]).to_csv(choice_path, sep="\t", index=False)
    print(f"Wrote {out}")
    print(f"Chose {chosen}")


if __name__ == "__main__":
    main()
