#!/usr/bin/env python3
"""Analysis 6 (light): transcriptional coupling of SCZ meQTL-supported VMRs.

Reuses Phase 5 / regulatory_context association tables. No colocalization in v1.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
REGCTX = (
    PROJECT / "heritability/elastic_net_model/all_individuals/"
    "tissue_comparison/regulatory_context/_m"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", default="caudate")
    p.add_argument("--outdir", default="")
    return p.parse_args()


def arch_vmr_to_coord(vmr_id: str) -> str:
    parts = str(vmr_id).split("_")
    if len(parts) >= 3 and parts[0].startswith("chr"):
        return f"{parts[0].replace('chr', '')}:{parts[1]}-{parts[2]}"
    if ":" in str(vmr_id) and "-" in str(vmr_id):
        return str(vmr_id).replace("chr", "")
    return str(vmr_id)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir) if args.outdir else (
        PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m" / args.region
    )
    vmr_path = outdir / "scz_vmr_level_summary.tsv.gz"
    if not vmr_path.exists():
        raise SystemExit(f"Missing {vmr_path}")
    vmr = pd.read_csv(vmr_path, sep="\t")
    vmr["vmr_id"] = vmr["vmr_id"].astype(str)
    vmr["any_sig"] = vmr["any_sig"].astype(bool)

    rows = []
    exemplars = []
    for modality, rel in [
        ("expression", f"{args.region}/AA/expression/nearest_gene_window_250kb/architecture_model_input.tsv"),
        ("psi", f"{args.region}/AA/psi/window_250kb/architecture_model_input.tsv"),
    ]:
        path = REGCTX / rel
        if not path.exists():
            rows.append({"modality": modality, "error": f"missing:{path}"})
            continue
        tx = pd.read_csv(path, sep="\t")
        tx["coord_id"] = tx["vmr_id"].map(arch_vmr_to_coord)
        tx["tx_associated"] = pd.to_numeric(tx["any_sig_fdr_05"], errors="coerce").fillna(0).eq(1)

        # Map SCZ vmr_id → coord
        pred = PROJECT / "heritability/elastic_net_model/all_individuals" / args.region / "_m" / f"{args.region}_summary_elastic-net_AA.tsv"
        p = pd.read_csv(pred, sep="\t")
        p["task_id"] = p["task_id"].astype(str)
        p["coord_id"] = (
            p["chrom"].astype(str).str.replace("^chr", "", regex=True)
            + ":" + p["start"].astype(str) + "-" + p["end"].astype(str)
        )
        vmr2 = vmr.merge(p[["task_id", "coord_id"]], left_on="vmr_id", right_on="task_id", how="left")
        # also try vmr_id already as coord
        miss = vmr2["coord_id"].isna()
        vmr2.loc[miss, "coord_id"] = vmr2.loc[miss, "vmr_id"].map(arch_vmr_to_coord)

        merged = vmr2.merge(
            tx[["coord_id", "tx_associated"]].drop_duplicates("coord_id"),
            on="coord_id",
            how="left",
        )
        merged["tx_associated"] = merged["tx_associated"].fillna(False)
        a = int(((merged["any_sig"]) & (merged["tx_associated"])).sum())
        b = int(((merged["any_sig"]) & (~merged["tx_associated"])).sum())
        c = int(((~merged["any_sig"]) & (merged["tx_associated"])).sum())
        d = int(((~merged["any_sig"]) & (~merged["tx_associated"])).sum())
        oddsratio, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "modality": modality,
            "n_scz_tested_vmrs": len(merged),
            "n_scz_sig_vmrs": int(merged["any_sig"].sum()),
            "n_scz_sig_and_tx": a,
            "odds_ratio": float(oddsratio),
            "fisher_pvalue": float(pval),
            "path": str(path),
        })
        hit = merged[merged["any_sig"] & merged["tx_associated"]].copy()
        hit["modality"] = modality
        exemplars.append(hit)

    write_tsv(outdir / "tx_integration_results.tsv", rows)
    if exemplars:
        ex = pd.concat(exemplars, ignore_index=True)
        ex.to_csv(outdir / "scz_meqtl_tx_coupled_vmrs.tsv.gz", sep="\t", index=False, compression="gzip")
    n_tx = int(sum(r.get("n_scz_sig_and_tx", 0) for r in rows if "n_scz_sig_and_tx" in r))
    write_tsv(outdir / "tx_decision_snapshot.tsv", [{
        "region": args.region,
        "n_sig_vmr_with_tx_coupling": n_tx,
        "criterion_ge1_locus_with_tx": n_tx >= 1,
    }])
    print(f"TX integration written under {outdir}; coupled VMRs={n_tx}")


if __name__ == "__main__":
    main()
