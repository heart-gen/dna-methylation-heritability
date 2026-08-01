#!/usr/bin/env python3
"""Compare covariate-sensitivity TensorQTL models (λ_NS, n_sig, external overlap)."""

from __future__ import annotations

import argparse, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
PHASE1 = PROJECT / "meqtl-validation" / "01_cpg_meqtl_mapping"
PHASE3 = PROJECT / "meqtl-validation" / "03_external_meqtl_validation" / "_m" / "harmonized"

PREFERRED_EXTERNAL = {
    "caudate": "jaffe_dlpfc_450k_meqtl",  # secondary tissue match
    "dlpfc": "jaffe_dlpfc_450k_meqtl",
    "hippocampus": "schulz_hippocampus_array_meqtl",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True)
    p.add_argument("--models", nargs="+", default=["M0", "M1", "M2", "M3a", "M3b", "M3c", "M4", "M5"])
    p.add_argument("--fdr", type=float, default=0.05)
    return p.parse_args()


def genomic_inflation(pvals: np.ndarray) -> float:
    from scipy.stats import chi2

    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if p.size == 0:
        return float("nan")
    return float(np.median(chi2.ppf(1 - p, 1)) / chi2.ppf(0.5, 1))


def load_external_supported_cpgs(resource: str, region: str) -> set[str]:
    path = PHASE3 / f"{resource}.{region}.vmr_support.tsv.gz"
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t", usecols=["phenotype_id", "external_meqtl_support"])
    return set(df.loc[df["external_meqtl_support"].astype(int) == 1, "phenotype_id"].astype(str))


def enrichment(sig_ids: set[str], all_ids: set[str], ext_ids: set[str]) -> dict:
    # 2x2 among tested CpGs that are in the annotated set (all_ids)
    ext = ext_ids & all_ids
    a = len(sig_ids & ext)
    b = len(sig_ids - ext)
    c = len(ext - sig_ids)
    d = len(all_ids - sig_ids - ext)
    # OR with Haldane-Anscombe correction
    or_ = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return {
        "n_sig": len(sig_ids),
        "n_ext_in_tested": len(ext),
        "n_sig_and_ext": a,
        "enrichment_or": float(or_),
        "frac_sig_with_ext": a / len(sig_ids) if sig_ids else np.nan,
    }


def resolve_cis_path(region: str, model: str) -> Path | None:
    sens = PHASE1 / region / "_m" / "covariate_sensitivity" / "tensorqtl" / model
    cand = sens / f"cpg_meqtl_{region}_{model}.cis_qtl.txt.gz"
    if cand.exists():
        return cand
    # M0 may reuse primary outputs
    if model == "M0":
        primary = PHASE1 / region / "_m" / "tensorqtl" / f"cpg_meqtl_{region}.cis_qtl.txt.gz"
        if primary.exists():
            return primary
    return None


def main() -> None:
    args = parse_args()
    region = args.region
    outdir = PHASE1 / region / "_m" / "covariate_sensitivity"
    outdir.mkdir(parents=True, exist_ok=True)
    resource = PREFERRED_EXTERNAL[region]
    ext_ids = load_external_supported_cpgs(resource, region)

    rows = []
    m0_or = None
    for model in args.models:
        cis = resolve_cis_path(region, model)
        if cis is None:
            rows.append({
                "region": region, "model_id": model, "status": "missing_cis_qtl",
                "external_resource": resource,
            })
            continue
        df = pd.read_csv(cis, sep="\t", index_col=0)
        if "phenotype_id" not in df.columns:
            df = df.reset_index().rename(columns={df.index.name or "index": "phenotype_id"})
            if "phenotype_id" not in df.columns:
                df = df.rename(columns={df.columns[0]: "phenotype_id"})
        # when index was phenotype id
        if "qval" not in df.columns and cis.name.endswith(".cis_qtl.txt.gz"):
            df = pd.read_csv(cis, sep="\t")
            if df.columns[0] not in {"phenotype_id", "pheno_id"}:
                # tensorqtl writes index unnamed sometimes
                pass
            df = pd.read_csv(cis, sep="\t", index_col=0)
            pheno = df.index.astype(str)
            q = pd.to_numeric(df["qval"], errors="coerce")
            p = pd.to_numeric(df["pval_beta"] if "pval_beta" in df.columns else df["pval_perm"], errors="coerce")
        else:
            # re-read cleanly
            df = pd.read_csv(cis, sep="\t", index_col=0)
            pheno = df.index.astype(str)
            q = pd.to_numeric(df["qval"], errors="coerce")
            pcol = "pval_beta" if "pval_beta" in df.columns else "pval_perm"
            p = pd.to_numeric(df[pcol], errors="coerce")

        sig_mask = q <= args.fdr
        ns_mask = q > args.fdr
        all_ids = set(pheno)
        sig_ids = set(pheno[sig_mask.fillna(False)])
        enr = enrichment(sig_ids, all_ids, ext_ids)
        lam_all = genomic_inflation(p.to_numpy())
        lam_ns = genomic_inflation(p[ns_mask.fillna(False)].to_numpy())
        if model == "M0":
            m0_or = enr["enrichment_or"]
        row = {
            "region": region,
            "model_id": model,
            "status": "ok",
            "cis_qtl_path": str(cis),
            "n_tested": len(df),
            "n_significant_fdr": int(sig_mask.fillna(False).sum()),
            "lambda_gc_all": lam_all,
            "lambda_gc_nonsignificant": lam_ns,
            "external_resource": resource,
            **enr,
            "enrichment_or_vs_M0": (
                float(enr["enrichment_or"] / m0_or) if (m0_or and model != "M0") else (1.0 if model == "M0" else np.nan)
            ),
        }
        rows.append(row)

    # Ensure M0 first for OR ratio — recompute ratios after
    tab = pd.DataFrame(rows)
    if (tab["model_id"] == "M0").any() and "enrichment_or" in tab.columns:
        m0_or = float(tab.loc[tab["model_id"] == "M0", "enrichment_or"].iloc[0])
        tab["enrichment_or_vs_M0"] = tab["enrichment_or"] / m0_or
    out = outdir / "comparison_summary.tsv"
    tab.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")
    if not tab.empty and "lambda_gc_nonsignificant" in tab.columns:
        print(tab[["model_id", "n_significant_fdr", "lambda_gc_nonsignificant", "enrichment_or"]].to_string(index=False))


if __name__ == "__main__":
    main()
