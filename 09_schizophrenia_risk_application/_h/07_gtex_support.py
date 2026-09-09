#!/usr/bin/env python3
"""External shared-genetic support from GTEx v11 brain (09).

Usage:
    python _h/07_gtex_support.py --run-id scz-AA-caudate-YYYYMMDD

AGENTS.md 7.8: "use external GTEx eQTL evidence as SUPPORT, not proof of
mediation", and the retention criteria require external shared genetic support
for at least one prioritized locus. Support here means one specific, checkable
thing: a prespecified schizophrenia risk variant (index SNP or LD proxy, from
stage 03) is itself a significant eQTL or sQTL for some gene in a GTEx brain
tissue. That is a shared-variant observation. It is not colocalization -- stages
08 to 10 do that properly, with priors, on full regional statistics -- and it is
not evidence of a causal chain through methylation.

Variants are matched on hg38 position and unordered allele pair. GTEx keys are
chr_pos_ref_alt_b38; the analysis panel keys are chr_pos_ref_alt_rsID. Matching
on position alone would join multi-allelic sites to the wrong allele, so the
allele pair is required and REF/ALT order is ignored (an allele swap is the same
variant with a sign flip, which does not matter for presence/absence support).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import pyarrow.dataset as ds
import pyarrow.compute as pc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_config import load_run_config  # noqa: E402


def repo_root() -> Path:
    d = Path(__file__).resolve()
    while d != d.parent:
        if (d / ".git").is_dir():
            return d
        d = d.parent
    raise SystemExit("Could not locate repository root")


def read_manifest(run_dir: Path) -> dict:
    m = pd.read_csv(run_dir / "manifest.tsv", sep="\t", dtype=str)
    return dict(zip(m["field"], m["value"]))


def variant_key(variant_id: str) -> str:
    """chrom_pos_{sorted allele pair} -- build- and order-independent."""
    f = str(variant_id).split("_")
    if len(f) < 4:
        return ""
    chrom, pos, ref, alt = f[0], f[1], f[2].upper(), f[3].upper()
    a, b = sorted((ref, alt))
    return f"{chrom}_{pos}_{a}_{b}"


def collect(allpairs_dir: Path, tissues: list[str], suffix: str,
            wanted: set[str], modality: str) -> pd.DataFrame:
    """Significant GTEx pairs whose variant is one of the risk variants."""
    rows = []
    for tissue in tissues:
        f = allpairs_dir / f"{tissue}.v11.{suffix}.signif_pairs.parquet"
        if not f.is_file():
            print(f"  {tissue}: no {suffix} signif_pairs file; skipped",
                  file=sys.stderr)
            continue
        tbl = ds.dataset(f, format="parquet").to_table(
            columns=["phenotype_id", "variant_id", "pval_nominal", "slope",
                     "slope_se", "pval_beta"])
        df = tbl.to_pandas()
        df["variant_key"] = df["variant_id"].map(variant_key)
        df = df[df["variant_key"].isin(wanted)]
        if df.empty:
            continue
        df["gtex_tissue"] = tissue
        df["modality"] = modality
        rows.append(df)
    if not rows:
        return pd.DataFrame(columns=[
            "phenotype_id", "variant_id", "pval_nominal", "slope", "slope_se",
            "pval_beta", "variant_key", "gtex_tissue", "modality"])
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "09_schizophrenia_risk_application" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")
    man = read_manifest(run_dir)
    cfg = load_run_config("schizophrenia", run_dir, root)
    gtex = cfg["colocalization"]["gtex"]
    tissues = list(gtex["tissues"])

    tests_f = run_dir / "results" / "risk-variant-cpg-tests.tsv.gz"
    if not tests_f.is_file():
        raise SystemExit(
            f"Missing {tests_f}; run 03b_combine_risk_variant_tests.R first.")
    tests = pd.read_csv(tests_f, sep="\t")
    if tests.empty:
        raise SystemExit("No risk-variant tests to look up in GTEx.")

    risk = tests[["locus_id", "index_snp", "risk_variant_id", "ld_r2"]].copy()
    risk = risk.drop_duplicates(["locus_id", "risk_variant_id"])
    risk["variant_key"] = risk["risk_variant_id"].map(variant_key)
    risk = risk[risk["variant_key"].ne("")]
    wanted = set(risk["variant_key"])
    print(f"Looking up {len(wanted)} risk variants across {len(tissues)} tissues")

    eqtl = collect(Path(gtex["eqtl_dir"]), tissues, "eQTLs", wanted, "eqtl")
    sqtl = collect(Path(gtex["sqtl_dir"]), tissues, "sQTLs", wanted, "sqtl")
    hits = pd.concat([eqtl, sqtl], ignore_index=True)

    if hits.empty:
        pairs = pd.DataFrame(columns=[
            "locus_id", "index_snp", "risk_variant_id", "ld_r2", "variant_key",
            "phenotype_id", "gtex_tissue", "modality", "pval_nominal", "slope",
            "slope_se", "pval_beta"])
    else:
        pairs = risk.merge(hits, on="variant_key", how="inner")
    pairs.to_csv(run_dir / "results" / "gtex-support-pairs.tsv.gz",
                 sep="\t", index=False, compression="gzip")

    # Locus-level rollup: the unit the retention criterion is stated in.
    loci = pd.read_csv(run_dir / "results" / "scz-loci.tsv", sep="\t")
    if pairs.empty:
        summary = pd.DataFrame({"locus_id": loci["locus_id"]})
        summary["n_risk_variants_with_gtex_support"] = 0
        summary["n_gtex_genes"] = 0
        summary["n_gtex_tissues"] = 0
        summary["gtex_modalities"] = ""
        summary["gtex_tissues"] = ""
        summary["gtex_genes"] = ""
    else:
        g = pairs.groupby("locus_id")
        summary = pd.DataFrame({
            "n_risk_variants_with_gtex_support": g["risk_variant_id"].nunique(),
            "n_gtex_genes": g["phenotype_id"].nunique(),
            "n_gtex_tissues": g["gtex_tissue"].nunique(),
            "gtex_modalities": g["modality"].apply(
                lambda s: ",".join(sorted(set(s)))),
            "gtex_tissues": g["gtex_tissue"].apply(
                lambda s: ",".join(sorted(set(s)))),
            # Capped: a locus can carry hundreds of genes and the full list
            # makes the table unreadable without adding information. The
            # complete set is in gtex-support-pairs.tsv.gz.
            "gtex_genes": g["phenotype_id"].apply(
                lambda s: ",".join(sorted(set(s))[:20])),
        }).reset_index()
        summary = loci[["locus_id"]].merge(summary, on="locus_id", how="left")
        for c, fill in [("n_risk_variants_with_gtex_support", 0),
                        ("n_gtex_genes", 0), ("n_gtex_tissues", 0),
                        ("gtex_modalities", ""), ("gtex_tissues", ""),
                        ("gtex_genes", "")]:
            summary[c] = summary[c].fillna(fill)

    summary["has_external_genetic_support"] = (
        summary["n_risk_variants_with_gtex_support"] > 0)
    summary["run_id"] = args.run_id
    summary["region"] = man["region"]
    summary["evidence_type"] = "shared_variant_presence"
    summary["permitted_claim"] = (
        "a schizophrenia risk variant at this locus is also a significant GTEx "
        "brain eQTL or sQTL")
    summary["forbidden_claim"] = (
        "methylation mediates the risk variant's effect on gene expression")
    summary.to_csv(run_dir / "results" / "gtex-support.tsv", sep="\t", index=False)

    n = int(summary["has_external_genetic_support"].sum())
    print(f"{n}/{len(summary)} loci carry external GTEx support")


if __name__ == "__main__":
    main()
