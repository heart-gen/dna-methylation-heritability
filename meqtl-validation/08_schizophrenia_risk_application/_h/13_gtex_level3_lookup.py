#!/usr/bin/env python3
"""GTEx v11 brain eQTL/sQTL lookup for Level-3 prioritized gene targets × index SNPs."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

ROOT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
GTEX_EQTL = Path("/projects/b1213/resources/public_data/gtex_v11/GTEx_Analysis_v11_eQTL")
GTEX_SQTL = Path("/projects/b1213/resources/public_data/gtex_v11/GTEx_Analysis_v11_sQTL")
LEVEL3 = ROOT / "meqtl-validation/08_schizophrenia_risk_application/_m/level3"
PRIOR = (
    ROOT
    / "meqtl-validation/08_schizophrenia_risk_application/_m/prioritized/prioritized_loci.tsv"
)
INDEX = (
    ROOT
    / "meqtl-validation/08_schizophrenia_risk_application/_m/caudate/scz_index_snps_hg38.tsv"
)

TISSUES = {
    "caudate": "Brain_Caudate_basal_ganglia",
    "dlpfc": "Brain_Frontal_Cortex_BA9",
    "hippocampus": "Brain_Hippocampus",
}


def strip_ensg(x: str) -> str:
    s = str(x)
    return s.split(".")[0] if s.startswith("ENSG") else s


def variant_key_from_gtex(vid: str) -> str:
    # chr1_985133_G_A_b38
    m = re.match(r"^(chr[^_]+)_(\d+)_([ACGT]+)_([ACGT]+)", str(vid))
    if not m:
        return ""
    return f"{m.group(1)}_{m.group(2)}_{m.group(3)}_{m.group(4)}"


def variant_key_from_index(row: pd.Series) -> str:
    chrom = str(row["chrom"])
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return f"{chrom}_{int(row['pos_hg38'])}_{row['REF']}_{row['ALT']}"


def load_signif(path: Path, genes: set[str], var_keys: set[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_parquet(path)
    df["gene_id"] = df["phenotype_id"].astype(str).map(strip_ensg)
    # sQTL phenotype_ids are often intron/cluster; keep gene from eGenes-style if present
    if "gene_id" not in df.columns or df["gene_id"].eq(df["phenotype_id"].map(strip_ensg)).all():
        pass
    df["variant_key"] = df["variant_id"].map(variant_key_from_gtex)
    hit = df[df["gene_id"].isin(genes) | df["variant_key"].isin(var_keys)].copy()
    # Prefer gene∩variant, else gene-level support, else variant-level
    hit["match_type"] = "other"
    both = hit["gene_id"].isin(genes) & hit["variant_key"].isin(var_keys)
    hit.loc[both, "match_type"] = "gene_and_variant"
    hit.loc[~both & hit["gene_id"].isin(genes), "match_type"] = "gene_any_variant"
    hit.loc[~both & hit["variant_key"].isin(var_keys), "match_type"] = "variant_any_gene"
    return hit


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default=str(LEVEL3 / "gtex"))
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    targets = pd.read_csv(LEVEL3 / "level3_gene_targets.tsv", sep="\t")
    genes = set(targets["gene_id"].dropna().astype(str).map(strip_ensg))
    prior = pd.read_csv(PRIOR, sep="\t")
    idx = pd.read_csv(INDEX, sep="\t")
    idx = idx[idx["index_snp"].isin(prior["index_snp"])].copy()
    idx["variant_key"] = idx.apply(variant_key_from_index, axis=1)
    var_keys = set(idx["variant_key"])
    snp_by_key = idx.drop_duplicates("variant_key").set_index("variant_key")["index_snp"].to_dict()

    rows = []
    for tissue_key, tissue_name in TISSUES.items():
        for qtl, base in (("eQTL", GTEX_EQTL), ("sQTL", GTEX_SQTL)):
            path = base / f"{tissue_name}.v11.{qtl}s.signif_pairs.parquet"
            # GTEx naming: eQTLs / sQTLs
            if not path.exists():
                path = base / f"{tissue_name}.v11.{qtl}.signif_pairs.parquet"
            hit = load_signif(path, genes, var_keys)
            if hit.empty:
                continue
            hit["tissue"] = tissue_key
            hit["gtex_tissue"] = tissue_name
            hit["qtl_type"] = qtl
            hit["index_snp"] = hit["variant_key"].map(snp_by_key)
            # Attach locus/gene symbols of interest
            gene_map = (
                targets.dropna(subset=["gene_id"])
                .assign(gene_id=lambda d: d["gene_id"].map(strip_ensg))
                .drop_duplicates("gene_id")
                .set_index("gene_id")["gene_symbol"]
                .to_dict()
            )
            hit["gene_symbol"] = hit["gene_id"].map(gene_map)
            keep = [
                "tissue",
                "gtex_tissue",
                "qtl_type",
                "match_type",
                "phenotype_id",
                "gene_id",
                "gene_symbol",
                "variant_id",
                "variant_key",
                "index_snp",
                "pval_nominal",
                "slope",
                "slope_se",
                "af",
            ]
            rows.append(hit[[c for c in keep if c in hit.columns]])

    all_hits = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    all_hits.to_csv(outdir / "gtex_level3_signif_hits.tsv.gz", sep="\t", index=False)

    # Locus-level rollup: any gene×variant match in caudate eQTL preferred
    locus_gene = targets[["locus_id", "index_snp", "gene_id", "gene_symbol", "modality"]].copy()
    locus_gene["gene_id"] = locus_gene["gene_id"].map(strip_ensg)
    if all_hits.empty:
        summary = locus_gene.assign(
            gtex_caudate_eqtl_gene_variant=False,
            gtex_any_eqtl_gene=False,
            gtex_any_sqtl_gene=False,
        )
    else:
        caud_gv = set(
            zip(
                all_hits.loc[
                    (all_hits["tissue"] == "caudate")
                    & (all_hits["qtl_type"] == "eQTL")
                    & (all_hits["match_type"] == "gene_and_variant"),
                    "gene_id",
                ],
                all_hits.loc[
                    (all_hits["tissue"] == "caudate")
                    & (all_hits["qtl_type"] == "eQTL")
                    & (all_hits["match_type"] == "gene_and_variant"),
                    "index_snp",
                ],
            )
        )
        any_eqtl_gene = set(
            all_hits.loc[all_hits["qtl_type"] == "eQTL", "gene_id"].astype(str)
        )
        any_sqtl_gene = set(
            all_hits.loc[all_hits["qtl_type"] == "sQTL", "gene_id"].astype(str)
        )
        summary = locus_gene.copy()
        summary["gtex_caudate_eqtl_gene_variant"] = [
            (g, s) in caud_gv
            for g, s in zip(summary["gene_id"].astype(str), summary["index_snp"].astype(str))
        ]
        summary["gtex_any_eqtl_gene"] = summary["gene_id"].astype(str).isin(any_eqtl_gene)
        summary["gtex_any_sqtl_gene"] = summary["gene_id"].astype(str).isin(any_sqtl_gene)

    summary.to_csv(outdir / "gtex_level3_locus_gene_summary.tsv", sep="\t", index=False)

    claim = pd.DataFrame(
        [
            {
                "n_signif_pair_hits": int(len(all_hits)),
                "n_gene_and_variant_hits": int(
                    (all_hits["match_type"] == "gene_and_variant").sum()
                )
                if len(all_hits)
                else 0,
                "n_loci_with_caudate_eqtl_gene_variant": int(
                    summary.loc[summary["gtex_caudate_eqtl_gene_variant"], "locus_id"].nunique()
                )
                if len(summary)
                else 0,
                "n_genes_with_any_brain_eqtl": int(summary["gtex_any_eqtl_gene"].sum())
                if len(summary)
                else 0,
            }
        ]
    )
    claim.to_csv(outdir / "gtex_level3_claim_snapshot.tsv", sep="\t", index=False)
    print(claim.to_string(index=False))
    print(f"Wrote {outdir}")


if __name__ == "__main__":
    main()
