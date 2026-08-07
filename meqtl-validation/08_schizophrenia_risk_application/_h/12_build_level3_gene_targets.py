#!/usr/bin/env python3
"""Build Level-3 gene/PSI targets for prioritized SCZ loci.

Links prioritized VMRs → expression/PSI associations + FINEMAP gene symbols.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
REG = (
    ROOT
    / "heritability/elastic_net_model/all_individuals/tissue_comparison/regulatory_context/_m"
)
PRED = (
    ROOT
    / "heritability/elastic_net_model/all_individuals/caudate/_m/caudate_summary_elastic-net_AA.tsv"
)
PRIOR = (
    ROOT
    / "meqtl-validation/08_schizophrenia_risk_application/_m/prioritized/prioritized_loci.tsv"
)
OUTDIR = (
    ROOT / "meqtl-validation/08_schizophrenia_risk_application/_m/level3"
)


def strip_ensg(x: str) -> str:
    s = str(x)
    return s.split(".")[0] if s.startswith("ENSG") else s


def load_task_coord_map() -> pd.DataFrame:
    pred = pd.read_csv(PRED, sep="\t")
    pred["task_id"] = pred["task_id"].astype(str)
    chrom = pred["chrom"].astype(str).str.replace("^chr", "", regex=True)
    pred["vmr_coord"] = "chr" + chrom + "_" + pred["start"].astype(str) + "_" + pred["end"].astype(str)
    return pred[["task_id", "vmr_coord", "chrom", "start", "end"]]


def expand_locus_vmrs(prior: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in prior.iterrows():
        vmrs = set()
        for field in ("best_vmr_id", "tx_vmrs", "scored_vmr_id"):
            raw = str(r.get(field, "") or "")
            if raw and raw.lower() != "nan":
                vmrs |= {x.strip() for x in raw.split(",") if x.strip()}
        for vid in sorted(vmrs):
            rows.append(
                {
                    "locus_id": r["locus_id"],
                    "index_snp": r["index_snp"],
                    "task_id": str(vid),
                    "tx_expression": bool(r.get("tx_expression", False)),
                    "tx_psi": bool(r.get("tx_psi", False)),
                    "finemap_genes": r.get("finemap_genes", ""),
                    "gwas_or": r.get("gwas_or"),
                    "gwas_p": r.get("gwas_p"),
                }
            )
    return pd.DataFrame(rows)


def load_assoc(path: Path, modality: str, coords: set[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df[df["vmr_id"].astype(str).isin(coords)].copy()
    if df.empty:
        return df
    df["modality"] = modality
    df["gene_id"] = df["feature_id"].astype(str).map(strip_ensg)
    if "feature_label" in df.columns:
        df["gene_symbol"] = df["feature_label"].astype(str)
    else:
        df["gene_symbol"] = df.get("target_name", pd.Series(dtype=str)).astype(str)
    if modality == "psi":
        # feature_label like p34722|NMRAL1|A5
        parts = df["feature_label"].astype(str).str.split("|", expand=True)
        if parts.shape[1] >= 2:
            df["gene_symbol"] = parts[1]
            df["gene_id"] = df["target_id"].astype(str).map(strip_ensg)
    keep = [
        "vmr_id",
        "modality",
        "feature_id",
        "gene_id",
        "gene_symbol",
        "beta",
        "p",
        "fdr",
        "sig_fdr_05",
        "distance",
        "link_type",
    ]
    return df[[c for c in keep if c in df.columns]]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=str(OUTDIR))
    parser.add_argument("--fdr", type=float, default=0.05)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prior = pd.read_csv(PRIOR, sep="\t")
    locus_vmrs = expand_locus_vmrs(prior)
    coord_map = load_task_coord_map()
    locus_vmrs = locus_vmrs.merge(coord_map, on="task_id", how="left")
    coords = set(locus_vmrs["vmr_coord"].dropna().astype(str))

    expr_path = REG / "caudate/AA/expression/nearest_gene_window_250kb/vmr_feature_associations.tsv.gz"
    psi_path = REG / "caudate/AA/psi/window_250kb/vmr_feature_associations.tsv.gz"
    expr = load_assoc(expr_path, "expression", coords)
    psi = load_assoc(psi_path, "psi", coords)
    assoc = pd.concat([expr, psi], ignore_index=True)

    linked = locus_vmrs.merge(assoc, left_on="vmr_coord", right_on="vmr_id", how="left")

    # FINEMAP symbols as additional gene targets (no ENSG until annotation join)
    fine_rows = []
    for _, r in prior.iterrows():
        genes = str(r.get("finemap_genes", "") or "")
        if not genes or genes.lower() == "nan":
            continue
        for g in genes.split(","):
            g = g.strip()
            if not g or g == "-":
                continue
            fine_rows.append(
                {
                    "locus_id": r["locus_id"],
                    "index_snp": r["index_snp"],
                    "task_id": "",
                    "vmr_coord": "",
                    "modality": "finemap_symbol",
                    "feature_id": "",
                    "gene_id": "",
                    "gene_symbol": g,
                    "beta": pd.NA,
                    "p": pd.NA,
                    "fdr": pd.NA,
                    "sig_fdr_05": False,
                    "distance": pd.NA,
                    "link_type": "finemap",
                    "tx_expression": bool(r.get("tx_expression", False)),
                    "tx_psi": bool(r.get("tx_psi", False)),
                    "finemap_genes": r.get("finemap_genes", ""),
                    "gwas_or": r.get("gwas_or"),
                    "gwas_p": r.get("gwas_p"),
                }
            )
    if fine_rows:
        linked = pd.concat([linked, pd.DataFrame(fine_rows)], ignore_index=True)

    # Resolve FINEMAP symbols → ENSG via gene annotation
    annot = pd.read_csv(
        ROOT / "inputs/counts/gene-annotation.tsv",
        sep="\t",
        usecols=["gene_id", "gene_name", "chrom", "start", "end", "strand"],
    )
    annot["gene_id_nov"] = annot["gene_id"].map(strip_ensg)
    sym2ensg = (
        annot.dropna(subset=["gene_name"])
        .drop_duplicates("gene_name")
        .set_index("gene_name")["gene_id_nov"]
        .to_dict()
    )
    miss = linked["gene_id"].isna() | (linked["gene_id"].astype(str) == "") | (linked["gene_id"].astype(str) == "nan")
    linked.loc[miss, "gene_id"] = linked.loc[miss, "gene_symbol"].map(sym2ensg)

    linked["sig_tx_fdr"] = linked["fdr"].notna() & (linked["fdr"] <= args.fdr)
    linked["source_priority"] = linked["modality"].map(
        {"expression": 0, "psi": 1, "finemap_symbol": 2}
    ).fillna(9)

    # Deduplicate primary gene list per locus×gene
    gene_targets = (
        linked[linked["gene_id"].notna() & (linked["gene_id"].astype(str) != "")]
        .sort_values(["locus_id", "gene_id", "source_priority", "fdr"], na_position="last")
        .drop_duplicates(["locus_id", "gene_id", "modality"], keep="first")
    )

    linked.to_csv(outdir / "level3_vmr_feature_links.tsv", sep="\t", index=False)
    gene_targets.to_csv(outdir / "level3_gene_targets.tsv", sep="\t", index=False)

    # Compact unique ENSG list for eQTL subsetting / lookups
    uniq = (
        gene_targets[["gene_id", "gene_symbol"]]
        .drop_duplicates("gene_id")
        .sort_values("gene_id")
    )
    uniq.to_csv(outdir / "level3_unique_genes.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "n_prioritized_loci": prior.shape[0],
                "n_vmr_links": locus_vmrs.shape[0],
                "n_expression_assoc_rows": int((linked["modality"] == "expression").sum()),
                "n_psi_assoc_rows": int((linked["modality"] == "psi").sum()),
                "n_finemap_symbol_rows": int((linked["modality"] == "finemap_symbol").sum()),
                "n_unique_ensg": uniq.shape[0],
                "n_sig_expression": int(
                    ((linked["modality"] == "expression") & linked["sig_tx_fdr"]).sum()
                ),
                "n_sig_psi": int(((linked["modality"] == "psi") & linked["sig_tx_fdr"]).sum()),
            }
        ]
    )
    summary.to_csv(outdir / "level3_gene_targets_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"Wrote targets under {outdir}")


if __name__ == "__main__":
    main()
