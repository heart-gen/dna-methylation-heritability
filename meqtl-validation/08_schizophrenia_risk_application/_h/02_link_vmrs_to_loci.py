#!/usr/bin/env python3
"""Analysis 2: link VMRs to SCZ-risk loci via meQTL / overlap / proximity.

Primary evidence: CpG-level meQTL connection (risk variant is genome-wide lead,
or risk variant lies within cis window of a significant meQTL CpG).
Proximity-only links are flagged exploratory.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
WINDOW = 500_000
FDR = 0.05


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", default="caudate")
    p.add_argument("--population", default="AA", choices=["AA", "EA"])
    p.add_argument("--window", type=int, default=WINDOW)
    p.add_argument("--fdr", type=float, default=FDR)
    p.add_argument("--outdir", default="")
    return p.parse_args()


def _rsid(x: str) -> str:
    m = re.search(r"(rs\d+)$", str(x))
    return m.group(1) if m else ""


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir) if args.outdir else (
        PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m" / args.region
    )
    outdir.mkdir(parents=True, exist_ok=True)

    index = pd.read_csv(outdir / "scz_index_snps_hg38.tsv", sep="\t")
    index = index[index["in_genotype_panel"].fillna(False).astype(bool)].copy()
    if index.empty:
        raise SystemExit("No index SNPs in genotype panel; run 01_define_scz_loci.py first")

    loci = pd.read_csv(outdir / "scz_risk_loci_hg38.tsv", sep="\t")
    phase1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping" / args.region / "_m"
    lead_path = phase1 / "tensorqtl" / "qc" / "lead_snp_per_cpg.tsv.gz"
    if args.population != "AA":
        lead_path = phase1 / "tensorqtl" / args.population / "qc" / "lead_snp_per_cpg.tsv.gz"
    lead = pd.read_csv(lead_path, sep="\t", compression="infer")
    if "phenotype_id" not in lead.columns:
        lead = lead.rename(columns={lead.columns[0]: "phenotype_id"})
    lead["sig_meqtl"] = lead["qval"].le(args.fdr) if "qval" in lead.columns else False
    lead["lead_rsid"] = lead["variant_id"].map(_rsid) if "variant_id" in lead.columns else ""

    prepared = phase1 / "prepared"
    maps = sorted(prepared.glob("cpg_vmr_map.chr*.tsv"))
    cpg_map = pd.concat([pd.read_csv(p, sep="\t") for p in maps], ignore_index=True)
    cpg_map["chrom"] = cpg_map["chrom"].astype(str)
    if not str(cpg_map["chrom"].iloc[0]).startswith("chr"):
        cpg_map["chrom"] = "chr" + cpg_map["chrom"].str.replace("^chr", "", regex=True)
    cpg_map["pos"] = pd.to_numeric(cpg_map["pos_1based"], errors="coerce")

    burden = pd.read_csv(
        PROJECT / "meqtl-validation/02_vmr_meqtl_burden/_m" / args.region / "vmr_meqtl_burden.tsv.gz",
        sep="\t",
    )
    burden["vmr_id"] = burden["vmr_id"].astype(str)

    # Parse VMR coordinates when available from interval-style ids
    def parse_vmr(vid: str) -> tuple[str, int, int] | tuple[None, None, None]:
        s = str(vid)
        m = re.match(r"^(?:chr)?(\d+|X|Y):(\d+)-(\d+)$", s)
        if m:
            return f"chr{m.group(1)}", int(m.group(2)), int(m.group(3))
        m = re.match(r"^chr([0-9XY]+)_(\d+)_(\d+)$", s)
        if m:
            return f"chr{m.group(1)}", int(m.group(2)), int(m.group(3))
        return None, None, None

    vmr_coords = []
    for vid in burden["vmr_id"].unique():
        chrom, start, end = parse_vmr(vid)
        if chrom is not None:
            vmr_coords.append({"vmr_id": str(vid), "vmr_chrom": chrom, "vmr_start": start, "vmr_end": end})
    # Prefer coordinates from cpg_map aggregation
    vmr_from_cpg = (
        cpg_map.groupby("vmr_id", as_index=False)
        .agg(vmr_chrom=("chrom", "first"), vmr_start=("pos", "min"), vmr_end=("pos", "max"))
    )
    vmr_from_cpg["vmr_id"] = vmr_from_cpg["vmr_id"].astype(str)
    vmr_ann = vmr_from_cpg.copy()
    if vmr_coords:
        extra = pd.DataFrame(vmr_coords)
        vmr_ann = pd.concat([vmr_ann, extra], ignore_index=True).drop_duplicates("vmr_id", keep="first")

    index["index_snp"] = index["index_snp"].astype(str)
    index["index_rsid"] = index["index_snp"].where(index["index_snp"].str.startswith("rs"), "")
    risk_rsids = set(index["index_rsid"].replace("", np.nan).dropna())
    risk_ids = set(index["genotype_variant_id"].dropna().astype(str))

    # --- Link type A: risk variant is genome-wide lead for a significant CpG ---
    lead_sig = lead[lead["sig_meqtl"]].copy()
    lead_sig["link_direct_lead"] = (
        lead_sig["variant_id"].astype(str).isin(risk_ids)
        | lead_sig["lead_rsid"].isin(risk_rsids)
    )
    direct = lead_sig[lead_sig["link_direct_lead"]].merge(
        cpg_map[["phenotype_id", "vmr_id", "chrom", "pos"]],
        on="phenotype_id",
        how="inner",
    )
    # attach index SNP
    direct = direct.merge(
        index[["index_snp", "genotype_variant_id", "locus_id", "chrom", "pos_hg38"]].rename(
            columns={"chrom": "risk_chrom", "pos_hg38": "risk_pos"}
        ),
        left_on="variant_id",
        right_on="genotype_variant_id",
        how="left",
    )
    # rsid fallback
    miss = direct["index_snp"].isna() & direct["lead_rsid"].ne("")
    if miss.any():
        by_rs = index[index["index_rsid"].ne("")].drop_duplicates("index_rsid")
        fill = direct.loc[miss, ["lead_rsid"]].merge(
            by_rs[["index_rsid", "index_snp", "genotype_variant_id", "locus_id", "chrom", "pos_hg38"]].rename(
                columns={"index_rsid": "lead_rsid", "chrom": "risk_chrom", "pos_hg38": "risk_pos"}
            ),
            on="lead_rsid",
            how="left",
        )
        for c in ["index_snp", "genotype_variant_id", "locus_id", "risk_chrom", "risk_pos"]:
            direct.loc[miss, c] = fill[c].to_numpy()

    direct["link_type"] = "cpg_meqtl_lead_is_risk_variant"
    direct["evidence_tier"] = "primary"

    # --- Link type B: significant meQTL CpG within cis window of risk variant ---
    rows_b = []
    for _, r in index.dropna(subset=["chrom", "pos_hg38"]).iterrows():
        chrom = str(r["chrom"])
        pos = int(r["pos_hg38"])
        cpg_near = cpg_map[
            (cpg_map["chrom"] == chrom)
            & (cpg_map["pos"] >= pos - args.window)
            & (cpg_map["pos"] <= pos + args.window)
        ]
        if cpg_near.empty:
            continue
        merged = cpg_near.merge(lead_sig[["phenotype_id", "variant_id", "qval", "slope", "pval_nominal"]], on="phenotype_id", how="inner")
        if merged.empty:
            continue
        merged = merged.copy()
        merged["index_snp"] = r["index_snp"]
        merged["genotype_variant_id"] = r["genotype_variant_id"]
        merged["locus_id"] = r["locus_id"]
        merged["risk_chrom"] = chrom
        merged["risk_pos"] = pos
        merged["cpg_to_risk_distance"] = (merged["pos"] - pos).abs()
        merged["link_type"] = "sig_meqtl_cpg_within_cis_window"
        merged["evidence_tier"] = "primary"
        rows_b.append(merged)
    window_links = pd.concat(rows_b, ignore_index=True) if rows_b else pd.DataFrame()

    # --- Link type C: VMR overlaps risk variant (direct overlap) ---
    overlap_rows = []
    for _, r in index.dropna(subset=["chrom", "pos_hg38"]).iterrows():
        hits = vmr_ann[
            (vmr_ann["vmr_chrom"] == r["chrom"])
            & (vmr_ann["vmr_start"] <= int(r["pos_hg38"]))
            & (vmr_ann["vmr_end"] >= int(r["pos_hg38"]))
        ]
        for _, v in hits.iterrows():
            overlap_rows.append({
                "vmr_id": v["vmr_id"],
                "index_snp": r["index_snp"],
                "genotype_variant_id": r["genotype_variant_id"],
                "locus_id": r["locus_id"],
                "risk_chrom": r["chrom"],
                "risk_pos": int(r["pos_hg38"]),
                "link_type": "vmr_overlaps_risk_variant",
                "evidence_tier": "primary",
            })
    overlap = pd.DataFrame(overlap_rows)

    # --- Link type D: proximity only (exploratory) ---
    prox_rows = []
    for _, r in index.dropna(subset=["chrom", "pos_hg38"]).iterrows():
        hits = vmr_ann[
            (vmr_ann["vmr_chrom"] == r["chrom"])
            & (vmr_ann["vmr_end"] >= int(r["pos_hg38"]) - args.window)
            & (vmr_ann["vmr_start"] <= int(r["pos_hg38"]) + args.window)
        ]
        for _, v in hits.iterrows():
            prox_rows.append({
                "vmr_id": v["vmr_id"],
                "index_snp": r["index_snp"],
                "genotype_variant_id": r["genotype_variant_id"],
                "locus_id": r["locus_id"],
                "risk_chrom": r["chrom"],
                "risk_pos": int(r["pos_hg38"]),
                "link_type": "proximity_within_cis_window",
                "evidence_tier": "exploratory",
            })
    proximity = pd.DataFrame(prox_rows)

    parts = []
    for df, cols_extra in [
        (direct, ["phenotype_id", "variant_id", "qval", "slope", "pval_nominal", "chrom", "pos"]),
        (window_links, ["phenotype_id", "variant_id", "qval", "slope", "pval_nominal", "chrom", "pos", "cpg_to_risk_distance"]),
        (overlap, []),
        (proximity, []),
    ]:
        if df is None or df.empty:
            continue
        base = ["vmr_id", "index_snp", "genotype_variant_id", "locus_id", "risk_chrom", "risk_pos", "link_type", "evidence_tier"]
        use = [c for c in base + cols_extra if c in df.columns]
        parts.append(df[use].copy())

    links = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if not links.empty:
        links["vmr_id"] = links["vmr_id"].astype(str)
        links = links.merge(
            burden[["vmr_id", "local_predictability", "proportion_cpgs_with_sig_meqtl", "n_tested_cpgs"]],
            on="vmr_id",
            how="left",
        )
        links = links.drop_duplicates(
            subset=[c for c in ["vmr_id", "index_snp", "phenotype_id", "link_type"] if c in links.columns]
        )

    links_path = outdir / "vmr_locus_links.tsv.gz"
    links.to_csv(links_path, sep="\t", index=False, compression="gzip")

    primary = links[links["evidence_tier"] == "primary"] if not links.empty else links
    write_tsv(outdir / "link_summary.tsv", [{
        "region": args.region,
        "population": args.population,
        "window_bp": args.window,
        "n_index_in_panel": int(len(index)),
        "n_link_rows": int(len(links)),
        "n_primary_link_rows": int(len(primary)),
        "n_vmrs_primary": int(primary["vmr_id"].nunique()) if len(primary) else 0,
        "n_loci_primary": int(primary["locus_id"].nunique()) if len(primary) and primary["locus_id"].notna().any() else 0,
        "n_direct_lead_links": int((links["link_type"] == "cpg_meqtl_lead_is_risk_variant").sum()) if len(links) else 0,
        "n_window_sig_meqtl_links": int((links["link_type"] == "sig_meqtl_cpg_within_cis_window").sum()) if len(links) else 0,
        "n_overlap_links": int((links["link_type"] == "vmr_overlaps_risk_variant").sum()) if len(links) else 0,
        "lead_path": str(lead_path),
        "output": str(links_path),
    }])
    print(f"Wrote {links_path}; primary VMRs={int(primary['vmr_id'].nunique()) if len(primary) else 0}")


if __name__ == "__main__":
    main()
