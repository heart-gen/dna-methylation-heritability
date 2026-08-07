#!/usr/bin/env python3
"""Prepare tidy tables for prioritized SCZ locus multi-panel figures."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
MODULE = PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m"
OUT = MODULE / "locus_panels"
PGC3 = Path(
    "/projects/b1213/resources/gwas/PGC/SCZ/PGC3/"
    "PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz"
)
WINDOW_PAD = 250_000  # bp around index SNP for GWAS track (hg19)


def read_pgc3_window(chrom_hg19: str, start: int, end: int) -> pd.DataFrame:
    chrom = str(chrom_hg19).replace("chr", "")
    rows = []
    with gzip.open(PGC3, "rt") as handle:
        header = None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("CHROM") or line.startswith("#CHROM"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            rec = dict(zip(header, parts))
            if str(rec["CHROM"]).replace("chr", "") != chrom:
                continue
            pos = int(float(rec["POS"]))
            if pos < start or pos > end:
                continue
            rows.append(
                {
                    "chrom_hg19": chrom,
                    "pos_hg19": pos,
                    "snp": rec["ID"],
                    "a1": rec["A1"],
                    "a2": rec["A2"],
                    "beta": float(rec["BETA"]) if rec["BETA"] not in ("", "NA") else float("nan"),
                    "se": float(rec["SE"]) if rec["SE"] not in ("", "NA") else float("nan"),
                    "pval": float(rec["PVAL"]) if rec["PVAL"] not in ("", "NA") else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pri = pd.read_csv(MODULE / "prioritized/prioritized_loci.tsv", sep="\t")
    idx = pd.read_csv(MODULE / "caudate/scz_index_snps_hg38.tsv", sep="\t")
    idx = idx.drop_duplicates("index_snp")

    # meQTL all regions
    meq = {}
    for region in ["caudate", "dlpfc", "hippocampus"]:
        p = MODULE / region / "risk_variant_cpg_meqtl.tsv.gz"
        meq[region] = pd.read_csv(p, sep="\t")

    gtex = pd.read_csv(MODULE / "level3/gtex/gtex_level3_signif_hits.tsv.gz", sep="\t")
    links = pd.read_csv(MODULE / "level3/level3_vmr_feature_links.tsv", sep="\t")
    l3 = pd.read_csv(MODULE / "level3/level3_locus_summary.tsv", sep="\t")
    tech = pd.read_csv(
        PROJECT
        / "meqtl-validation/07_repeat_mappability_sensitivity/_m/caudate/vmr_technical_annotations.tsv",
        sep="\t",
    )
    if "task_id" in tech.columns:
        tech["task_id"] = tech["task_id"].dropna().astype(float).astype(int).astype(str)
    en = pd.read_csv(
        PROJECT
        / "heritability/elastic_net_model/all_individuals/caudate/_m/caudate_summary_elastic-net_AA.tsv",
        sep="\t",
    )
    en["task_id"] = en["task_id"].astype(str)

    manifest = []
    for _, row in pri.iterrows():
        snp = str(row["index_snp"])
        locus_id = float(row["locus_id"])
        loc_dir = OUT / snp
        loc_dir.mkdir(parents=True, exist_ok=True)

        meta = idx[idx["index_snp"] == snp]
        if meta.empty:
            raise SystemExit(f"Missing index SNP metadata for {snp}")
        m = meta.iloc[0]
        chrom19 = str(int(float(m["chr_hg19"]))) if pd.notna(m["chr_hg19"]) else str(m["chrom"]).replace("chr", "")
        pos19 = int(float(m["bp_hg19"]))
        # prefer FINEMAP/locus range if present
        left = int(float(m["range_left_hg19"])) if pd.notna(m.get("range_left_hg19")) else pos19 - WINDOW_PAD
        right = int(float(m["range_right_hg19"])) if pd.notna(m.get("range_right_hg19")) else pos19 + WINDOW_PAD
        left = max(1, left)
        gwas = read_pgc3_window(chrom19, left, right)
        gwas.to_csv(loc_dir / "gwas_regional_hg19.tsv.gz", sep="\t", index=False, compression="gzip")

        # caudate meQTL CpG track
        c_meq = meq["caudate"]
        c_meq = c_meq[c_meq["index_snp"].astype(str) == snp].copy()
        c_meq.to_csv(loc_dir / "meqtl_caudate_cpgs.tsv.gz", sep="\t", index=False, compression="gzip")

        # cross-region forest: FDR pairs in any region for this index SNP, prefer shared
        forest_rows = []
        for region, df in meq.items():
            sub = df[(df["index_snp"].astype(str) == snp) & (df["significant_fdr"].astype(bool))].copy()
            if sub.empty:
                sub = df[df["index_snp"].astype(str) == snp].nsmallest(5, "pval_nominal")
            for _, r in sub.iterrows():
                forest_rows.append(
                    {
                        "region": region,
                        "index_snp": snp,
                        "phenotype_id": r["phenotype_id"],
                        "vmr_id": r["vmr_id"],
                        "beta": r["beta"],
                        "se": r["se"],
                        "pval_nominal": r["pval_nominal"],
                        "qval": r.get("qval", float("nan")),
                        "significant_fdr": bool(r.get("significant_fdr", False)),
                        "local_predictability": r.get("local_predictability", float("nan")),
                        "cpg_pos": r.get("cpg_pos", float("nan")),
                    }
                )
        forest = pd.DataFrame(forest_rows)
        # keep top CpGs by min p across regions for readability
        if not forest.empty:
            forest["abs_z"] = forest["beta"].abs() / forest["se"].clip(lower=1e-12)
            top_cpgs = (
                forest.groupby("phenotype_id")["abs_z"].max().sort_values(ascending=False).head(8).index
            )
            forest = forest[forest["phenotype_id"].isin(top_cpgs)]
        forest.to_csv(loc_dir / "cross_region_forest.tsv.gz", sep="\t", index=False, compression="gzip")

        # VMR predictability / tech for tx_vmrs
        vmrs = [str(int(float(row["best_vmr_id"])))]
        if pd.notna(row.get("tx_vmrs")) and str(row["tx_vmrs"]).strip():
            vmrs.extend([v.strip() for v in str(row["tx_vmrs"]).split(",") if v.strip()])
        vmrs = sorted(set(vmrs))
        vrows = []
        for vid in vmrs:
            er = en[en["task_id"] == vid]
            tr = tech[tech["task_id"].astype(str) == vid] if "task_id" in tech.columns else pd.DataFrame()
            vrows.append(
                {
                    "vmr_id": vid,
                    "is_best": vid == str(int(float(row["best_vmr_id"]))),
                    "local_predictability": float(er.iloc[0]["h2_unscaled"]) if len(er) else float("nan"),
                    "chrom": str(er.iloc[0]["chrom"]) if len(er) else "",
                    "start": int(er.iloc[0]["start"]) if len(er) else None,
                    "end": int(er.iloc[0]["end"]) if len(er) else None,
                    "umap_k24_mean": float(tr.iloc[0]["umap_k24_mean"]) if len(tr) and "umap_k24_mean" in tr else float("nan"),
                    "line_l1_frac": float(tr.iloc[0]["line_l1_frac"]) if len(tr) and "line_l1_frac" in tr else float("nan"),
                }
            )
        pd.DataFrame(vrows).to_csv(loc_dir / "vmr_predictability.tsv", sep="\t", index=False)

        # TX + Level3 annotations
        tx = links[(links["locus_id"] == locus_id) & (links["sig_tx_fdr"].astype(bool))].copy()
        tx.to_csv(loc_dir / "tx_links_fdr.tsv", sep="\t", index=False)
        gt = gtex[gtex["index_snp"].astype(str) == snp].copy() if "index_snp" in gtex.columns else gtex
        if "index_snp" not in gtex.columns and "variant_id" in gtex.columns:
            # already filtered in signif hits file sometimes by locus join upstream
            gt = gtex.copy()
        # filter gtex hits to this snp when column exists
        if "index_snp" in gtex.columns:
            gt = gtex[gtex["index_snp"].astype(str) == snp]
        gt.to_csv(loc_dir / "gtex_level3_hits.tsv", sep="\t", index=False)

        l3row = l3[l3["locus_id"] == locus_id]
        write_tsv(
            loc_dir / "locus_meta.tsv",
            [
                {
                    "rank": int(row["rank"]),
                    "locus_id": locus_id,
                    "index_snp": snp,
                    "chrom_hg38": m["chrom"],
                    "pos_hg38": int(m["pos_hg38"]),
                    "chrom_hg19": chrom19,
                    "pos_hg19": pos19,
                    "gwas_window_start_hg19": left,
                    "gwas_window_end_hg19": right,
                    "gwas_p": m["gwas_p"],
                    "gwas_or": m["gwas_or"],
                    "risk_allele": m["risk_allele"],
                    "best_vmr_id": int(float(row["best_vmr_id"])),
                    "max_predictability": row.get("max_predictability", ""),
                    "tx_expression": row.get("tx_expression", ""),
                    "tx_psi": row.get("tx_psi", ""),
                    "level3_pass": bool(l3row.iloc[0]["level3_pass"]) if len(l3row) and "level3_pass" in l3row else "",
                    "n_gwas_points": int(len(gwas)),
                    "n_caudate_meqtl_tests": int(len(c_meq)),
                    "n_caudate_meqtl_fdr": int(c_meq["significant_fdr"].astype(bool).sum()) if len(c_meq) else 0,
                }
            ],
        )
        manifest.append(
            {
                "rank": int(row["rank"]),
                "index_snp": snp,
                "outdir": str(loc_dir),
                "n_gwas_points": int(len(gwas)),
                "n_forest_rows": int(len(forest)),
                "hero": snp in {"rs8048039", "rs13331198"},
            }
        )
        print(f"Prepared {snp}: GWAS n={len(gwas)} meQTL n={len(c_meq)} forest n={len(forest)}")

    write_tsv(OUT / "locus_panel_manifest.tsv", manifest)
    print(f"Wrote manifest {OUT / 'locus_panel_manifest.tsv'}")


if __name__ == "__main__":
    main()
