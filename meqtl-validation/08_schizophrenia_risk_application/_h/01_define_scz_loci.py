#!/usr/bin/env python3
"""Analysis 1: prespecify PGC3 schizophrenia-risk loci (hg38) + index SNPs.

Does not use methylation results. Index SNPs are matched to the analysis
genotype panel by rsID (or chr:pos_A1_A2) to obtain hg38 coordinates.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
DEFAULT_LOCI = Path(
    "/projects/b1213/resources/gwas/pgc3/fine_mapped_loci/to_hg38/_m/gwas_loci_hg38.tsv"
)
DEFAULT_INDEX = (
    PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m/raw"
    / "pgc3_primary_gwas_index_snps.tsv"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", default="caudate")
    p.add_argument("--population", default="AA", choices=["AA", "EA"])
    p.add_argument("--loci-tsv", default=str(DEFAULT_LOCI))
    p.add_argument("--index-tsv", default=str(DEFAULT_INDEX))
    p.add_argument("--outdir", default="")
    return p.parse_args()


def _parse_alleles(a1a2: str) -> tuple[str, str]:
    s = str(a1a2)
    if "/" in s:
        a1, a2 = s.split("/", 1)
        return a1.strip().upper(), a2.strip().upper()
    return "", ""


def _extract_rsid(variant_id: str) -> str:
    m = re.search(r"(rs\d+)$", str(variant_id))
    return m.group(1) if m else ""


def load_pvar(prefix: Path) -> pd.DataFrame:
    pvar = prefix.with_suffix(".pvar")
    df = pd.read_csv(
        pvar, sep="\t", comment="#", header=None,
        names=["CHROM", "POS", "ID", "REF", "ALT"],
    )
    df["chrom"] = "chr" + df["CHROM"].astype(str).str.replace("^chr", "", regex=True)
    df["rsid"] = df["ID"].map(_extract_rsid)
    df["pos_hg38"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
    return df


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir) if args.outdir else (
        PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m" / args.region
    )
    outdir.mkdir(parents=True, exist_ok=True)

    loci = pd.read_csv(args.loci_tsv, sep="\t")
    for c in ["chrom", "start", "end", "locus"]:
        if c not in loci.columns:
            raise SystemExit(f"loci table missing column {c}")
    loci["chrom"] = loci["chrom"].astype(str)
    if not str(loci["chrom"].iloc[0]).startswith("chr"):
        loci["chrom"] = "chr" + loci["chrom"].str.replace("^chr", "", regex=True)
    loci["start"] = pd.to_numeric(loci["start"], errors="coerce").astype(int)
    loci["end"] = pd.to_numeric(loci["end"], errors="coerce").astype(int)
    loci["locus_id"] = loci["locus"].astype(str)

    index = pd.read_csv(args.index_tsv, sep="\t")
    index = index.rename(columns={
        "SNP": "index_snp",
        "CHR": "chr_hg19",
        "BP": "bp_hg19",
        "P": "gwas_p",
        "OR": "gwas_or",
        "SE": "gwas_se",
        "A1A2": "a1a2",
        "range.left": "range_left_hg19",
        "range.right": "range_right_hg19",
        "span(kb)": "span_kb",
    })
    index["index_snp"] = index["index_snp"].astype(str)
    alleles = index["a1a2"].map(_parse_alleles)
    index["risk_allele"] = [a[0] for a in alleles]
    index["other_allele"] = [a[1] for a in alleles]

    geno = (
        PROJECT / "meqtl-validation/01_cpg_meqtl_mapping" / args.region / "_m"
        / "genotypes" / f"meqtl_{args.population}"
    )
    if not geno.with_suffix(".pvar").exists():
        raise SystemExit(f"Missing genotype pvar: {geno}.pvar")
    pvar = load_pvar(geno)

    # Match by rsID first
    by_rs = pvar[pvar["rsid"].ne("")].drop_duplicates("rsid", keep="first")
    matched = index.merge(
        by_rs[["rsid", "ID", "chrom", "pos_hg38", "REF", "ALT"]].rename(
            columns={"rsid": "index_snp", "ID": "genotype_variant_id"}
        ),
        on="index_snp",
        how="left",
    )

    # Fallback: chr:pos_A1_A2 style index IDs against pvar ID / chrom-pos
    miss = matched["genotype_variant_id"].isna()
    if miss.any():
        coord_re = re.compile(r"^(\d+|chr\d+|chrX):(\d+)_([ACGT]+)_([ACGT]+)$", re.I)
        for i in matched.index[miss]:
            snp = matched.at[i, "index_snp"]
            m = coord_re.match(snp)
            if not m:
                continue
            chrom = "chr" + m.group(1).replace("chr", "")
            pos = int(m.group(2))
            # These coordinates in the index table are typically hg19; skip direct pos match
            # unless an exact ID string appears in the pvar.
            hit = pvar[pvar["ID"].astype(str).str.contains(re.escape(snp), regex=True)]
            if hit.empty:
                continue
            row = hit.iloc[0]
            matched.at[i, "genotype_variant_id"] = row["ID"]
            matched.at[i, "chrom"] = row["chrom"]
            matched.at[i, "pos_hg38"] = row["pos_hg38"]
            matched.at[i, "REF"] = row["REF"]
            matched.at[i, "ALT"] = row["ALT"]

    # Assign locus_id by hg38 position within published locus intervals
    matched["locus_id"] = pd.NA
    for i, r in matched.dropna(subset=["chrom", "pos_hg38"]).iterrows():
        sub = loci[
            (loci["chrom"] == r["chrom"])
            & (loci["start"] <= int(r["pos_hg38"]))
            & (loci["end"] >= int(r["pos_hg38"]))
        ]
        if len(sub) == 1:
            matched.at[i, "locus_id"] = sub.iloc[0]["locus_id"]
        elif len(sub) > 1:
            # Prefer smallest spanning interval
            sub = sub.assign(span=sub["end"] - sub["start"]).sort_values("span")
            matched.at[i, "locus_id"] = sub.iloc[0]["locus_id"]

    matched["in_genotype_panel"] = matched["genotype_variant_id"].notna()
    matched["analysis_region"] = args.region
    matched["analysis_population"] = args.population

    keep = [
        "locus_id", "index_snp", "genotype_variant_id", "chrom", "pos_hg38",
        "REF", "ALT", "risk_allele", "other_allele", "gwas_p", "gwas_or", "gwas_se",
        "chr_hg19", "bp_hg19", "range_left_hg19", "range_right_hg19", "span_kb",
        "in_genotype_panel", "analysis_region", "analysis_population",
    ]
    keep = [c for c in keep if c in matched.columns]
    out_index = matched[keep].sort_values(["chrom", "pos_hg38"], na_position="last")
    out_index.to_csv(outdir / "scz_index_snps_hg38.tsv", sep="\t", index=False)

    # Locus table with optional index annotation
    loci_out = loci.copy()
    idx_one = (
        out_index.dropna(subset=["locus_id"])
        .drop_duplicates("locus_id", keep="first")
        [["locus_id", "index_snp", "genotype_variant_id", "chrom", "pos_hg38", "gwas_p", "gwas_or"]]
        .rename(columns={"chrom": "index_chrom", "pos_hg38": "index_pos_hg38"})
    )
    loci_out = loci_out.merge(idx_one, on="locus_id", how="left")
    loci_out.to_csv(outdir / "scz_risk_loci_hg38.tsv", sep="\t", index=False)

    write_tsv(outdir / "define_loci_summary.tsv", [{
        "region": args.region,
        "population": args.population,
        "n_loci": int(len(loci_out)),
        "n_index_snps": int(len(out_index)),
        "n_index_in_genotype_panel": int(out_index["in_genotype_panel"].sum()),
        "n_index_assigned_to_locus": int(out_index["locus_id"].notna().sum()),
        "loci_tsv": str(args.loci_tsv),
        "index_tsv": str(args.index_tsv),
        "genotype_prefix": str(geno),
    }])
    print(
        f"Wrote {outdir}/scz_risk_loci_hg38.tsv and scz_index_snps_hg38.tsv; "
        f"panel matches={int(out_index['in_genotype_panel'].sum())}/{len(out_index)}"
    )


if __name__ == "__main__":
    main()
