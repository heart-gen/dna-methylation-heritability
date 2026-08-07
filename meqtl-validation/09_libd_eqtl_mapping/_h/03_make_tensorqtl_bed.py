#!/usr/bin/env python3
"""Convert prepared LIBD expression + annotation to TensorQTL phenotype BED."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prepared-dir", required=True)
    p.add_argument("--outdir", default="")
    p.add_argument("--prefix", default="genes")
    p.add_argument("--keep-interval", action="store_true")
    args = p.parse_args()

    prepared = Path(args.prepared_dir)
    outdir = Path(args.outdir) if args.outdir else prepared
    outdir.mkdir(parents=True, exist_ok=True)

    expr = pd.read_csv(prepared / "normalized_expression.tsv.gz", sep="\t")
    annot = pd.read_csv(prepared / "feature.bed", sep="\t", dtype={"chrom": str})
    chroms = set(pd.read_csv(prepared / "vcf_chr_list.txt", header=None)[0].astype(str))

    annot = annot.rename(columns={"chrom": "chr"})
    bed = annot[["chr", "start", "end", "feature_id", "strand"]].copy()
    bed = bed[bed["chr"].isin(chroms)].copy()
    if not args.keep_interval:
        plus = bed["strand"].astype(str) != "-"
        bed.loc[plus, "end"] = bed.loc[plus, "start"] + 1
        bed.loc[~plus, "start"] = bed.loc[~plus, "end"]
        bed.loc[~plus, "end"] = bed.loc[~plus, "start"] + 1

    expr = expr.set_index("feature_id")
    bed = bed.merge(expr, left_on="feature_id", right_index=True, how="inner")
    bed = bed.drop(columns=["strand"])
    bed["chr_order"] = (
        bed["chr"].str.replace("chr", "", regex=False).replace({"X": "23"}).astype(int)
    )
    bed = bed.sort_values(["chr_order", "start", "end", "feature_id"]).drop(columns=["chr_order"])

    out_bed = outdir / f"{args.prefix}.expression.bed"
    header = list(bed.columns)
    header[0] = "#chr"
    bed.to_csv(out_bed, sep="\t", index=False, header=header)
    for exe in ("bgzip", "tabix"):
        if shutil.which(exe) is None:
            raise SystemExit(f"Required executable not found: {exe}")
    subprocess.check_call(["bgzip", "-f", str(out_bed)])
    subprocess.check_call(["tabix", "-f", str(out_bed) + ".gz"])
    print(str(out_bed) + ".gz")


if __name__ == "__main__":
    main()
