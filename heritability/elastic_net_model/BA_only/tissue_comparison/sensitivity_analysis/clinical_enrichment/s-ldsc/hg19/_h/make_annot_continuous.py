#!/usr/bin/env python
from __future__ import print_function
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from pybedtools import BedTool


def _read_bim(bimfile):
    # FIX: delim_whitespace is deprecated in pandas 2.2, use sep=r"\s+" instead
    return pd.read_csv(
        bimfile,
        sep=r"\s+",
        usecols=[0, 1, 2, 3],
        names=["CHR", "SNP", "CM", "BP"],
    )


def _detect_bed_value_column(bed_file):
    # Read one row to determine number of columns (3=positional, 4+=continuous)
    sample = pd.read_csv(bed_file, sep="\t", header=None, nrows=1)
    return sample.shape[1] >= 4


def _prepare_bed(bed_file, nomerge, has_value, merge_op):
    bed = BedTool(bed_file).sort()
    if not nomerge:
        if has_value:
            # Preserve 4th column via merge; aggregate overlapping intervals
            bed = bed.merge(c=4, o=merge_op)
        else:
            bed = bed.merge()
    return bed


def _write_annot(df_annot, annot_file):
    if annot_file.endswith(".gz"):
        with gzip.open(annot_file, "wb") as f:
            df_annot.to_csv(f, sep="\t", index=False)
    else:
        df_annot.to_csv(annot_file, sep="\t", index=False)


def make_annot_files(args, beds, has_value_list):
    print("making combined annot file")

    df_bim = _read_bim(args.bimfile)

    iter_bim = [["chr" + str(x1), x2 - 1, x2] for (x1, x2) in np.array(df_bim[["CHR", "BP"]])]
    bimbed = BedTool(iter_bim)

    df_annot = pd.DataFrame(index=range(len(df_bim)))

    for i, (bed, has_value) in enumerate(zip(beds, has_value_list), start=1):
        print(f"Processing Q{i}")

        if has_value:
            mapped = bimbed.map(bed, c=4, o="sum")
            values = []

            for row in mapped:
                val = row[3]
                if val in (".", "", None):
                    values.append(0.0)
                else:
                    values.append(float(val))

            df_annot[f"ANNOT_Q{i}"] = values

        else:
            annotbed = bimbed.intersect(bed)
            bp = [x.start + 1 for x in annotbed]

            df_int = pd.DataFrame({"BP": bp, f"ANNOT_Q{i}": 1})

            df_tmp = pd.merge(df_bim[["BP"]], df_int, how="left", on="BP")
            df_tmp.fillna(0, inplace=True)

            df_annot[f"ANNOT_Q{i}"] = df_tmp[f"ANNOT_Q{i}"].astype(int)

    _write_annot(df_annot, args.annot_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
    "--bed-files",
    type=str,
    nargs="+",
    required=True,
    help="List of BED files (e.g. Q1.bed Q2.bed ... Q5.bed)"
    )
    parser.add_argument("--windowsize", type=int, default=0, help="(ignored for bed files; kept for compatibility)")
    parser.add_argument("--nomerge", action="store_true", default=False, help="do not merge the bed file")
    parser.add_argument("--merge-op", type=str, default="mean", choices=["mean", "sum", "max"], help="aggregation for overlapping intervals when bed has values")
    parser.add_argument("--bimfile", type=str, required=True, help="plink bim file used to compute LD scores")
    parser.add_argument("--annot-file", type=str, required=True, help="output annot file")

    args = parser.parse_args()

    beds = []
    has_value_list = []

    for bed_file in args.bed_files:
        has_value = _detect_bed_value_column(bed_file)
        bed = _prepare_bed(bed_file, args.nomerge, has_value, args.merge_op)

        beds.append(bed)
        has_value_list.append(has_value)

    make_annot_files(args, beds, has_value_list)


if __name__ == "__main__":
    main()
