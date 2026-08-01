#!/usr/bin/env python3
"""Merge per-chromosome CpG phenotype BEDs into one sorted autosomal BED."""

from __future__ import annotations

import gzip
import argparse
import subprocess
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-dir", required=True)
    parser.add_argument("--out-prefix", default="cpg_phenotypes.all_autosomes")
    args = parser.parse_args()

    prep = Path(args.prepared_dir)
    beds = sorted(prep.glob("cpg_phenotypes.chr*.bed")) + sorted(
        prep.glob("cpg_phenotypes.chr*.bed.gz")
    )
    # Prefer .bed.gz over .bed for the same chrom when both exist
    by_stem: dict[str, Path] = {}
    for bed in beds:
        stem = bed.name.replace(".bed.gz", "").replace(".bed", "")
        if stem.endswith(".chr"):
            continue
        if bed.name.endswith(".bed.gz") or stem not in by_stem:
            by_stem[stem] = bed
    beds = [by_stem[k] for k in sorted(by_stem)]
    if not beds:
        raise SystemExit(f"No prepared BEDs found under {prep}")

    out = prep / f"{args.out_prefix}.bed"
    header = None
    rows: list[str] = []
    for bed in beds:
        with open_text(bed) as handle:
            h = handle.readline()
            if header is None:
                header = h
            elif h != header:
                raise SystemExit(f"Sample header mismatch in {bed}")
            rows.extend(handle.readlines())

    rows.sort(key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
    with out.open("w") as handle:
        handle.write(header or "")
        handle.writelines(rows)
    print(f"Wrote {out} with {len(rows)} CpGs")

    try:
        subprocess.check_call(["bgzip", "-f", str(out)])
        gz = Path(str(out) + ".gz")
        subprocess.check_call(["tabix", "-p", "bed", str(gz)])
        print(f"Wrote {gz} + tabix index")
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"WARNING: bgzip/tabix unavailable or failed ({exc}); left {out}")


if __name__ == "__main__":
    main()
