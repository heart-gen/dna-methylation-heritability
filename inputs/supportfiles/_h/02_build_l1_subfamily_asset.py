#!/usr/bin/env python3
"""Build the LINE/L1 subfamily asset for Module 04's subfamily resolution.

Usage:
    python -I inputs/supportfiles/_h/02_build_l1_subfamily_asset.py \
        --rmsk <path to UCSC hg38 rmsk.txt.gz>

The project's RepeatMasker asset (`repeat-masker-hg38.gz`) is six columns:
chrom, start, end, repName, swScore, strand. That is enough for age classes,
which are a property of the subfamily name, but not for completeness: whether
an element is full-length is a property of WHERE IN ITS CONSENSUS it aligns,
and that lives in the repStart / repEnd / repLeft columns of the UCSC table the
asset was cut from. This stage reads that table, proves it is the same release
as the project asset, and writes one row per L1 element with the classes that
config/l1_subfamilies.yml declares.

Identity is proven, not assumed. The asset is re-derived from the table
(chrom, genoStart, genoEnd, repName, swScore, strand), sorted, and its MD5 must
equal the sorted MD5 of `repeat-masker-hg38.gz`; the LINE/L1 subset
(repClass LINE, repFamily L1, repName ^L1) must likewise equal
`repeat-masker.LINE_L1.hg38.bed.gz`, the file every accepted Module 04 run
read. The second check is what lets the subfamily classes partition the
accepted `line_l1_frac` exactly, which Module 04 stage 10 then verifies on
every VMR.

Consensus coordinates. UCSC stores them strand-dependently:
    + strand: consensus span repStart..repEnd, bases left after = -repLeft
    - strand: consensus span repLeft..repEnd,  bases left after = -repStart
so consensus length = repEnd + bases-left-after. It is not constant within a
subfamily (a handful of rows disagree), so the subfamily's consensus length is
its MODAL value, and that choice is recorded in the report.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml


def repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").is_dir():
            return parent
    raise SystemExit("Could not locate repository root")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sorted_md5(lines) -> str:
    """MD5 of the C-locale sorted lines, computed by `sort` to bound memory."""
    p = subprocess.run(["sort", "-S", "2G"], input="".join(lines), text=True,
                       capture_output=True, check=True, env={"LC_ALL": "C"})
    return hashlib.md5(p.stdout.encode()).hexdigest()


def file_sorted_md5(path: Path) -> str:
    with gzip.open(path, "rt") as fh:
        return sorted_md5(fh.readlines())


def classify(name: str, rules: list) -> str:
    for r in rules:
        if re.search(r["pattern"], name):
            return r["name"]
    raise SystemExit(f"repName {name!r} matches no age class; the classes in "
                     "config/l1_subfamilies.yml must partition the L1 asset")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rmsk", required=True, type=Path)
    args = ap.parse_args()

    root = repo_root()
    cfg = yaml.safe_load((root / "config" / "l1_subfamilies.yml").read_text())
    src = cfg["source"]
    out_dir = root / "inputs" / "supportfiles" / "_m"
    asset = root / src["project_asset"]
    l1_asset = root / src["project_l1_asset"]

    if sha256(args.rmsk) != src["rmsk_sha256"]:
        raise SystemExit(f"{args.rmsk} does not have the SHA-256 recorded in "
                         "config/l1_subfamilies.yml:source.rmsk_sha256")

    all6, l1_6, rows = [], [], []
    with gzip.open(args.rmsk, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, gs, ge, strand, name, rclass, rfam = (
                f[5], int(f[6]), int(f[7]), f[9], f[10], f[11], f[12])
            six = f"{chrom}\t{gs}\t{ge}\t{name}\t{f[1]}\t{strand}\n"
            all6.append(six)
            if rclass == "LINE" and rfam == "L1" and name.startswith("L1"):
                l1_6.append(six)
                rs, re_, rl = int(f[13]), int(f[14]), int(f[15])
                if strand == "+":
                    cstart, cend, cleft = rs, re_, -rl
                else:
                    cstart, cend, cleft = rl, re_, -rs
                rows.append((chrom, gs, ge, name, strand, cstart, cend, cend + cleft))

    md5_asset = file_sorted_md5(asset)
    md5_l1 = file_sorted_md5(l1_asset)
    md5_all6 = sorted_md5(all6)
    md5_l16 = sorted_md5(l1_6)
    if md5_all6 != md5_asset:
        raise SystemExit(f"rmsk table != {asset.name} (sorted md5 {md5_all6} vs "
                         f"{md5_asset}); not the release the project used")
    if md5_l16 != md5_l1:
        raise SystemExit(f"LINE/L1 subset != {l1_asset.name} (sorted md5 "
                         f"{md5_l16} vs {md5_l1})")
    del all6, l1_6

    # Modal consensus length per subfamily.
    by_name = defaultdict(Counter)
    for r in rows:
        by_name[r[3]][r[7]] += 1
    cons_len = {n: c.most_common(1)[0][0] for n, c in by_name.items()}

    age_rules = cfg["age_classes"]
    comp = cfg["completeness"]
    min_cons = int(comp["min_consensus_length_bp"])
    min_cov = float(comp["full_length_min_consensus_fraction"])
    five_max = int(comp["five_prime_max_consensus_start"])

    out = out_dir / "repeat-masker.LINE_L1.subfamily.hg38.tsv.gz"
    tally = Counter()
    bp = Counter()
    with gzip.open(out, "wt") as fh:
        fh.write("chrom\tstart\tend\trepName\tstrand\tage_class\tconsensus_length"
                 "\tconsensus_start\tconsensus_end\tconsensus_fraction"
                 "\tfull_length\tretains_5prime\n")
        for chrom, gs, ge, name, strand, cs, ce, _ in rows:
            L = cons_len[name]
            frac = (ce - cs + 1) / L if L > 0 else 0.0
            evaluable = L >= min_cons
            full = evaluable and frac >= min_cov
            five = evaluable and cs <= five_max
            age = classify(name, age_rules)
            fh.write(f"{chrom}\t{gs}\t{ge}\t{name}\t{strand}\t{age}\t{L}\t{cs}"
                     f"\t{ce}\t{frac:.4f}\t{str(full).upper()}\t{str(five).upper()}\n")
            key = (age, "full_length" if full else "fragment")
            tally[key] += 1
            bp[key] += ge - gs

    short = sorted(n for n, L in cons_len.items() if L < min_cons)
    rep = out_dir / "l1-subfamily-asset-report.tsv"
    with open(rep, "w") as fh:
        fh.write("field\tvalue\n")
        for k, v in [
            ("rmsk_source_url", src["rmsk_url"]),
            ("rmsk_last_modified", src["rmsk_last_modified"]),
            ("rmsk_sha256", src["rmsk_sha256"]),
            ("sorted_md5_rmsk_six_column", md5_all6),
            ("sorted_md5_project_asset", md5_asset),
            ("sorted_md5_rmsk_l1_subset", md5_l16),
            ("sorted_md5_project_l1_asset", md5_l1),
            ("n_l1_elements", len(rows)),
            ("n_subfamilies", len(cons_len)),
            ("consensus_length_rule", "modal (repEnd + bases-left-after) per repName"),
            ("n_subfamilies_consensus_below_min", len(short)),
            ("subfamilies_consensus_below_min", ",".join(short)),
            ("output_sha256", sha256(out)),
        ]:
            fh.write(f"{k}\t{v}\n")
        for (age, c), n in sorted(tally.items()):
            fh.write(f"n_elements.{age}.{c}\t{n}\n")
            fh.write(f"bp.{age}.{c}\t{bp[(age, c)]}\n")
    print(f"wrote {out} ({len(rows)} elements) and {rep}", file=sys.stderr)


if __name__ == "__main__":
    main()
