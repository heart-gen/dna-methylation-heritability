#!/usr/bin/env python3
"""06_partitioned_heritability -- munge the frozen trait list to LDSC format.

Usage:
    python 04_munge_sumstats.py --run-id <id> [--trait scz]

The legacy step_1.sh hardcoded each trait as its own shell block, with the
column mapping inline and a `continue` on any missing file. That made the tested
set a property of the script rather than of a declared configuration, and let a
trait vanish from the FDR family without anything recording that it had.

Here the trait list, its column mappings and its sample sizes all live in
config/partitioned_heritability.yml, and a trait that fails to munge stops the
run: the FDR family is frozen at run open and may not shrink afterwards.
"""
from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

# GWAS Catalog harmonised files ship both `hm_*` and raw columns; LDSC's header
# cleaner maps several of each onto the same internal name and then aborts on
# the collision. Dropping the raw twin of any harmonised column it duplicates is
# what the legacy per-trait heredocs did, generalised.
HARMONISED_PAIRS = [
    ("hm_rsid", "variant_id"), ("hm_chrom", "chromosome"),
    ("hm_pos", "base_pair_location"), ("hm_other_allele", "other_allele"),
    ("hm_effect_allele", "effect_allele"), ("hm_beta", "beta"),
    ("hm_odds_ratio", "odds_ratio"), ("hm_ci_lower", "ci_lower"),
    ("hm_ci_upper", "ci_upper"),
    ("hm_effect_allele_frequency", "effect_allele_frequency"),
    ("rsid", "variant_id"),
]

FLAG_FOR = {
    "snp": "--snp", "a1": "--a1", "a2": "--a2",
    "signed_sumstats": "--signed-sumstats", "p": "--p", "frq": "--frq",
    "n": "--N", "n_col": "--N-col", "n_cas": "--N-cas", "n_con": "--N-con",
    "n_cas_col": "--N-cas-col", "n_con_col": "--N-con-col",
    "skip_rows": "--skip-rows",
}


def repo_root() -> Path:
    root = os.environ.get("V2_REPO_ROOT")
    if root:
        return Path(root)
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").is_dir():
            return parent
    raise SystemExit("Could not locate repository root")


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def drop_duplicate_harmonised(src: Path, dst: Path) -> Path:
    """Write `src` without the raw twin of any harmonised column present."""
    if dst.exists():
        return dst
    opener = gzip.open if str(dst).endswith(".gz") else open
    with _open(src) as fin, opener(dst, "wt") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        drop = {raw for keep, raw in HARMONISED_PAIRS
                if keep in header and raw in header}
        keep_idx = [i for i, c in enumerate(header) if c not in drop]
        fout.write("\t".join(header[i] for i in keep_idx) + "\n")
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            fout.write("\t".join(parts[i] for i in keep_idx) + "\n")
    return dst


def munge_one(trait: dict, run_dir: Path, cfg: dict, script_dir: Path) -> Path:
    out_dir = run_dir / "sumstats"
    out_dir.mkdir(parents=True, exist_ok=True)
    name = trait["name"]
    out_prefix = out_dir / name
    final = out_dir / f"{name}.sumstats.gz"
    if final.exists():
        print(f"[06] {name}: already munged")
        return final

    src = Path(trait["file"])
    if not src.exists():
        raise SystemExit(f"{name}: declared GWAS file missing: {src}")

    if trait.get("drop_duplicate_harmonised_columns"):
        suffix = ".tsv.gz" if str(src).endswith(".gz") else ".tsv"
        src = drop_duplicate_harmonised(src, out_dir / f"{name}.deduped{suffix}")

    ref = cfg["ld_references"][cfg["ld_reference_arm"]]
    cmd = [sys.executable, str(script_dir / "ldsc_wrapper.py"),
           cfg["ldsc_dir"], str(script_dir / "munge_sumstats.py"),
           "--sumstats", str(src), "--out", str(out_prefix),
           "--merge-alleles", ref["hapmap3_snplist"],
           "--chunksize", str(cfg["munge_chunksize"])]
    for key, value in trait["munge"].items():
        flag = FLAG_FOR.get(key)
        if flag is None:
            raise SystemExit(f"{name}: unknown munge key '{key}'")
        cmd += [flag, str(value)]

    print(f"[06] munging {name}", flush=True)
    res = subprocess.run(cmd, cwd=str(run_dir))
    if res.returncode != 0 or not final.exists():
        raise SystemExit(
            f"{name}: munge_sumstats failed (exit {res.returncode}). The trait "
            f"list is the FDR family and may not be reduced to work around a "
            f"failure.")
    return final


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--trait", default=None,
                    help="munge only this trait (default: all)")
    args = ap.parse_args()

    root = repo_root()
    script_dir = Path(__file__).resolve().parent
    run_dir = root / "06_partitioned_heritability" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")

    cfg = yaml.safe_load((root / "config" / "partitioned_heritability.yml").read_text())
    traits = cfg["traits"]
    if args.trait:
        traits = [t for t in traits if t["name"] == args.trait]
        if not traits:
            raise SystemExit(f"Trait '{args.trait}' is not in the frozen list")

    rows = []
    for trait in traits:
        path = munge_one(trait, run_dir, cfg, script_dir)
        with gzip.open(path, "rt") as f:
            n_snps = sum(1 for _ in f) - 1
        rows.append({"trait": trait["name"], "label": trait["label"],
                     "class": trait["class"], "n_snps": n_snps,
                     "sumstats": path.name})
        print(f"[06] {trait['name']}: {n_snps} SNPs")

    pd.DataFrame(rows).to_csv(run_dir / "sumstats" / "munge-summary.tsv",
                              sep="\t", index=False)


if __name__ == "__main__":
    main()
