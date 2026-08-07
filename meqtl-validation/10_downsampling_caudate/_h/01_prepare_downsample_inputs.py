#!/usr/bin/env python3
"""
Prepare per-replicate inputs for official caudate TensorQTL downsampling.

Reuses Experiment 3 sample lists (shared3-aware, N=111, 30 reps, seed 20260805).
Writes per-rep covariates and phenotype BED.gz subsets under
meqtl-validation/10_downsampling_caudate/_m/repXXX/.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
PHASE1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping/caudate/_m"
EXP3_LISTS = (
    PROJECT / "meqtl-validation/04_cross_region_sharing/_m/caudate_downsample"
)
OUTDIR = PROJECT / "meqtl-validation/10_downsampling_caudate/_m"
SEED = 20260805


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=str(OUTDIR))
    p.add_argument("--source-lists", default=str(EXP3_LISTS))
    p.add_argument("--max-reps", type=int, default=0, help="0 = all replicates in source manifest")
    p.add_argument("--skip-phenotypes", action="store_true", help="Only lists + covariates")
    return p.parse_args()


def load_covariates(path: Path) -> pd.DataFrame:
    cov = pd.read_csv(path, sep="\t", index_col=0)
    # Ensure samples as index
    if not cov.index.astype(str).str.startswith("Br").any() and cov.columns.astype(str).str.startswith("Br").any():
        cov = cov.T
    cov.index = cov.index.astype(str)
    return cov


def subset_phenotype_bed(src: Path, samples: list[str], dest: Path) -> None:
    """Stream-subset phenotype BED to selected sample columns; write bgzip-ready BED then gzip."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix("")  # .bed if dest ends with .bed.gz
    if dest.name.endswith(".bed.gz"):
        tmp = dest.with_name(dest.name.replace(".bed.gz", ".bed"))

    opener = gzip.open if str(src).endswith(".gz") else open
    sample_set = set(samples)
    with opener(src, "rt") as fin, tmp.open("w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        # Expect #chr start end phenotype_id samples...
        if len(header) < 5:
            raise SystemExit(f"Unexpected phenotype BED header in {src}")
        meta = header[:4]
        sample_cols = header[4:]
        keep_idx = [i for i, s in enumerate(sample_cols) if str(s) in sample_set]
        kept_names = [sample_cols[i] for i in keep_idx]
        missing = [s for s in samples if s not in set(kept_names)]
        if missing:
            raise SystemExit(f"Phenotype BED missing {len(missing)} samples e.g. {missing[:5]}")
        # Preserve requested sample order
        order = {s: j for j, s in enumerate(kept_names)}
        ordered_idx = [keep_idx[order[s]] for s in samples]
        fout.write("\t".join(meta + samples) + "\n")
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            vals = [parts[4 + i] for i in ordered_idx]
            fout.write("\t".join(parts[:4] + vals) + "\n")

    # bgzip + tabix via shell tools if available; else gzip
    import subprocess

    if dest.exists():
        dest.unlink()
    try:
        subprocess.run(["bgzip", "-f", str(tmp)], check=True)
        gz = Path(str(tmp) + ".gz")
        if gz != dest:
            gz.rename(dest)
        subprocess.run(["tabix", "-f", "-p", "bed", str(dest)], check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        # Fallback gzip without tabix
        with tmp.open("rb") as fin, gzip.open(dest, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        tmp.unlink(missing_ok=True)
        print(f"WARNING: bgzip/tabix unavailable; wrote gzip-only {dest}")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    src = Path(args.source_lists)
    lists_src = src / "sample_lists"
    lists_dst = outdir / "sample_lists"
    lists_dst.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(src / "downsample_replicate_manifest.tsv", sep="\t")
    design = pd.read_csv(src / "downsample_design_summary.tsv", sep="\t")
    if args.max_reps > 0:
        manifest = manifest.head(args.max_reps).copy()

    pheno_src = PHASE1 / "prepared/cpg_phenotypes.all_autosomes.bed.gz"
    cov_src = PHASE1 / "prepared/covariates.txt"
    if not pheno_src.exists():
        raise SystemExit(f"Missing phenotype BED: {pheno_src}")
    if not cov_src.exists():
        raise SystemExit(f"Missing covariates: {cov_src}")
    cov_all = load_covariates(cov_src)

    rows = []
    for _, mrow in manifest.iterrows():
        rep = int(mrow["replicate"])
        src_list = Path(mrow["path"])
        if not src_list.exists():
            # fallback to relative under source lists
            src_list = lists_src / f"caudate_downsample_rep{rep:03d}.txt"
        samples = [s.strip() for s in src_list.read_text().splitlines() if s.strip()]
        if len(samples) != int(mrow["n_samples"]):
            raise SystemExit(f"rep {rep}: expected {mrow['n_samples']} samples, got {len(samples)}")

        dst_list = lists_dst / f"caudate_downsample_rep{rep:03d}.txt"
        shutil.copy2(src_list, dst_list)

        rep_dir = outdir / f"rep{rep:03d}"
        rep_dir.mkdir(parents=True, exist_ok=True)
        # Covariates: samples as rows (id), matching prepared orientation
        cov = cov_all.loc[samples].copy()
        cov.index.name = "id"
        cov_path = rep_dir / "covariates.txt"
        cov.to_csv(cov_path, sep="\t")

        pheno_path = rep_dir / "cpg_phenotypes.bed.gz"
        if not args.skip_phenotypes:
            print(f"rep {rep:03d}: writing phenotype BED ({len(samples)} samples)...")
            subset_phenotype_bed(pheno_src, samples, pheno_path)
        elif not pheno_path.exists():
            print(f"WARNING rep {rep}: --skip-phenotypes and missing {pheno_path}")

        rows.append({
            "replicate": rep,
            "n_samples": len(samples),
            "n_shared3_kept": int(mrow["n_shared3_kept"]),
            "n_fill": int(mrow["n_fill"]),
            "sample_list": str(dst_list),
            "covariates": str(cov_path),
            "phenotype_bed": str(pheno_path),
            "outdir": str(rep_dir / "tensorqtl"),
            "prefix": f"cpg_meqtl_caudate_ds_rep{rep:03d}",
            "seed": SEED,
            "source_list": str(src_list),
        })
        print(f"rep {rep:03d}: prepared under {rep_dir}")

    write_tsv(outdir / "downsample_replicate_manifest.tsv", rows)
    # Refresh design with official-remap note
    d = design.iloc[0].to_dict()
    d["method_note"] = (
        "Official TensorQTL cis permutation FDR remap of Exp3 sample lists "
        "(module 10_downsampling_caudate); covariates=locked M3a"
    )
    d["n_reps"] = len(rows)
    d["module"] = "10_downsampling_caudate"
    write_tsv(outdir / "downsample_design_summary.tsv", [d])
    # Shared3 audit copy
    shared = src / "shared3_donors.txt"
    if shared.exists():
        shutil.copy2(shared, outdir / "shared3_donors.txt")
    print(f"Prepared {len(rows)} replicates under {outdir}")


if __name__ == "__main__":
    main()
