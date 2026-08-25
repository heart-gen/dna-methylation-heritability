#!/usr/bin/env python3
"""Prepare per-chromosome tensorqtl inputs for 05_cpg_meqtl_burden.

Usage:
    python _h/01b_prepare_meqtl_inputs.py --run-id cmb-AA-caudate-20260823 --chrom 22

01_prepare_cpg_set.R decides WHICH CpGs are tested and writes the denominator
table that every burden fraction must use. This stage turns that decision into
the three files tensorqtl needs, for one autosome:

    inputs/chr{N}.phenotype.bed.gz   tested CpGs x donors, with positions
    inputs/chr{N}.covariates.tsv     covariates x donors
    inputs/chr{N}.{pgen,pvar,psam}   cis genotypes under the locked genotype QC

Splitting by chromosome is not just parallelism: a cis-meQTL scan never crosses
a chromosome, so per-chromosome genotype subsets keep peak memory proportional
to one chromosome rather than to the whole genome.

Every QC constant comes from config/meqtl_parameters.yml, which is prespecified.
Nothing here may be tuned after results are seen.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


def repo_root() -> Path:
    d = Path(__file__).resolve()
    while d != d.parent:
        if (d / ".git").is_dir():
            return d
        d = d.parent
    raise SystemExit("Could not locate repository root")


def read_manifest(run_dir: Path) -> dict:
    m = pd.read_csv(run_dir / "manifest.tsv", sep="\t", dtype=str)
    return dict(zip(m["field"], m["value"]))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--threads", default="1")
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "05_cpg_meqtl_burden" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")
    man = read_manifest(run_dir)
    cohort, region = man["cohort"], man["region"]

    cfg = yaml.safe_load((root / "config" / "meqtl_parameters.yml").read_text())
    gq = cfg["genotype_qc"]
    chrom = str(args.chrom)
    chrom_label = f"chr{chrom}"

    out_dir = run_dir / "inputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------ phenotypes
    member = pd.read_csv(run_dir / "results" / "tested-cpg-membership.tsv", sep="\t")
    member = member[member["chr"].astype(str) == chrom_label]
    if member.empty:
        # An autosome with no tested CpG is a real, recordable outcome (small
        # chromosome, heavy masking). Leave a marker so 02 skips it explicitly
        # rather than failing, and so the reconciliation can see why.
        (out_dir / f"{chrom_label}.no-tested-cpgs").write_text(
            "no tested CpG on this chromosome\n")
        print(f"{chrom_label}: no tested CpGs; nothing to prepare")
        return

    # 01_prepare_cpg_set.R already emitted the TESTED methylation for this
    # chromosome, donors x CpGs, keyed chr:pos. Read that, not the ~2M-column
    # Module 01 matrix -- reading the wide file here was OOM-killed, and it
    # would also duplicate the chr:pos identifier logic in a second language.
    meth_f = run_dir / "results" / "tested_meth" / f"{chrom_label}.tsv"
    if not meth_f.is_file():
        raise SystemExit(
            f"Missing tested methylation for {chrom_label}: {meth_f}. "
            "Run 01_prepare_cpg_set.R first.")

    meth = pd.read_csv(meth_f, sep="\t", dtype={"FID": str})
    donors = meth["FID"].astype(str).tolist()
    keep = [c for c in member["cpg_id"].astype(str) if c in meth.columns]
    if not keep:
        raise SystemExit(
            f"{chrom_label}: no tested CpG column matched the membership table")

    pos = dict(zip(member["cpg_id"].astype(str), member["cpg_pos"].astype(int)))
    mat = meth[keep].T
    mat.columns = donors
    bed = pd.DataFrame({
        "#chr": chrom_label,
        "start": [pos[c] - 1 for c in keep],     # BED is 0-based half-open
        "end": [pos[c] for c in keep],
        "phenotype_id": keep,
    })
    bed = pd.concat([bed.reset_index(drop=True), mat.reset_index(drop=True)], axis=1)
    bed = bed.sort_values(["start", "end"]).reset_index(drop=True)
    bed_f = out_dir / f"{chrom_label}.phenotype.bed"
    bed.to_csv(bed_f, sep="\t", index=False, na_rep="NA")
    # bgzip/tabix ship in the same conda env as tensorqtl, but PATH is not
    # necessarily set when the interpreter is invoked by absolute path (as the
    # step script does via `conda run`). Resolve them next to the interpreter.
    bindir = Path(sys.executable).parent
    bgzip = str(bindir / "bgzip") if (bindir / "bgzip").exists() else "bgzip"
    tabix = str(bindir / "tabix") if (bindir / "tabix").exists() else "tabix"
    subprocess.run([bgzip, "-f", str(bed_f)], check=True)
    subprocess.run([tabix, "-p", "bed", "-f", str(bed_f) + ".gz"], check=True)

    cat_dir = (root / "01_vmr_catalog" / "_m" / "runs"
               / man["upstream_vmr_catalog_run_id"])

    # ------------------------------------------------------------ covariates
    # The locked donor and covariate model (AGENTS.md 7.5: "preserve the locked
    # donor and covariate models"), taken from the SAME Module 01 files
    # 00_shared/locus_io.R reads, so 02, 03 and 05 share one covariate design.
    cov_dir = cat_dir / "covs" / f"chr_{chrom}"
    prefix = "TOPMed_LIBD.AA" if cohort == "AA" else "TOPMed_LIBD"
    covar = pd.read_csv(cov_dir / f"{prefix}.covar", sep=r"\s+", header=None,
                        names=["FID", "IID", "sex", "diagnosis"], dtype={0: str, 1: str})
    qcovar = pd.read_csv(cov_dir / f"{prefix}.qcovar", sep=r"\s+", header=None,
                         names=["FID", "IID", "age"], dtype={0: str, 1: str})
    cov = covar.merge(qcovar, on=["FID", "IID"])

    # Genotype PCs. Ancestry structure is a covariate for meQTL mapping exactly
    # as it is for the local-variance model.
    gen = yaml.safe_load((root / "config" / "paths.yml").read_text())["genotype"]
    arm = gen["AA"] if cohort == "AA" else gen["all_individuals"]
    evec_f = root / arm["eigenvec"]
    if evec_f.is_file():
        evec = pd.read_csv(evec_f, sep=r"\s+", dtype={0: str, 1: str})
        evec.columns = ["FID", "IID"] + [f"snpPC{i}" for i in
                                         range(1, evec.shape[1] - 1)]
        n_pc = 3
        cov = cov.merge(evec[["FID", "IID"] + [f"snpPC{i}" for i in
                                               range(1, n_pc + 1)]],
                        on=["FID", "IID"], how="left")
    else:
        print(f"WARNING: no eigenvec at {evec_f}; proceeding without genotype PCs",
              file=sys.stderr)

    cov = cov.set_index("FID").drop(columns=["IID"])
    cov = cov.loc[[d for d in donors if d in cov.index]]

    # tensorqtl builds its residualizer straight from the covariate values, so
    # the design must be fully numeric -- `sex` and `diagnosis` arrive as
    # strings ("F"/"M", "Control"/...). Dummy-code them here, dropping the
    # first level, which is the same treatment coding
    # 00_shared/locus_io.R applies for modules 02 and 03.
    categorical = [c for c in cov.columns if cov[c].dtype == object]
    if categorical:
        cov = pd.get_dummies(cov, columns=categorical, drop_first=True,
                             dtype=float)
    cov = cov.apply(pd.to_numeric, errors="coerce")

    if cov.isnull().any().any():
        bad = cov.columns[cov.isnull().any()].tolist()
        raise SystemExit(f"{chrom_label}: null covariate values in {bad}")
    # A constant column makes the residualizer rank-deficient.
    constant = [c for c in cov.columns if cov[c].nunique() < 2]
    if constant:
        print(f"{chrom_label}: dropping constant covariate(s) {constant}")
        cov = cov.drop(columns=constant)
    cov.T.to_csv(out_dir / f"{chrom_label}.covariates.tsv", sep="\t")

    # ------------------------------------------------------------- genotypes
    src = root / arm["pgen"]
    src_prefix = str(src)[: -len(".pgen")]
    geno_prefix = out_dir / f"{chrom_label}"
    plink2 = "/projects/p32505/opt/bin/plink2"

    # The distributed .psam has no header, so plink2 refuses it ("Line 1 has
    # fewer tokens than expected"). Stage a headered copy alongside symlinked
    # .pgen/.pvar rather than touching the shared source file, which other
    # analyses read.
    stage = out_dir / "geno_stage"
    stage.mkdir(exist_ok=True)
    staged_prefix = stage / f"{chrom_label}_src"
    for ext in ("pgen", "pvar"):
        link = staged_prefix.with_suffix(f".{ext}")
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(f"{src_prefix}.{ext}")
    psam_lines = ["#FID\tIID\tSEX"]
    for line in Path(f"{src_prefix}.psam").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        psam_lines.append(f"{f[0]}\t{f[1]}\t{f[2] if len(f) > 2 else 'NA'}")
    staged_prefix.with_suffix(".psam").write_text("\n".join(psam_lines) + "\n")

    keep_f = out_dir / f"{chrom_label}.keep"
    keep_f.write_text("".join(f"{d}\t{d}\n" for d in cov.index))

    cmd = [
        plink2, "--pfile", str(staged_prefix),
        "--chr", chrom,
        "--maf", str(gq["maf_min"]),
        "--geno", str(gq["missingness_max"]),
        "--hwe", str(gq["hwe_p_min"]),
        "--threads", str(args.threads),
        "--make-pgen", "--out", str(geno_prefix),
    ]
    print(" ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f"plink2 failed for {chrom_label}:\n{res.stderr[-4000:]}")

    print(f"{chrom_label}: {len(keep)} tested CpGs, {cov.shape[0]} donors, "
          f"{cov.shape[1]} covariates")


if __name__ == "__main__":
    main()
