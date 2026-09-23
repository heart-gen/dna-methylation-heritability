#!/usr/bin/env python3
"""Decision pilot: executed covariate model vs the locked M3a, one autosome.

NOT A PIPELINE STAGE. This script is deliberately not named `NN_*` and is not
referenced by `submit_meqtl_burden.sh` or any `step_*.sh`. It creates no run ID,
writes nothing under `_m/`, and seals nothing. Its only output is evidence for
the PI decision described below.

Why it exists
-------------
`config/covariates.yml:primary_meqtl` is `pi_locked` and states

    locked_model: M3a
    M3a = agedeath + sex + primarydx + snpPC1-5 + methPC1-5

with the lock rationale "M3a improves lambda_NS vs M0 and increases
FDR-significant CpG discoveries in all three regions".

`_h/01b_prepare_meqtl_inputs.py:145` sets `n_pc = 3` and adds no methylation
PCs, so the executed model is

    executed = agedeath + sex + primarydx + snpPC1-3

No gate compares the executed design against the lock, so the three accepted
runs `cmb-AA-{caudate,dlpfc,hippocampus}-20260825` were mapped under `executed`,
not under `M3a`. AGENTS.md 7.5 requires genomic inflation to be resolved before
the figure freeze and 12 makes the covariate model a PI decision, so the choice
between the two is not an agent's to make. This pilot supplies the numbers.

A second defect in the same file is measured here as well: lines 203-204 write
`{chrom}.keep` from the analysis donors, and the plink2 call at 206-214 never
passes `--keep`, so `--maf/--geno/--hwe` are evaluated over every donor in the
source pfile (526 AA donors) rather than the region's analysis set.

Design
------
Five arms on one autosome, all sharing the sealed run's tested-CpG set and
phenotype BED, so the only things that move are the covariate design and the
donor set used for genotype QC:

    executed        snpPC1-3, no methPC        sealed genotypes (no --keep)
    m0              snpPC1-5, no methPC        sealed genotypes (no --keep)
    m3a             snpPC1-5 + methPC1-5       sealed genotypes (no --keep)
    executed_keep   snpPC1-3, no methPC        genotypes rebuilt with --keep
    m3a_keep        snpPC1-5 + methPC1-5       genotypes rebuilt with --keep

`executed` is a reproduction check: its covariate file must be byte-identical to
the sealed run's, and its lambda must match the lambda recomputed from the
sealed nominal parquet.

methPC provenance
-----------------
`config/covariates.yml` names no file for methPC1-5; it specifies a method,
"PCA on M0-residualized CpG phenotypes". The v1 implementation of that method is
`meqtl-validation/01_cpg_meqtl_mapping/_h/09_estimate_latent_factors.py`, and
this script reimplements it unchanged in substance: M0 residualization, a
50,000-CpG subsample under the legacy seed, CpG standardization, PCA on
samples x CpGs. The v1 factors themselves are NOT reused, because they were
estimated on the v1 VMR catalog's CpGs and the v2 catalog is a different CpG
set.

Statistics
----------
The lambda is the module's own: `_h/04_qc_plots.py::genomic_inflation` over
nominal cis pairs with |CpG-to-SNP distance| > `distal_min_distance_bp`. It is
imported from that file rather than reimplemented. Significance is Storey
q <= `fdr_threshold` on the permutation p-values via `_h/storey_qvalue.py`, the
same helper `02b_combine_meqtl.R` uses; the FDR family here is the one
chromosome, identically in every arm.

Usage (run each stage under the `genomics` env, which holds tensorqtl):
    python _h/pilot_covariate_model.py methpcs   --out DIR [--run-id ...]
    python _h/pilot_covariate_model.py prepare   --out DIR --arm ARM --chrom N
    python _h/pilot_covariate_model.py map       --out DIR --arm ARM --chrom N
    python _h/pilot_covariate_model.py summarize --out DIR --chrom N
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

SEALED_RUN = "cmb-AA-caudate-20260825"
# Seed of the v1 latent-factor implementation, kept so the recipe is the locked
# one rather than a new one.
LATENT_SEED = 20260730
LATENT_SUBSAMPLE = 50_000
LATENT_MAX_K = 15

ARMS = {
    "executed":      {"n_snp_pc": 3, "meth_pc": 0, "keep": False},
    # M0 is the config's own named baseline, and it decomposes the executed-vs-
    # M3a gap, which moves TWO things at once: snpPC1-3 -> snpPC1-5 and no
    # methPC -> methPC1-5. Without it a lambda difference cannot be attributed
    # to the latent factors rather than to the two extra ancestry PCs.
    "m0":            {"n_snp_pc": 5, "meth_pc": 0, "keep": False},
    "m3a":           {"n_snp_pc": 5, "meth_pc": 5, "keep": False},
    "executed_keep": {"n_snp_pc": 3, "meth_pc": 0, "keep": True},
    "m3a_keep":      {"n_snp_pc": 5, "meth_pc": 5, "keep": True},
}

PLINK2 = "/projects/p32505/opt/bin/plink2"
GENOMICS_PY = "/projects/p32505/opt/envs/genomics/bin/python"


def repo_root() -> Path:
    d = Path(__file__).resolve()
    while d != d.parent:
        if (d / ".git").is_dir() or (d / ".git").is_file():
            return d
        d = d.parent
    raise SystemExit("Could not locate repository root")


ROOT = repo_root()
MODULE = ROOT / "05_cpg_meqtl_burden"
# Code and config come from this checkout; sealed runs and large inputs live in
# the canonical tree, which `config/paths.yml` names. In a git worktree they are
# not the same directory, and `_m/runs/` exists only in the canonical one.
DATA_ROOT = Path(yaml.safe_load(
    (ROOT / "config" / "paths.yml").read_text())["project_root"])


def load_h_module(name: str):
    """Import a stage from `_h/` so its helpers are the ones under test."""
    path = MODULE / "_h" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"h_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def sealed_dir(run_id: str) -> Path:
    d = DATA_ROOT / "05_cpg_meqtl_burden" / "_m" / "runs" / run_id
    if not d.is_dir():
        raise SystemExit(f"No sealed run at {d}")
    return d


def read_manifest(run_dir: Path) -> dict:
    m = pd.read_csv(run_dir / "manifest.tsv", sep="\t", dtype=str)
    return dict(zip(m["field"], m["value"]))


def base_covariates(run_dir: Path, man: dict, chrom: str,
                    donors: list[str]) -> pd.DataFrame:
    """age/sex/diagnosis exactly as `01b_prepare_meqtl_inputs.py` reads them."""
    cat_dir = (DATA_ROOT / "01_vmr_catalog" / "_m" / "runs"
               / man["upstream_vmr_catalog_run_id"])
    cov_dir = cat_dir / "covs" / f"chr_{chrom}"
    prefix = "TOPMed_LIBD.AA" if man["cohort"] == "AA" else "TOPMed_LIBD"
    covar = pd.read_csv(cov_dir / f"{prefix}.covar", sep=r"\s+", header=None,
                        names=["FID", "IID", "sex", "diagnosis"],
                        dtype={0: str, 1: str})
    qcovar = pd.read_csv(cov_dir / f"{prefix}.qcovar", sep=r"\s+", header=None,
                         names=["FID", "IID", "age"], dtype={0: str, 1: str})
    return covar.merge(qcovar, on=["FID", "IID"])


def snp_pcs(man: dict, n_pc: int) -> pd.DataFrame:
    gen = yaml.safe_load((ROOT / "config" / "paths.yml").read_text())["genotype"]
    arm = gen["AA"] if man["cohort"] == "AA" else gen["all_individuals"]
    evec = pd.read_csv(DATA_ROOT / arm["eigenvec"], sep=r"\s+", dtype={0: str, 1: str})
    evec.columns = ["FID", "IID"] + [f"snpPC{i}" for i in range(1, evec.shape[1] - 1)]
    return evec[["FID", "IID"] + [f"snpPC{i}" for i in range(1, n_pc + 1)]]


def build_covariates(run_dir: Path, man: dict, chrom: str, donors: list[str],
                     n_snp_pc: int, meth_pc: int,
                     methpc_file: Path | None) -> pd.DataFrame:
    """The covariate matrix of `01b_prepare_meqtl_inputs.py`, parameterized.

    With (n_snp_pc=3, meth_pc=0) this reproduces the executed design; with
    (5, 5) it is the locked M3a. Everything else -- the covar/qcovar source,
    the merge keys, the donor reindex, treatment dummy coding, the null and
    constant-column checks -- is copied from that stage so the two arms differ
    only in the covariate set.
    """
    cov = base_covariates(run_dir, man, chrom, donors)
    cov = cov.merge(snp_pcs(man, n_snp_pc), on=["FID", "IID"], how="left")
    cov = cov.set_index("FID").drop(columns=["IID"])
    cov = cov.loc[[d for d in donors if d in cov.index]]

    if meth_pc:
        if methpc_file is None or not methpc_file.is_file():
            raise SystemExit("methPC arm requested but no methPC table; run `methpcs` first")
        mp = pd.read_csv(methpc_file, sep="\t", index_col=0)
        mp.index = mp.index.astype(str)
        want = [f"methPC{i}" for i in range(1, meth_pc + 1)]
        missing = [d for d in cov.index if d not in mp.index]
        if missing:
            raise SystemExit(f"methPC table missing donors: {missing[:5]}")
        cov = cov.join(mp.loc[cov.index, want])

    categorical = [c for c in cov.columns if cov[c].dtype == object]
    if categorical:
        cov = pd.get_dummies(cov, columns=categorical, drop_first=True, dtype=float)
    cov = cov.apply(pd.to_numeric, errors="coerce")
    if cov.isnull().any().any():
        bad = cov.columns[cov.isnull().any()].tolist()
        raise SystemExit(f"null covariate values in {bad}")
    constant = [c for c in cov.columns if cov[c].nunique() < 2]
    if constant:
        print(f"dropping constant covariate(s) {constant}")
        cov = cov.drop(columns=constant)
    return cov


# --------------------------------------------------------------- stage: methpcs
def stage_methpcs(args) -> None:
    """methPC1-15 by the v1 recipe, on the v2 tested-CpG phenotypes.

    Reimplements meqtl-validation/01_cpg_meqtl_mapping/_h/09_estimate_latent_factors.py:
    residualize CpG methylation on M0 (age + sex + diagnosis + snpPC1-5), drop
    zero-variance CpGs, mean-impute, standardize CpGs, PCA on samples x CpGs.
    The CpG matrix is the sealed run's own per-chromosome phenotype BEDs pooled
    over the autosomes, so the factors describe the v2 CpG set.
    """
    run_dir = sealed_dir(args.run_id)
    man = read_manifest(run_dir)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    beds = sorted((run_dir / "inputs").glob("chr*.phenotype.bed.gz"))
    if not beds:
        raise SystemExit(f"No phenotype BEDs under {run_dir / 'inputs'}")

    frames, donors = [], None
    rng = np.random.default_rng(LATENT_SEED)
    n_total = 0
    for b in beds:
        df = pd.read_csv(b, sep="\t")
        df = df.set_index("phenotype_id").drop(columns=["#chr", "start", "end"])
        df.columns = df.columns.astype(str)
        if donors is None:
            donors = df.columns.tolist()
        elif df.columns.tolist() != donors:
            raise SystemExit(f"{b.name}: donor columns differ from {beds[0].name}")
        n_total += df.shape[0]
        frames.append(df.astype(np.float32))
    phen = pd.concat(frames)
    del frames
    print(f"pooled {len(beds)} autosomes: {phen.shape[0]} CpGs x {phen.shape[1]} donors")

    if LATENT_SUBSAMPLE and phen.shape[0] > LATENT_SUBSAMPLE:
        idx = rng.choice(phen.shape[0], size=LATENT_SUBSAMPLE, replace=False)
        phen = phen.iloc[idx]
        print(f"subsampled {LATENT_SUBSAMPLE} / {n_total} CpGs (seed {LATENT_SEED})")

    sd = phen.std(axis=1, ddof=0)
    phen = phen.loc[sd > 0]
    Y = phen.to_numpy(dtype=float)
    if np.isnan(Y).any():
        col_means = np.nanmean(Y, axis=1)
        inds = np.where(np.isnan(Y))
        Y[inds] = np.take(col_means, inds[0])

    # M0 design, in the same treatment coding the mapping stage uses.
    m0 = build_covariates(run_dir, man, "10", donors, n_snp_pc=5, meth_pc=0,
                          methpc_file=None)
    m0 = m0.loc[donors]
    C = m0.to_numpy(dtype=float)
    X = np.column_stack([np.ones(C.shape[0]), C])
    beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None)
    resid = (Y.T - X @ beta)          # samples x cpgs

    Z = (resid - resid.mean(axis=0)) / np.clip(resid.std(axis=0, ddof=0), 1e-8, None)
    from sklearn.decomposition import PCA
    pca = PCA(n_components=LATENT_MAX_K, random_state=LATENT_SEED)
    scores = pca.fit_transform(Z)
    factors = pd.DataFrame(scores, index=donors,
                           columns=[f"methPC{i}" for i in range(1, LATENT_MAX_K + 1)])
    factors.index.name = "FID"
    f_out = out / "methpcs.tsv"
    factors.to_csv(f_out, sep="\t")

    summary = {
        "source_run": args.run_id,
        "method": "PCA_on_M0_residuals",
        "m0_covariates": list(m0.columns),
        "n_donors": len(donors),
        "n_cpgs_available": n_total,
        "n_cpgs_used": int(phen.shape[0]),
        "subsample": LATENT_SUBSAMPLE,
        "seed": LATENT_SEED,
        "var_explained_k5": float(pca.explained_variance_ratio_[:5].sum()),
        "var_explained_k10": float(pca.explained_variance_ratio_[:10].sum()),
        "var_explained_k15": float(pca.explained_variance_ratio_[:15].sum()),
        "var_explained_each": [float(v) for v in pca.explained_variance_ratio_],
    }
    (out / "methpcs_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"wrote {f_out}")


# --------------------------------------------------------------- stage: prepare
def stage_prepare(args) -> None:
    run_dir = sealed_dir(args.run_id)
    man = read_manifest(run_dir)
    spec = ARMS[args.arm]
    chrom = str(args.chrom)
    lab = f"chr{chrom}"
    out = Path(args.out)
    arm_dir = out / args.arm
    arm_dir.mkdir(parents=True, exist_ok=True)

    sealed_bed = run_dir / "inputs" / f"{lab}.phenotype.bed.gz"
    for ext in ("", ".tbi"):
        dst = arm_dir / f"{lab}.phenotype.bed.gz{ext}"
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(str(sealed_bed) + ext)

    bed_head = pd.read_csv(sealed_bed, sep="\t", nrows=0)
    donors = [c for c in bed_head.columns if c not in ("#chr", "start", "end",
                                                       "phenotype_id")]
    cov = build_covariates(run_dir, man, chrom, donors,
                           spec["n_snp_pc"], spec["meth_pc"],
                           out / "methpcs.tsv")
    cov_f = arm_dir / f"{lab}.covariates.tsv"
    cov.T.to_csv(cov_f, sep="\t")
    print(f"{args.arm}: {cov.shape[0]} donors x {cov.shape[1]} covariates "
          f"({list(cov.columns)})")

    if not spec["keep"]:
        # Reuse the sealed genotype triple bit-for-bit: the no-keep arms must
        # differ from the accepted run in the covariate design alone.
        for ext in ("pgen", "pvar", "psam"):
            dst = arm_dir / f"{lab}.{ext}"
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(run_dir / "inputs" / f"{lab}.{ext}")
        print(f"{args.arm}: reusing sealed genotypes (no --keep), "
              f"{sum(1 for _ in open(run_dir / 'inputs' / f'{lab}.psam')) - 1} "
              "pfile samples")
        return

    # --keep arms: rebuild from the source pfile with the donor restriction the
    # sealed run computed and then discarded.
    gen = yaml.safe_load((ROOT / "config" / "paths.yml").read_text())["genotype"]
    g = gen["AA"] if man["cohort"] == "AA" else gen["all_individuals"]
    src_prefix = str(DATA_ROOT / g["pgen"])[: -len(".pgen")]
    stage = arm_dir / "geno_stage"
    stage.mkdir(exist_ok=True)
    staged = stage / f"{lab}_src"
    for ext in ("pgen", "pvar"):
        link = staged.with_suffix(f".{ext}")
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(f"{src_prefix}.{ext}")
    psam_lines = ["#FID\tIID\tSEX"]
    fid_to_iid: dict[str, str] = {}
    for line in Path(f"{src_prefix}.psam").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        psam_lines.append(f"{f[0]}\t{f[1]}\t{f[2] if len(f) > 2 else 'NA'}")
        fid_to_iid[f[0]] = f[1]
    staged.with_suffix(".psam").write_text("\n".join(psam_lines) + "\n")

    # The sealed stage wrote `{donor}\t{donor}`, i.e. the FID twice. plink2 reads
    # a two-column --keep file as FID/IID, and the real IID is a chip barcode
    # (e.g. 3998646007_R01C01), so that file matches no sample at all: had
    # --keep been passed, plink2 would have exited with no samples remaining.
    # A usable restriction has to pair each FID with its psam IID.
    missing = [d for d in cov.index if d not in fid_to_iid]
    if missing:
        raise SystemExit(f"donors absent from source psam: {missing[:5]}")
    keep_f = arm_dir / f"{lab}.keep"
    keep_f.write_text("".join(f"{d}\t{fid_to_iid[d]}\n" for d in cov.index))

    cfg = yaml.safe_load((ROOT / "config" / "meqtl_parameters.yml").read_text())
    gq = cfg["genotype_qc"]
    cmd = [
        PLINK2, "--pfile", str(staged),
        "--keep", str(keep_f),
        "--chr", chrom,
        "--maf", str(gq["maf_min"]),
        "--geno", str(gq["missingness_max"]),
        "--hwe", str(gq["hwe_p_min"]),
        "--threads", str(args.threads),
        "--make-pgen", "--out", str(arm_dir / lab),
    ]
    print(" ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f"plink2 failed:\n{res.stderr[-4000:]}")
    print(res.stdout[-2000:])


# ------------------------------------------------------------------- stage: map
def stage_map(args) -> None:
    """map_cis + map_nominal, using `_h/02_map_cpg_meqtl.py`'s own helpers."""
    h = load_h_module("02_map_cpg_meqtl")
    cfg = yaml.safe_load((ROOT / "config" / "meqtl_parameters.yml").read_text())
    window = int(cfg["cis_window_bp"])
    maf = float(cfg["genotype_qc"]["maf_min"])
    seed = int(cfg["mapping"]["seed"])

    lab = f"chr{args.chrom}"
    arm_dir = Path(args.out) / args.arm
    bed_f = arm_dir / f"{lab}.phenotype.bed.gz"
    cov_f = arm_dir / f"{lab}.covariates.tsv"
    geno_prefix = str(arm_dir / lab)

    import tensorqtl
    from tensorqtl import cis, genotypeio, pgen

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(str(bed_f))
    phenotype_df.columns = phenotype_df.columns.astype(str)

    probe = pgen.PgenReader(geno_prefix)
    geno_ids = set(map(str, probe.sample_ids))
    sample_ids = [s for s in phenotype_df.columns if s in geno_ids]
    if not sample_ids:
        raise SystemExit("No shared donors between phenotypes and genotypes")

    cov = pd.read_csv(cov_f, sep="\t", index_col=0).T
    cov.index = cov.index.astype(str)
    covariates_df = cov.loc[sample_ids].apply(pd.to_numeric, errors="raise")
    phenotype_df = phenotype_df[sample_ids]
    assert covariates_df.index.equals(phenotype_df.columns)

    pgr = pgen.PgenReader(geno_prefix, select_samples=sample_ids)
    pgr = h.prepare_chr_matched_genotypes(pgr)
    n_prepared = phenotype_df.shape[0]
    phenotype_df, phenotype_pos_df, dropped = h.filter_phenotypes_with_cis_variants(
        genotypeio, pgr, phenotype_df, phenotype_pos_df, window)
    print(f"{args.arm} {lab}: donors={len(sample_ids)} "
          f"covariates={covariates_df.shape[1]} prepared_cpgs={n_prepared} "
          f"tested_cpgs={phenotype_df.shape[0]} variants={pgr.num_variants}")

    genotype_df = pgr.load_genotypes()
    variant_df = h.normalize_variant_chrom(pgr.variant_df.copy())

    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                         covariates_df=covariates_df, maf_threshold=maf,
                         window=window, seed=seed)
    cis_df.index.name = "cpg_id"
    cis_df.to_csv(arm_dir / f"{lab}.cis_qtl.tsv.gz", sep="\t", float_format="%.6g")

    nominal_dir = arm_dir / "nominal"
    nominal_dir.mkdir(exist_ok=True)
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                    prefix=lab, covariates_df=covariates_df, maf_threshold=maf,
                    window=window, output_dir=str(nominal_dir))

    (arm_dir / f"{lab}.mapping_meta.json").write_text(json.dumps({
        "arm": args.arm, "chrom": lab, "n_donors": len(sample_ids),
        "n_covariates": int(covariates_df.shape[1]),
        "covariates": list(covariates_df.columns),
        "n_prepared_cpgs": int(n_prepared),
        "n_tested_cpgs": int(phenotype_df.shape[0]),
        "n_untested_no_cis_variant": len(dropped),
        "n_variants_in_pfile": int(pgr.num_variants),
    }, indent=2))


# ------------------------------------------------------------- stage: summarize
def storey_q(p: np.ndarray) -> tuple[np.ndarray, float]:
    script = MODULE / "_h" / "storey_qvalue.py"
    txt = "\n".join(format(x, ".17e") for x in p)
    res = subprocess.run([GENOMICS_PY, str(script)], input=txt,
                         capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f"storey_qvalue.py failed: {res.stderr[-2000:]}")
    lines = res.stdout.splitlines()
    pi0 = float(lines[0].split("\t")[1])
    return np.array([float(x) for x in lines[1:]]), pi0


def stage_summarize(args) -> None:
    qc = load_h_module("04_qc_plots")
    cfg = yaml.safe_load((ROOT / "config" / "meqtl_parameters.yml").read_text())
    distal_bp = int(cfg["genomic_inflation"]["distal_min_distance_bp"])
    fdr = float(cfg["mapping"]["fdr_threshold"])
    lab = f"chr{args.chrom}"
    out = Path(args.out)

    import pyarrow.parquet as pq
    rows = []
    for arm in ARMS:
        arm_dir = out / arm
        nom = arm_dir / "nominal" / f"{lab}.cis_qtl_pairs.{lab}.parquet"
        cisf = arm_dir / f"{lab}.cis_qtl.tsv.gz"
        if not nom.is_file() or not cisf.is_file():
            print(f"skipping {arm}: not mapped")
            continue
        meta = json.loads((arm_dir / f"{lab}.mapping_meta.json").read_text())

        pf = pq.ParquetFile(nom)
        allp, distp = [], []
        for b in pf.iter_batches(batch_size=4_000_000,
                                 columns=["pval_nominal", "start_distance"]):
            p = b.column("pval_nominal").to_numpy().astype(np.float32)
            d = np.abs(b.column("start_distance").to_numpy())
            allp.append(p)
            distp.append(p[d > distal_bp])
        allp = np.concatenate(allp)
        distp = np.concatenate(distp)
        lam = qc.genomic_inflation(distp.astype(np.float64))
        lam_all = qc.genomic_inflation(allp.astype(np.float64))

        cisd = pd.read_csv(cisf, sep="\t")
        cisd["pval_source"] = np.where(cisd["pval_beta"].notna(), "pval_beta",
                                       "pval_perm")
        cisd["pvalue"] = np.where(cisd["pval_beta"].notna(), cisd["pval_beta"],
                                  cisd["pval_perm"])
        cisd = cisd[np.isfinite(cisd["pvalue"])]
        q, pi0 = storey_q(cisd["pvalue"].to_numpy())
        n = len(cisd)
        order = np.argsort(cisd["pvalue"].to_numpy())
        pv = cisd["pvalue"].to_numpy()[order]
        bhq = np.minimum.accumulate((pv * n / np.arange(1, n + 1))[::-1])[::-1]

        rows.append({
            "arm": arm,
            "n_covariates": meta["n_covariates"],
            "covariates": ";".join(meta["covariates"]),
            "n_donors": meta["n_donors"],
            "n_variants_in_pfile": meta["n_variants_in_pfile"],
            "n_prepared_cpgs": meta["n_prepared_cpgs"],
            "n_tested_cpgs": meta["n_tested_cpgs"],
            "n_untested_no_cis_variant": meta["n_untested_no_cis_variant"],
            "n_nominal_pairs": int(allp.size),
            "n_distal_pairs": int(distp.size),
            "lambda_distal": lam,
            "lambda_all_pairs": lam_all,
            "storey_pi0": pi0,
            "n_significant_storey": int((q <= fdr).sum()),
            "n_significant_bh": int((bhq <= fdr).sum()),
            "median_nominal_p": float(np.median(allp)),
        })
        del allp, distp

    tab = pd.DataFrame(rows)
    f = out / f"pilot_summary_{lab}.tsv"
    tab.to_csv(f, sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 50):
        print(tab.drop(columns=["covariates"]).to_string(index=False))
    print(f"\nwrote {f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["methpcs", "prepare", "map", "summarize"])
    ap.add_argument("--out", required=True, help="scratch directory; never _m/")
    ap.add_argument("--arm", choices=sorted(ARMS))
    ap.add_argument("--chrom", default="10")
    ap.add_argument("--run-id", default=SEALED_RUN)
    ap.add_argument("--threads", default="8")
    args = ap.parse_args()

    out = Path(args.out).resolve()
    if "_m/runs" in str(out) or (MODULE / "_m") in out.parents:
        raise SystemExit("refusing to write under _m/: this is not a production run")

    if args.stage in ("prepare", "map") and not args.arm:
        raise SystemExit("--arm is required for prepare/map")
    {"methpcs": stage_methpcs, "prepare": stage_prepare,
     "map": stage_map, "summarize": stage_summarize}[args.stage](args)


if __name__ == "__main__":
    main()
