#!/usr/bin/env python3
"""Estimate the locked methylation latent factors (methPC) for one run.

Usage:
    python _h/01a_estimate_latent_factors.py --run-id cmb-AA-caudate-20260825

`config/covariates.yml:primary_meqtl` locks the primary meQTL model to

    M3a = agedeath + sex + primarydx + snpPC1-5 + methPC1-5

and specifies methPC1-5 as "PCA on M0-residualized CpG phenotypes". It names no
file, no CpG set, no subsample, no seed and no standardization, and it does not
say whether the M0 residualization carries the ancestry PCs. A method description
is not a specification, and AGENTS.md 9 requires a run to carry its configuration,
so every one of those choices is resolved HERE, explicitly, and written to
`results/latent-factor-provenance.tsv` with the factors themselves.

The resolution is the 2026-09-23 decision pilot's (PILOT_COVARIATE_MODEL.md):
the v1 implementation of the same method,
`meqtl-validation/01_cpg_meqtl_mapping/_h/09_estimate_latent_factors.py`, under
its own seed and 50,000-CpG subsample, re-fitted on the **v2** tested-CpG
phenotypes. The v1 factor tables are deliberately not reused -- they were
estimated on the v1 VMR catalog's CpGs, and the v2 catalog is a different CpG
set, so they would carry a different catalog's structure into this one.

Two things about these factors that must not be lost downstream:

  * methPC1 is largely collinear with this region's estimated cell composition
    (R^2 = 0.72 on RNA MuSiC proportions; Oligo rho = +0.77). Adopting M3a
    therefore imports a substantial cell-composition adjustment into the PRIMARY
    meQTL model. That is a bulk-tissue correlation between a methylation PC and an
    RNA-derived proportion estimate: it establishes collinearity, NOT a cell type
    of origin, which AGENTS.md 2.3 forbids inferring. Nothing in this module may
    describe the primary model as free of cell composition.
  * the factors are estimated from the same CpGs that are then tested, so they are
    not an independent covariate set. They are a latent-structure adjustment, and
    the lock's rationale is calibration (lambda) plus discovery, not independence.

Runs once per run, on the submit host, after `01_prepare_cpg_set.R` and before
the mapping array: every array task reads the one factor table, and 22 tasks
each estimating their own would be both wasteful and non-deterministic.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from meqtl_covariates import (  # noqa: E402
    CovariateLockError, build_covariate_matrix, latent_factor_recipe,
    load_covariates_config, locked_meqtl_design, read_manifest, repo_root,
    run_directory, vmr_catalog_dir,
)


def ordered_member_cpgs(run_dir: Path) -> pd.DataFrame:
    """Tested CpGs in numeric chromosome order, then position.

    AGENTS.md 10.1: "Chromosome ordering is numeric and deterministic." It matters
    here and not only as hygiene -- the subsample is drawn by INDEX, so the row
    order of the pooled matrix decides which CpGs are selected.
    """
    f = Path(run_dir) / "results" / "tested-cpg-membership.tsv"
    if not f.is_file():
        raise SystemExit(f"Missing {f}. Run 01_prepare_cpg_set.R first.")
    m = pd.read_csv(f, sep="\t", dtype={"chr": str, "cpg_id": str})
    m["chrom_n"] = pd.to_numeric(m["chr"].str.replace("^chr", "", regex=True),
                                 errors="coerce")
    off = m["chrom_n"].isna()
    if off.any():
        raise SystemExit(
            f"{int(off.sum())} tested CpG(s) on non-numeric chromosomes "
            f"{sorted(m.loc[off, 'chr'].unique())}; the chromosome policy is "
            "autosomal (config/meqtl_parameters.yml).")
    m["cpg_pos"] = m["cpg_pos"].astype(int)
    return m.sort_values(["chrom_n", "cpg_pos", "cpg_id"]).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    root = repo_root()
    run_dir = run_directory(root, args.run_id)
    man = read_manifest(run_dir)
    cohort, region = man["cohort"], man["region"]
    smoke = str(man.get("smoke_run", "FALSE")).upper() == "TRUE"

    cfg = load_covariates_config(root)
    design = locked_meqtl_design(cfg)
    recipe = latent_factor_recipe(cfg)
    print(f"[05.01a] locked design: {design.describe()}")
    print(f"[05.01a] latent-factor recipe from {recipe.source}: {recipe.values}")

    out_dir = run_dir / "results"
    factors_f = out_dir / "latent-factors.tsv"
    prov_f = out_dir / "latent-factor-provenance.tsv"

    if not design.latent_factors:
        ## A lock without latent factors is a legitimate configuration (M0), and
        ## the right behaviour is to record that no factor table is needed rather
        ## than to leave a stale one from an earlier design in place.
        if factors_f.exists():
            raise SystemExit(
                f"The lock ({design.model_id}) names no latent factor, but "
                f"{factors_f} exists. Refusing to leave a factor table that the "
                "design does not use where 01b would read it.")
        pd.DataFrame({"field": ["locked_model", "n_latent_factors"],
                      "value": [design.model_id, 0]}).to_csv(
            prov_f, sep="\t", index=False)
        print(f"[05.01a] {design.model_id} names no latent factor; nothing to do")
        return

    member = ordered_member_cpgs(run_dir)
    meth_dir = out_dir / "tested_meth"

    ## Pass 1: headers only, to fix the global CpG order and pick the subsample
    ## without ever holding the whole matrix. The tested_meth tables are written
    ## before the blacklist exclusion, so the membership table -- which 01b also
    ## uses -- is the authority on which CpGs are tested.
    per_chrom: dict[str, list[str]] = {}
    ordered: list[tuple[str, str]] = []
    for chrom, grp in member.groupby("chr", sort=False):
        f = meth_dir / f"{chrom}.tsv"
        if not f.is_file():
            print(f"[05.01a] {chrom}: no tested methylation table; skipped")
            continue
        have = set(pd.read_csv(f, sep="\t", nrows=0).columns)
        ids = [c for c in grp["cpg_id"] if c in have]
        if not ids:
            continue
        per_chrom[chrom] = ids
        ordered.extend((chrom, c) for c in ids)
    n_total = len(ordered)
    if n_total == 0:
        raise SystemExit("No tested CpG has methylation; cannot estimate factors")

    seed = recipe.seed
    n_sub = recipe.n_cpg_subsample
    rng = np.random.default_rng(seed)
    if n_sub and n_total > n_sub:
        ## Sorted so each chromosome is read once, contiguously. The SET is what
        ## the seed fixes; sorting the indices cannot change the factors, because
        ## permuting CpGs permutes the COLUMNS of the sample x CpG matrix and the
        ## sample-space principal components are invariant to that.
        sel = np.sort(rng.choice(n_total, size=n_sub, replace=False))
        chosen = [ordered[i] for i in sel]
        print(f"[05.01a] subsampled {n_sub} / {n_total} tested CpGs (seed {seed})")
    else:
        chosen = ordered
        print(f"[05.01a] using all {n_total} tested CpGs (below the subsample size)")

    by_chrom: dict[str, list[str]] = {}
    for chrom, cpg in chosen:
        by_chrom.setdefault(chrom, []).append(cpg)

    ## Pass 2: read only the selected columns.
    donors: list[str] | None = None
    blocks, kept_ids = [], []
    for chrom in per_chrom:
        ids = by_chrom.get(chrom)
        if not ids:
            continue
        df = pd.read_csv(meth_dir / f"{chrom}.tsv", sep="\t",
                         usecols=["FID"] + ids, dtype={"FID": str})
        d = df["FID"].astype(str).tolist()
        if donors is None:
            donors = d
        elif d != donors:
            raise SystemExit(
                f"{chrom}.tsv donor order differs from the first chromosome "
                "read; tested_meth tables must share one donor order.")
        blocks.append(df[ids].to_numpy(dtype=np.float64))
        kept_ids.extend(ids)
    assert donors is not None
    Y = np.concatenate(blocks, axis=1)          # samples x cpgs
    del blocks
    print(f"[05.01a] matrix {Y.shape[0]} donors x {Y.shape[1]} CpGs")

    sd = np.nanstd(Y, axis=0, ddof=0)
    keep = np.isfinite(sd) & (sd > 0)
    if not keep.all():
        print(f"[05.01a] dropping {int((~keep).sum())} zero-variance CpG(s)")
        Y = Y[:, keep]
        kept_ids = [c for c, k in zip(kept_ids, keep) if k]
    if Y.shape[1] == 0:
        raise SystemExit("Every selected CpG had zero variance")

    n_missing = int(np.isnan(Y).sum())
    if n_missing:
        col_means = np.nanmean(Y, axis=0)
        idx = np.where(np.isnan(Y))
        Y[idx] = np.take(col_means, idx[1])
        print(f"[05.01a] mean-imputed {n_missing} missing value(s)")

    ## M0 -- base terms + every locked ancestry PC -- is what the factors are
    ## residual to. build_covariate_matrix asserts that design too, so an M0 that
    ## quietly lost a PC cannot produce factors that then look fine.
    cat_dir = vmr_catalog_dir(root, man)
    first_chrom = str(int(member["chrom_n"].iloc[0]))
    m0, m0_record = build_covariate_matrix(
        root, cat_dir, cohort, first_chrom, donors, design,
        terms="m0", context=f"{args.run_id} latent-factor M0 design",
        allow_unlocked=smoke)
    X = np.column_stack([np.ones(m0.shape[0]), m0.to_numpy(dtype=float)])
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    del Y

    Z = resid - resid.mean(axis=0)
    Z = Z / np.clip(resid.std(axis=0, ddof=0), 1e-8, None)
    k = int(min(recipe.n_factors_estimated, Z.shape[0] - 1, Z.shape[1]))
    if k < design.n_latent_factors:
        raise SystemExit(
            f"Only {k} factor(s) are estimable from {Z.shape[0]} donors x "
            f"{Z.shape[1]} CpGs, but the lock names "
            f"{design.n_latent_factors}.")
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    ## Deterministic signs: a sign flip is free in an SVD and would otherwise
    ## differ between BLAS builds. Fix the largest-magnitude loading positive.
    for j in range(Vt.shape[0]):
        if Vt[j, np.argmax(np.abs(Vt[j]))] < 0:
            U[:, j] *= -1
            Vt[j] *= -1
    scores = U[:, :k] * S[:k]
    var_ratio = (S ** 2 / np.sum(S ** 2))[:k]

    factors = pd.DataFrame(scores, index=donors,
                           columns=[f"methPC{i}" for i in range(1, k + 1)])
    factors.index.name = "FID"
    factors.to_csv(factors_f, sep="\t", float_format="%.10g")

    used = list(design.latent_factors)
    prov = {
        "run_id": args.run_id,
        "region": region,
        "cohort": cohort,
        "locked_model": design.model_id,
        "lock_status": design.lock_status,
        "lock_date": design.lock_date,
        "config_covariates_sha256": design.config_sha256,
        "latent_factors_in_lock": ";".join(used),
        "n_factors_estimated": k,
        "recipe_source": recipe.source,
        "m0_design": ";".join(m0.columns),
        "m0_matches_lock": m0_record["matches_lock"],
        "n_donors": len(donors),
        "n_tested_cpgs_available": n_total,
        "n_cpgs_used": Z.shape[1],
        "n_missing_imputed": n_missing,
        "var_explained_locked_k": float(var_ratio[:len(used)].sum()),
        "var_explained_each": ";".join(f"{v:.6f}" for v in var_ratio),
        "cell_composition_note": (
            "methPC1 is largely collinear with RNA MuSiC cell proportions "
            "(R2 ~ 0.72 in caudate); the primary model is NOT free of cell "
            "composition. Collinearity only -- no cell type of origin is implied "
            "(AGENTS.md 2.3)."),
    }
    prov.update({f"recipe_{k2}": v for k2, v in recipe.values.items()})
    pd.DataFrame({"field": list(prov), "value": [str(v) for v in prov.values()]}
                 ).to_csv(prov_f, sep="\t", index=False)
    (out_dir / "latent-factor-summary.json").write_text(
        json.dumps({k2: (v if not isinstance(v, np.generic) else v.item())
                    for k2, v in prov.items()}, indent=2) + "\n")

    print(f"[05.01a] wrote {factors_f} ({k} factors; "
          f"{prov['var_explained_locked_k']:.4f} of residual variance in "
          f"{len(used)})")
    print(f"[05.01a] provenance: {prov_f}")


if __name__ == "__main__":
    try:
        main()
    except CovariateLockError as e:
        raise SystemExit(f"COVARIATE LOCK: {e}")
