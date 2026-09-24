#!/usr/bin/env python3
"""The Module 05 covariate matrix, assembled once and used by two stages.

`01a_estimate_latent_factors.py` needs the M0 design to residualize methylation
on; `01b_prepare_meqtl_inputs.py` needs the full locked design to hand to
tensorqtl. Those are the same assembly with a different term list, so it lives
here rather than in two files that can drift -- which is how the module came to
fit `snpPC1-3` while `config/covariates.yml` locked `snpPC1-5 + methPC1-5`.

Nothing here decides which terms are in the design. That comes from
`00_shared/covariate_lock.py`, i.e. from `config/covariates.yml`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

def repo_root(start: Path | None = None) -> Path:
    d = (start or Path(__file__)).resolve()
    while d != d.parent:
        if (d / ".git").is_dir() or (d / ".git").is_file():
            return d
        d = d.parent
    raise SystemExit("Could not locate repository root")


## 00_shared is located by walking up to the repository root, NOT by a path
## relative to this file. submit_meqtl_burden.sh snapshots _h/ into
## {run_dir}/code/_h and every stage executes from the snapshot, where a
## "../../00_shared" would resolve to {run_dir}/00_shared and not exist. Walking
## up finds the real root from either location, and it matches what the R stages
## already do through V2_REPO_ROOT.
sys.path.insert(0, str(repo_root() / "00_shared"))
from covariate_lock import (  # noqa: E402
    CovariateLockError, LockedDesign, assert_executed_design,
    latent_factor_recipe, load_covariates_config, locked_meqtl_design,
)

__all__ = [
    "CovariateLockError", "LockedDesign", "assert_executed_design",
    "build_covariate_matrix", "genotype_arm", "latent_factor_recipe",
    "load_covariates_config", "locked_meqtl_design", "read_manifest",
    "repo_root", "run_directory", "vmr_catalog_dir",
]


def run_directory(root: Path, run_id: str) -> Path:
    run_dir = Path(root) / "05_cpg_meqtl_burden" / "_m" / "runs" / run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")
    return run_dir


def read_manifest(run_dir: Path) -> dict:
    m = pd.read_csv(Path(run_dir) / "manifest.tsv", sep="\t", dtype=str)
    return dict(zip(m["field"], m["value"]))


def genotype_arm(root: Path, cohort: str) -> dict:
    gen = yaml.safe_load((Path(root) / "config" / "paths.yml").read_text())["genotype"]
    return gen["AA"] if cohort == "AA" else gen["all_individuals"]


def vmr_catalog_dir(root: Path, man: dict) -> Path:
    return (Path(root) / "01_vmr_catalog" / "_m" / "runs"
            / man["upstream_vmr_catalog_run_id"])


def _base_covariates(cat_dir: Path, cohort: str, chrom: str) -> pd.DataFrame:
    """agedeath / sex / primarydx from the SAME Module 01 files locus_io.R reads."""
    cov_dir = Path(cat_dir) / "covs" / f"chr_{chrom}"
    prefix = "TOPMed_LIBD.AA" if cohort == "AA" else "TOPMed_LIBD"
    covar = pd.read_csv(cov_dir / f"{prefix}.covar", sep=r"\s+", header=None,
                        names=["FID", "IID", "sex", "diagnosis"],
                        dtype={0: str, 1: str})
    qcovar = pd.read_csv(cov_dir / f"{prefix}.qcovar", sep=r"\s+", header=None,
                         names=["FID", "IID", "age"], dtype={0: str, 1: str})
    return covar.merge(qcovar, on=["FID", "IID"])


def _snp_pcs(root: Path, cohort: str, pcs: tuple[str, ...]) -> pd.DataFrame:
    """The ancestry PCs the lock names -- all of them, or a hard failure.

    The previous draft printed a warning and continued without genotype PCs when
    the eigenvec was missing. A silently ancestry-unadjusted meQTL scan in an
    admixed cohort is not a degraded run, it is a different analysis, so this is
    now fatal.
    """
    arm = genotype_arm(root, cohort)
    evec_f = Path(root) / arm["eigenvec"]
    if not evec_f.is_file():
        raise CovariateLockError(
            f"No eigenvec at {evec_f}, so the locked ancestry PCs {list(pcs)} "
            "cannot be fitted. Refusing to map an ancestry-unadjusted scan.")
    evec = pd.read_csv(evec_f, sep=r"\s+", dtype={0: str, 1: str})
    evec.columns = ["FID", "IID"] + [f"snpPC{i}"
                                     for i in range(1, evec.shape[1] - 1)]
    have = [p for p in pcs if p in evec.columns]
    if len(have) != len(pcs):
        raise CovariateLockError(
            f"{evec_f} carries {evec.shape[1] - 2} PCs; the lock names "
            f"{list(pcs)} and {sorted(set(pcs) - set(have))} are absent.")
    return evec[["FID", "IID"] + list(pcs)]


def build_covariate_matrix(root: Path, cat_dir: Path,
                           cohort: str, chrom: str, donors: list[str],
                           design: LockedDesign, terms: str = "full",
                           latent_file: Path | None = None,
                           context: str = "", allow_unlocked: bool = False,
                           ) -> tuple[pd.DataFrame, dict]:
    """The design matrix for one chromosome, donors x covariates.

    `terms="full"` is the locked model. `terms="m0"` is the same thing without the
    latent factors, which is what the factors themselves are estimated against --
    `primary_meqtl.sensitivity_models` defines M0 as
    `agedeath+sex+primarydx+snpPC1-5` and M3a as `M0 + methPC1-5`, so M0 carries
    EVERY locked ancestry PC. The lock does not say that; the 2026-09-23 pilot
    resolved it and it is recorded per run.

    Returns the matrix (indexed by donor, in `donors` order) and the record
    `assert_executed_design()` produced.
    """
    if terms not in ("full", "m0"):
        raise ValueError(f"terms must be 'full' or 'm0', got {terms!r}")

    cov = _base_covariates(cat_dir, cohort, chrom)
    cov = cov.merge(_snp_pcs(root, cohort, design.ancestry_pcs),
                    on=["FID", "IID"], how="left")
    cov = cov.set_index("FID").drop(columns=["IID"])
    missing_donors = [d for d in donors if d not in cov.index]
    if missing_donors:
        raise CovariateLockError(
            f"{len(missing_donors)} analysis donor(s) have no covariate row, "
            f"e.g. {missing_donors[:5]}. Donor alignment is by ID and a missing "
            "ID fails loudly (AGENTS.md 7.1, 10.1).")
    cov = cov.loc[list(donors)]

    if terms == "full" and design.latent_factors:
        if latent_file is None or not Path(latent_file).is_file():
            raise CovariateLockError(
                f"The locked model {design.model_id} includes "
                f"{list(design.latent_factors)} but no latent-factor table was "
                f"found at {latent_file}. Run _h/01a_estimate_latent_factors.py "
                "for this run first.")
        mp = pd.read_csv(latent_file, sep="\t", index_col=0)
        mp.index = mp.index.astype(str)
        want = list(design.latent_factors)
        absent_cols = [c for c in want if c not in mp.columns]
        if absent_cols:
            raise CovariateLockError(
                f"{latent_file} lacks {absent_cols}; it carries "
                f"{list(mp.columns)}.")
        absent_rows = [d for d in cov.index if d not in mp.index]
        if absent_rows:
            raise CovariateLockError(
                f"{latent_file} is missing donors {absent_rows[:5]}; it was "
                "estimated on a different donor set than this run.")
        cov = cov.join(mp.loc[cov.index, want])

    ## tensorqtl builds its residualizer straight from the values, so the design
    ## must be numeric. `sex` and `diagnosis` arrive as strings ("F"/"M",
    ## "Control"/"Schizo"); dummy-code them dropping the first level, the same
    ## treatment coding 00_shared/locus_io.R applies in Modules 02 and 03.
    categorical = [c for c in cov.columns if cov[c].dtype == object]
    if categorical:
        cov = pd.get_dummies(cov, columns=categorical, drop_first=True,
                             dtype=float)
    cov = cov.apply(pd.to_numeric, errors="coerce")
    if cov.isnull().any().any():
        bad = cov.columns[cov.isnull().any()].tolist()
        raise CovariateLockError(f"null covariate values in {bad}")

    ## A constant column makes the residualizer rank-deficient, so it still has
    ## to go -- but under the lock, dropping one means a locked term is no longer
    ## in the design, and the assertion below is what turns that into a refusal
    ## rather than a printed note nobody reads.
    constant = [c for c in cov.columns if cov[c].nunique() < 2]
    if constant:
        print(f"dropping constant covariate(s) {constant}")
        cov = cov.drop(columns=constant)

    expected = design
    if terms == "m0":
        expected = LockedDesign(
            model_id=f"{design.model_id}:M0", lock_status=design.lock_status,
            lock_date=design.lock_date, base_terms=design.base_terms,
            ancestry_pcs=design.ancestry_pcs, latent_factors=(),
            config_sha256=design.config_sha256, config_file=design.config_file)
    record = assert_executed_design(
        cov.columns, expected, context or f"chr{chrom} ({terms})",
        allow_unlocked=allow_unlocked)
    return cov, record
