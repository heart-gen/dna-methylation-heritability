#!/usr/bin/env python3
"""The locked primary meQTL covariate design, read from config rather than typed.

AGENTS.md 12 makes the covariate model a PI decision and 7.5 requires Module 05
to "preserve the locked donor and covariate models".
`config/covariates.yml:primary_meqtl` records that decision:

    lock_status: locked_M3a          lock_date: 2026-08-01
    M3a = agedeath + sex + primarydx + snpPC1-5 + methPC1-5

Until 2026-09-24 nothing in Module 05 read that file. `01b_prepare_meqtl_inputs.py`
set `n_pc = 3` as a literal and added no methylation PC, so the three accepted
runs `cmb-AA-{caudate,dlpfc,hippocampus}-20260825` were mapped under

    executed = agedeath + sex + primarydx + snpPC1-3

while their manifests checksummed the lock they did not follow. The PI decision
of 2026-09-24 is that the lock is authoritative, so this module exists to make
the config the single source of the design: `locked_meqtl_design()` expands it,
and `assert_executed_design()` refuses a matrix that does not match.

The assertion is deliberately applied to the columns of the covariate matrix
that is actually handed to tensorqtl, not to a config value echoed back to
itself. A run that reads the lock, builds something else, and then reports the
lock is exactly the failure this file is here to make impossible.

The latent-factor half of the lock is UNDERSPECIFIED, and that is recorded
rather than papered over: `latent_factor_policy` names a method ("PCA on
M0-residualized CpG phenotypes") and no file, CpG set, subsample, seed,
standardization, or answer to whether the M0 residualization carries snpPC1-5.
`latent_factor_recipe()` resolves every one of those, prefers an explicit
`primary_meqtl.latent_factor_recipe` block, which the PI added on 2026-09-24
(see 05_cpg_meqtl_burden/PROPOSED_CONFIG_CHANGE.md for what was decided and
why), and reports which source it used so the run manifest
records the resolution instead of implying the config supplied it.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# The latent-factor recipe the 2026-09-23 decision pilot resolved
# (05_cpg_meqtl_burden/PILOT_COVARIATE_MODEL.md, "What the methylation PCs
# actually are"). These are the v1 implementation's own choices, taken from
# meqtl-validation/01_cpg_meqtl_mapping/_h/09_estimate_latent_factors.py, so the
# recipe is the locked one rather than a newly invented one. The v1 factor TABLES
# are deliberately not reused: they were estimated on the v1 VMR catalog's CpGs
# and the v2 catalog is a different CpG set.
#
# `cpg_order` is the one choice the pilot did not have to make explicit and this
# file does. The subsample is drawn by index, so which CpGs it selects depends on
# the row order of the pooled matrix; the pilot's order was
# `sorted(glob("chr*.phenotype.bed.gz"))`, i.e. chr1, chr10, chr11, ..., chr2,
# which is not reproducible from a description. AGENTS.md 10.1 requires numeric
# chromosome ordering, so that is what is fixed here and recorded per run.
#
# `pca_solver` and `sign_convention` are likewise made explicit rather than left
# to a library default. The pilot used `sklearn.decomposition.PCA(n_components=15)`,
# whose `svd_solver="auto"` heuristic picks the RANDOMIZED solver at these shapes
# (153 x 50,000) -- reproducible only against one sklearn version, and a solver
# switch inside a minor release would silently change a production covariate.
LATENT_FACTOR_DEFAULTS: dict[str, object] = {
    "method": "pca_on_m0_residuals",
    "implementation": "v1 09_estimate_latent_factors.py, re-fitted on v2 tested CpGs",
    "seed": 20260730,
    "n_cpg_subsample": 50000,
    "n_factors_estimated": 15,
    "residualize_on": "M0",              # base terms + EVERY ancestry PC
    "cpg_standardization": "zscore_per_cpg",
    "cpg_order": "chrom_numeric_then_position",
    "missing_imputation": "cpg_mean",
    "pca_solver": "numpy_svd_full",
    "sign_convention": "max_abs_loading_positive",
}

_RECIPE_KEYS = tuple(LATENT_FACTOR_DEFAULTS)


class CovariateLockError(RuntimeError):
    """A design that does not match the locked model. Never caught internally."""


@dataclass(frozen=True)
class LockedDesign:
    model_id: str
    lock_status: str
    lock_date: str
    base_terms: tuple[str, ...]
    ancestry_pcs: tuple[str, ...]
    latent_factors: tuple[str, ...]
    config_sha256: str
    config_file: str

    @property
    def terms(self) -> tuple[str, ...]:
        return self.base_terms + self.ancestry_pcs + self.latent_factors

    @property
    def n_snp_pcs(self) -> int:
        return len(self.ancestry_pcs)

    @property
    def n_latent_factors(self) -> int:
        return len(self.latent_factors)

    def describe(self) -> str:
        return f"{self.model_id} = " + " + ".join(self.terms)


@dataclass(frozen=True)
class LatentFactorRecipe:
    source: str                      # "config" or "code_default"
    values: dict = field(default_factory=dict)

    def __getitem__(self, k):
        return self.values[k]

    @property
    def seed(self) -> int:
        return int(self.values["seed"])

    @property
    def n_cpg_subsample(self) -> int:
        return int(self.values["n_cpg_subsample"])

    @property
    def n_factors_estimated(self) -> int:
        return int(self.values["n_factors_estimated"])


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_covariates_config(root: Path) -> dict:
    f = Path(root) / "config" / "covariates.yml"
    if not f.is_file():
        raise CovariateLockError(f"Config not found: {f}")
    cfg = yaml.safe_load(f.read_text())
    cfg["_config_file"] = str(f)
    cfg["_config_sha256"] = file_sha256(f)
    return cfg


def locked_meqtl_design(cfg: dict) -> LockedDesign:
    """Expand `primary_meqtl` into the term list the mapping stage must fit.

    Every element comes from the config. `base_terms` is derived as
    `required_phenotype_columns` minus `ancestry_pcs` rather than typed out, so
    adding a phenotype covariate to the lock cannot leave this file behind.
    """
    pm = cfg.get("primary_meqtl")
    if not isinstance(pm, dict):
        raise CovariateLockError("config/covariates.yml has no primary_meqtl block")

    required = [str(x) for x in (pm.get("required_phenotype_columns") or [])]
    ancestry = [str(x) for x in (pm.get("ancestry_pcs") or [])]
    latent = [str(x) for x in (pm.get("locked_latent_factors") or [])]
    if not required:
        raise CovariateLockError(
            "primary_meqtl.required_phenotype_columns is empty; the locked design "
            "cannot be expanded and Module 05 will not guess it (AGENTS.md 12).")
    if not ancestry:
        raise CovariateLockError("primary_meqtl.ancestry_pcs is empty")

    base = [c for c in required if c not in set(ancestry)]
    missing_from_required = [pc for pc in ancestry if pc not in set(required)]
    if missing_from_required:
        raise CovariateLockError(
            "primary_meqtl.ancestry_pcs names PCs absent from "
            f"required_phenotype_columns: {missing_from_required}. The two lists "
            "disagree about the locked design, which is a PI matter, not a "
            "default this code may pick.")

    model_id = str(pm.get("locked_model") or "")
    if not model_id:
        raise CovariateLockError("primary_meqtl.locked_model is unset")

    ## The named model and the expanded term list must agree. M3a is defined in
    ## `sensitivity_models` as "M0 + methPC1-5", so a lock that says M3a while
    ## listing no latent factor is self-contradictory and is not silently taken
    ## to mean M0.
    if model_id == "M3a" and not latent:
        raise CovariateLockError(
            "locked_model is M3a but locked_latent_factors is empty. M3a is "
            "'M0 + methPC1-5' in primary_meqtl.sensitivity_models.")

    return LockedDesign(
        model_id=model_id,
        lock_status=str(pm.get("lock_status", "")),
        lock_date=str(pm.get("lock_date", "")),
        base_terms=tuple(base),
        ancestry_pcs=tuple(ancestry),
        latent_factors=tuple(latent),
        config_sha256=str(cfg.get("_config_sha256", "")),
        config_file=str(cfg.get("_config_file", "")),
    )


def latent_factor_recipe(cfg: dict) -> LatentFactorRecipe:
    """The methPC derivation, from config when pinned there and from code if not.

    Returning `source` matters as much as returning the values: the manifest must
    say whether the recipe was a PI-pinned config block or this file's resolution
    of an underspecified policy (AGENTS.md 9 -- a run carries its configuration,
    and a method description is not a specification).
    """
    pm = cfg.get("primary_meqtl") or {}
    block = pm.get("latent_factor_recipe")
    if block is None:
        return LatentFactorRecipe("code_default", dict(LATENT_FACTOR_DEFAULTS))
    if not isinstance(block, dict):
        raise CovariateLockError(
            "primary_meqtl.latent_factor_recipe must be a mapping")
    missing = [k for k in _RECIPE_KEYS if k not in block]
    if missing:
        raise CovariateLockError(
            "primary_meqtl.latent_factor_recipe is partial, which is worse than "
            f"absent -- it looks pinned but is not. Missing: {missing}. State "
            "every key or remove the block.")
    return LatentFactorRecipe("config", {k: block[k] for k in _RECIPE_KEYS})


# ---------------------------------------------------------------------------
# Executed column name -> locked term.
#
# The covariate matrix on disk does not use the config's names, for two reasons
# that are both older than this file: Module 01's `.qcovar` column is `age` where
# the config says `agedeath`, and `sex`/`diagnosis` are strings that get
# treatment-coded into `sex_M` / `diagnosis_Schizo` before tensorqtl sees them
# (tensorqtl builds its residualizer straight from the values, so the design must
# be numeric). Normalizing here means the gate compares designs rather than
# spellings, and an unmapped column is a failure rather than a silent pass.
_PC_RE = re.compile(r"^(snpPC|methPC)(\d+)$")
_DUMMY_PREFIX = {"sex": "sex", "diagnosis": "primarydx"}
_ALIAS = {"age": "agedeath", "agedeath": "agedeath",
          "sex": "sex", "diagnosis": "primarydx", "primarydx": "primarydx"}


def canonical_term(column: str) -> str | None:
    """The locked term an executed covariate column stands for, or None."""
    col = str(column)
    if col in _ALIAS:
        return _ALIAS[col]
    if _PC_RE.match(col):
        return col
    for prefix, term in _DUMMY_PREFIX.items():
        if col.startswith(prefix + "_"):
            return term
    return None


def design_terms(columns) -> tuple[list[str], list[str]]:
    """(locked terms represented, columns that map to nothing)."""
    terms, unmapped = [], []
    for c in columns:
        t = canonical_term(c)
        if t is None:
            unmapped.append(str(c))
        elif t not in terms:
            terms.append(t)
    return terms, unmapped


def assert_executed_design(columns, design: LockedDesign, context: str,
                           allow_unlocked: bool = False) -> dict:
    """Refuse a covariate matrix whose design is not the locked one.

    `columns` must be the columns of the matrix that will be (or was) fit, taken
    from the matrix itself. Returns a record for the run manifest.

    `allow_unlocked` exists only for smoke runs, which may legitimately have a
    constant covariate dropped in a handful of donors. It downgrades the refusal
    to a warning and says so in the returned record; it never makes the run
    citable (AGENTS.md 14).
    """
    cols = [str(c) for c in columns]
    got, unmapped = design_terms(cols)
    want = list(design.terms)
    missing = [t for t in want if t not in got]
    extra = [t for t in got if t not in want]

    record = {
        "locked_model": design.model_id,
        "locked_terms": ";".join(want),
        "executed_columns": ";".join(cols),
        "executed_terms": ";".join(got),
        "missing_terms": ";".join(missing),
        "extra_terms": ";".join(extra),
        "unmapped_columns": ";".join(unmapped),
        "config_sha256": design.config_sha256,
        "matches_lock": not (missing or extra or unmapped),
        "enforced": not allow_unlocked,
    }
    if record["matches_lock"]:
        return record

    parts = [f"{context}: executed covariate design does not match the locked "
             f"model {design.model_id} ({design.describe()})."]
    if missing:
        parts.append(f"Missing locked term(s): {missing}.")
    if extra:
        parts.append(f"Term(s) not in the lock: {extra}.")
    if unmapped:
        parts.append(f"Column(s) matching no locked term: {unmapped}.")
    parts.append(f"Executed columns were: {cols}.")
    parts.append("config/covariates.yml:primary_meqtl is the PI decision "
                 "(AGENTS.md 12); change the code, not the lock.")
    msg = " ".join(parts)
    if allow_unlocked:
        print("WARNING (smoke run, not citable): " + msg)
        return record
    raise CovariateLockError(msg)
