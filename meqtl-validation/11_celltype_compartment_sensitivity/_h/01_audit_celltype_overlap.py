#!/usr/bin/env python3
"""Audit MuSiC / DNAm scMD overlap with AA meQTL inclusion lists."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
OUT = PROJECT / "meqtl-validation/11_celltype_compartment_sensitivity/_m"
REGIONS = ["caudate", "dlpfc", "hippocampus"]


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for region in REGIONS:
        cov = pd.read_csv(
            PROJECT / f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/prepared/covariates.txt",
            sep="\t",
        )
        ids = set(cov["id"].astype(str))
        music = pd.read_csv(
            PROJECT / f"inputs/cell_proportions/_m/music-proportions-{region}.tsv", sep="\t"
        )
        music_ids = set(music["sample_id"].astype(str))
        dnam_path = PROJECT / f"inputs/cell_proportions/_m/dnam-scmd-proportions-{region}.tsv"
        dnam = pd.read_csv(dnam_path, sep="\t") if dnam_path.exists() else pd.DataFrame()
        dnam_ids = set(dnam["sample_id"].astype(str)) if len(dnam) else set()

        val_path = PROJECT / f"inputs/cell_proportions/_m/dnam-scmd-validation-{region}.tsv"
        gate = False
        gate_note = "missing_validation"
        if val_path.exists():
            v = pd.read_csv(val_path, sep="\t")
            gate = len(v) == 1 and truthy(v.iloc[0].get("integration_pass", False))
            gate_note = str(v.iloc[0].to_dict()) if len(v) else "empty"

        m5 = PROJECT / (
            f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/"
            "covariate_sensitivity/covariates_M5.txt"
        )
        m6d = PROJECT / (
            f"meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/"
            "covariate_sensitivity/covariates_M6d.txt"
        )
        rows.append(
            {
                "region": region,
                "n_meqtl_aa": len(ids),
                "n_music_total": len(music_ids),
                "n_music_overlap": len(ids & music_ids),
                "frac_music_overlap": len(ids & music_ids) / len(ids) if ids else float("nan"),
                "n_dnam_total": len(dnam_ids),
                "n_dnam_overlap": len(ids & dnam_ids),
                "frac_dnam_overlap": len(ids & dnam_ids) / len(ids) if ids else float("nan"),
                "dnam_integration_pass": gate,
                "dnam_validation_note": gate_note[:200],
                "m5_covariates_exist": m5.exists(),
                "m6d_covariates_exist": m6d.exists(),
                "music_cell_types": ",".join(sorted(music["cell_type"].astype(str).unique())),
                "dnam_cell_types": (
                    ",".join(sorted(dnam["cell_type"].astype(str).unique())) if len(dnam) else ""
                ),
            }
        )
    write_tsv(OUT / "celltype_sample_overlap.tsv", rows)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
