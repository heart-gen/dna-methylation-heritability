#!/usr/bin/env python3
"""Build shared-donor-aware caudate downsample sample lists (Tier A).

Primary target N = min(N_DLPFC, N_hippocampus) = 111.
Each replicate keeps all donors shared across all three regions, then fills
randomly from remaining caudate donors (without replacement).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import write_tsv  # noqa: E402

PROJECT = Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
PHASE1 = PROJECT / "meqtl-validation/01_cpg_meqtl_mapping"
P7 = PROJECT / "meqtl-validation/08_schizophrenia_risk_application/_m"
SEED = 20260801
N_REPS = 30


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-reps", type=int, default=N_REPS)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--target-n", type=int, default=0, help="0 = min(DLPFC, hippocampus)")
    p.add_argument("--outdir", default=str(P7 / "caudate_downsample"))
    return p.parse_args()


def region_samples(region: str) -> list[str]:
    pref = PHASE1 / region / "_m" / "preflight" / "sample_inclusion_primary.tsv"
    if pref.exists():
        df = pd.read_csv(pref, sep="\t")
        col = "brnum" if "brnum" in df.columns else df.columns[0]
        return sorted(df[col].astype(str).unique())
    cov = pd.read_csv(PHASE1 / region / "_m/prepared/covariates.txt", sep="\t", index_col=0)
    if cov.index.astype(str).str.startswith("Br").any():
        return sorted(cov.index.astype(str))
    return sorted(cov.columns.astype(str))


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    lists_dir = outdir / "sample_lists"
    lists_dir.mkdir(parents=True, exist_ok=True)

    caud = region_samples("caudate")
    dlpfc = region_samples("dlpfc")
    hip = region_samples("hippocampus")
    shared3 = sorted(set(caud) & set(dlpfc) & set(hip))
    caud_only = sorted(set(caud) - set(shared3))
    target_n = args.target_n if args.target_n > 0 else min(len(dlpfc), len(hip))

    if len(shared3) > target_n:
        raise SystemExit(
            f"shared3={len(shared3)} exceeds target_n={target_n}; cannot keep all shared donors"
        )
    n_fill = target_n - len(shared3)
    if n_fill > len(caud_only):
        raise SystemExit(
            f"Need {n_fill} fill donors but only {len(caud_only)} non-shared caudate donors"
        )

    rng = np.random.default_rng(args.seed)
    rows = []
    for rep in range(1, args.n_reps + 1):
        fill = sorted(rng.choice(caud_only, size=n_fill, replace=False).tolist())
        samples = sorted(shared3 + fill)
        assert len(samples) == target_n
        path = lists_dir / f"caudate_downsample_rep{rep:03d}.txt"
        path.write_text("\n".join(samples) + "\n")
        rows.append({
            "replicate": rep,
            "n_samples": len(samples),
            "n_shared3_kept": len(shared3),
            "n_fill": n_fill,
            "path": str(path),
            "seed": args.seed,
        })

    write_tsv(outdir / "downsample_replicate_manifest.tsv", rows)
    write_tsv(outdir / "downsample_design_summary.tsv", [{
        "n_caudate_full": len(caud),
        "n_dlpfc": len(dlpfc),
        "n_hippocampus": len(hip),
        "n_shared_all3": len(shared3),
        "n_caudate_nonshared": len(caud_only),
        "target_n": target_n,
        "n_reps": args.n_reps,
        "seed": args.seed,
        "strategy": "keep_all_shared3_then_fill_random_without_replacement",
        "matched_to": "min(n_dlpfc, n_hippocampus)",
    }])
    # Save shared3 list for audit
    (outdir / "shared3_donors.txt").write_text("\n".join(shared3) + "\n")
    print(
        f"Wrote {args.n_reps} lists under {lists_dir} "
        f"(target_n={target_n}; shared3={len(shared3)}; fill={n_fill})"
    )


if __name__ == "__main__":
    main()
