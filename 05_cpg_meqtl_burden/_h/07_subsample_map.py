#!/usr/bin/env python3
"""Re-map cis-meQTLs on paired delete-d donor draws (05_cpg_meqtl_burden).

Usage (inside an array task):
    python _h/07_subsample_map.py --run-id cmb-AA-crossregion-YYYYMMDD \
        --chrom 22 --draw-start 1 --draw-end 10

For each contrast region, loads that accepted cell's PREPARED inputs for one
autosome (phenotype BED, the executed M3a covariates, genotypes) and, for each
draw, runs the permutation pass of 02_map_cpg_meqtl.py on the donors the draw
leaves in. Only the per-CpG permutation p-value is kept: that is the
CpG-level evidence the burden counts. No q-value is formed here, for the reason
02_map_cpg_meqtl.py gives -- the FDR family is the whole region, so it is
applied once per draw by 08a_subsample_counts.py.

What is held fixed, and why. The covariate values (genotype PCs, methylation
PCs, age, sex, diagnosis) are each donor's values in the accepted cell. The
methylation PCs are therefore not re-estimated per draw. Re-estimating them
would also re-estimate the outcome's adjustment, which is the right thing for a
full end-to-end bootstrap and is not what this one claims to be: it measures the
donor-sampling variance of the slope GIVEN the locked design, as 09b's bootstrap
does. Module 02's score is likewise held fixed (it is the predictor, and the
burden model conditions on it).

No donor is ever duplicated. 06_slope_inference_new_run.R says why a bootstrap
was abandoned; the guard below refuses a draw that would map a donor twice.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

H_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(H_DIR))
from meqtl_covariates import read_manifest, repo_root  # noqa: E402

# The chunk-alignment helpers live in 02_map_cpg_meqtl.py and must not be
# re-implemented: their docstrings say why. A leading digit makes it
# unimportable by name, so load it by path.
_spec = importlib.util.spec_from_file_location("map02", H_DIR / "02_map_cpg_meqtl.py")
map02 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(map02)


def load_cell(cell_dir: Path, chrom_label: str, window: int):
    import tensorqtl
    from tensorqtl import genotypeio, pgen

    in_dir = cell_dir / "inputs"
    if (in_dir / f"{chrom_label}.no-tested-cpgs").exists():
        return None
    ph, pos = tensorqtl.read_phenotype_bed(str(in_dir / f"{chrom_label}.phenotype.bed.gz"))
    ph.columns = ph.columns.astype(str)
    cov = pd.read_csv(in_dir / f"{chrom_label}.covariates.tsv", sep="\t", index_col=0).T
    cov.index = cov.index.astype(str)
    geno_prefix = str(in_dir / chrom_label)
    probe = pgen.PgenReader(geno_prefix)
    gids = set(map(str, probe.sample_ids))
    ids = [s for s in ph.columns if s in gids]
    cov = cov.loc[ids].apply(pd.to_numeric, errors="raise")
    ph = ph[ids]
    pgr = map02.prepare_chr_matched_genotypes(pgen.PgenReader(geno_prefix, select_samples=ids))
    ph, pos, _ = map02.filter_phenotypes_with_cis_variants(genotypeio, pgr, ph, pos, window)
    G = pgr.load_genotypes()
    V = map02.normalize_variant_chrom(pgr.variant_df.copy())
    return {"ph": ph, "pos": pos, "cov": cov, "G": G, "V": V, "ids": ids}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--draw-start", type=int, required=True)
    ap.add_argument("--draw-end", type=int, required=True)
    args = ap.parse_args()
    os.environ["CUDA_VISIBLE_DEVICES"] = ""

    root = repo_root()
    mod = root / "05_cpg_meqtl_burden" / "_m" / "runs"
    run_dir = mod / args.run_id
    man = read_manifest(run_dir)
    cfg = yaml.safe_load((root / "config" / "meqtl_parameters.yml").read_text())
    si = cfg["cross_region_slope_inference"]
    window = int(cfg["cis_window_bp"])
    maf = float(cfg["genotype_qc"]["maf_min"])
    nperm = int(si["draw_nperm"])
    seed0 = int(cfg["mapping"]["seed"])
    contrast = man["contrast"].split(",")

    from tensorqtl import cis

    deleted = pd.read_csv(run_dir / "inputs" / "subsample-deletions.tsv", sep="\t",
                          dtype={"donor": str})
    chrom_label = f"chr{args.chrom}"
    for region in contrast:
        cell_dir = mod / man[f"upstream_cpg_meqtl_burden_{region}"]
        out_dir = run_dir / "results" / "subsample" / region
        out_dir.mkdir(parents=True, exist_ok=True)
        out_f = out_dir / f"{chrom_label}.draws{args.draw_start:04d}-{args.draw_end:04d}.tsv.gz"
        if out_f.exists():
            print(f"{out_f.name} exists; skipping")
            continue
        cell = load_cell(cell_dir, chrom_label, window)
        if cell is None:
            pd.DataFrame(columns=["draw", "cpg_id", "pval_beta", "pval_perm"]).to_csv(
                out_f, sep="\t", index=False)
            continue
        parts = []
        for b in range(args.draw_start, args.draw_end + 1):
            drop = set(deleted.loc[deleted["draw"] == b, "donor"])
            if not drop:
                raise SystemExit(f"draw {b} has no deletions; the draw table is incomplete")
            picked = [d for d in cell["ids"] if d not in drop]
            if len(set(picked)) != len(picked):
                raise SystemExit(f"{region} draw {b}: a donor would be mapped twice")
            G = cell["G"][picked]
            ph = cell["ph"][picked]
            cov = cell["cov"].loc[picked]
            # The locked seed offset by draw: each draw gets its own permutation
            # stream, and a rerun of a draw reproduces it exactly.
            res = cis.map_cis(G, cell["V"], ph, cell["pos"], covariates_df=cov,
                              maf_threshold=maf, window=window, nperm=nperm,
                              seed=seed0 + b)
            res.index.name = "cpg_id"
            res = res.reset_index()[["cpg_id", "pval_beta", "pval_perm"]]
            res.insert(0, "draw", b)
            res["n_samples"] = len(picked)
            parts.append(res)
            print(f"{region} {chrom_label} draw {b}: {len(picked)} samples, "
                  f"{res.shape[0]} CpGs", flush=True)
        tmp = out_f.with_suffix(".partial")
        pd.concat(parts).to_csv(tmp, sep="\t", index=False, float_format="%.6g",
                                compression="gzip")
        tmp.rename(out_f)


if __name__ == "__main__":
    main()
