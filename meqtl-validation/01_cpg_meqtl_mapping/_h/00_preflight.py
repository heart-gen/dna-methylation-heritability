#!/usr/bin/env python3
"""
Phase 1 preflight: sample, CpG, genotype, and VMR readiness for cis-meQTL.

Writes inclusion lists and inventory counts. Does not run TensorQTL.

Population AA uses AA-only genotypes + AA-only CpG matrices.
Population EA uses all_individuals genotypes + all_individuals CpG matrices.
"""

from __future__ import annotations

import argparse, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import (  # noqa: E402
    cpg_matrix_relpath,
    cpg_samples_relpath,
    genotype_paths,
    load_paths,
    load_yaml,
    nonmissing,
    norm_region,
    read_psam_brnums,
    read_tsv,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True, choices=["caudate", "dlpfc", "hippocampus"])
    p.add_argument("--population", default="AA", choices=["AA", "EA"])
    return p.parse_args()


def load_cpg_ids(staff: Path, paths: dict, region: str, population: str) -> set[str]:
    # Prefer samples.txt (fast); fall back to phen FID column
    samples = staff / cpg_samples_relpath(paths, region, population)
    if samples.exists():
        ids: set[str] = set()
        with samples.open() as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if parts and parts[0] and not parts[0].startswith("#"):
                    ids.add(parts[0])
        return ids
    for chrom in ("22", "21", "1"):
        phen_path = staff / cpg_matrix_relpath(paths, region, chrom, population)
        if not phen_path.exists():
            continue
        ids = set()
        with phen_path.open() as handle:
            handle.readline()
            for line in handle:
                parts = line.rstrip("\n").split("\t", 1)
                if parts:
                    ids.add(parts[0])
        return ids
    return set()


def main() -> None:
    args = parse_args()
    paths = load_paths()
    meqtl = load_yaml("meqtl_parameters.yml")
    covars = load_yaml("covariates.yml")
    project = Path(paths["project_root"])
    staff = Path(paths["cpg_methylation_root"])
    # Keep AA primary outputs in preflight/ for backward compatibility;
    # EA (and any non-AA) under preflight/{population}/.
    if args.population == "AA":
        outdir = project / "meqtl-validation" / "01_cpg_meqtl_mapping" / args.region / "_m" / "preflight"
    else:
        outdir = (
            project / "meqtl-validation" / "01_cpg_meqtl_mapping" / args.region / "_m"
            / "preflight" / args.population
        )
    outdir.mkdir(parents=True, exist_ok=True)

    phen = read_tsv(project / paths["phenotype_table"])
    required = covars["primary_meqtl"]["required_phenotype_columns"]
    samples = []
    for row in phen:
        if row.get("race") != args.population:
            continue
        if norm_region(row.get("region", "")) != args.region:
            continue
        if not all(nonmissing(row, c) for c in required):
            continue
        samples.append(row)

    geno = genotype_paths(paths, args.population)
    geno_ids = read_psam_brnums(project / geno["psam"])
    cpg_ids = load_cpg_ids(staff, paths, args.region, args.population)

    included = []
    for row in samples:
        br = row["brnum"]
        ok_geno = br in geno_ids
        ok_cpg = br in cpg_ids
        included.append({
            "brnum": br,
            "region": args.region,
            "race": args.population,
            "primarydx": row.get("primarydx", ""),
            "agedeath": row.get("agedeath", ""),
            "sex": row.get("sex", ""),
            "has_genotype": ok_geno,
            "has_cpg_matrix": ok_cpg,
            "include_primary": ok_geno and ok_cpg,
        })

    keep = [r for r in included if r["include_primary"]]
    write_tsv(outdir / "sample_inclusion.tsv", included)
    write_tsv(
        outdir / "sample_inclusion_primary.tsv",
        [{"brnum": r["brnum"]} for r in keep],
        ["brnum"],
    )

    if args.population == "AA":
        pred_key = "local_predictability_summary_template"
    else:
        pred_key = "local_predictability_ea_template"
    vmr_bed = project / paths["vmr_bed_template"].format(region=args.region)
    pred = project / paths[pred_key].format(region=args.region)
    n_vmr = sum(1 for _ in vmr_bed.open()) if vmr_bed.exists() else 0
    n_pred = sum(1 for _ in pred.open()) - 1 if pred.exists() else 0

    cpg_rel = cpg_matrix_relpath(paths, args.region, "22", args.population)
    cpg_root = (staff / cpg_rel).parent.parent
    chr_dirs = sorted(cpg_root.glob("chr_*")) if cpg_root.exists() else []
    summary = [{
        "region": args.region,
        "population": args.population,
        "n_phenotype_complete_covariates": len(samples),
        "n_included_primary": len(keep),
        "n_vmrs": n_vmr,
        "n_predictability_rows": n_pred,
        "n_cpg_chr_dirs": len(chr_dirs),
        "cis_window_bp": meqtl["cis_window_bp"],
        "maf_min": meqtl["genotype_qc"]["maf_min"],
        "imputation_r2_min": meqtl["genotype_qc"]["imputation_r2_min"],
        "methylation_source": meqtl["cpg_qc"]["methylation_source"],
        "fdr_family": meqtl["mapping"]["fdr_family"],
        "vmr_bed": str(vmr_bed),
        "predictability_path": str(pred),
        "genotype_prefix": str(project / geno["pgen"]).replace(".pgen", ""),
        "cpg_root": str(cpg_root),
        "cpg_matrix_template": cpg_rel,
    }]
    write_tsv(outdir / "preflight_summary.tsv", summary)
    print(f"{args.region} {args.population}: {len(keep)} samples ready; wrote {outdir}")


if __name__ == "__main__":
    main()
