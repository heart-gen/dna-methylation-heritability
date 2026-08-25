#!/usr/bin/env python3
"""Cis-meQTL mapping for one autosome (05_cpg_meqtl_burden).

Usage (inside an array task):
    python _h/02_map_cpg_meqtl.py --run-id cmb-AA-caudate-20260823 --chrom 22

Runs both passes that config/meqtl_parameters.yml prescribes as
`mode_primary: cis_nominal_and_permutation`:

  * a permutation pass, giving each CpG one beta-approximated p-value and its
    lead SNP -- the CpG-level evidence the VMR burden fraction counts;
  * a nominal pass, giving every tested CpG x SNP pair, which the QC stage needs
    for the p-value histogram, QQ plot and genomic-inflation estimate.

This stage does NOT compute q-values. `fdr_family: per_brain_region` means the
FDR family spans all autosomes of a region, so calling Storey's method inside a
per-chromosome array task would produce 22 unrelated families and a burden
fraction that depends on how the work was parallelized. 02b_combine_meqtl.R
pools the chromosomes and applies the correction once.

The core mapping helpers are carried over from the validated legacy
implementation at meqtl-validation/01_cpg_meqtl_mapping/_h/04_tensorqtl_map.py:
the chr-prefix normalization and the paired-chunk phenotype filtering exist
because tensorqtl's chunk pairing silently mismatches phenotypes to genotype
windows without them. Do not simplify them away.
"""

from __future__ import annotations

import argparse
import os
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


def normalize_variant_chrom(variant_df: pd.DataFrame) -> pd.DataFrame:
    """Give every contig a single 'chr' prefix, matching the phenotype BED."""
    variant_df = variant_df.copy()
    variant_df.loc[:, "chrom"] = (
        "chr" + variant_df["chrom"].astype(str).str.replace("^chr", "", regex=True)
    )
    return variant_df


def prepare_chr_matched_genotypes(pgr):
    """Align PgenReader contig names with the chr-prefixed phenotype BED."""
    pgr.variant_df = normalize_variant_chrom(pgr.variant_df)
    if "index" not in pgr.variant_df.columns:
        pgr.variant_df = pgr.variant_df.copy()
        pgr.variant_df["index"] = range(pgr.variant_df.shape[0])
    pgr.variant_dfs = {
        c: g[["pos", "index"]]
        for c, g in pgr.variant_df.groupby("chrom", sort=False)
    }
    pvar = pgr.pvar_df.copy()
    pvar["chrom"] = "chr" + pvar["chrom"].astype(str).str.replace("^chr", "", regex=True)
    pgr.pvar_df = pvar
    return pgr


def filter_phenotypes_with_cis_variants(genotypeio, pgr, phenotype_df,
                                        phenotype_pos_df, window):
    """Drop phenotypes lacking cis variants so chunk row indices stay aligned.

    tensorqtl.generate_paired_chunks builds range_df only from phenotypes with
    cis variants, then selects phenotype_df[slice] POSITIONALLY. If any
    phenotypes were dropped from range_df, later chunks pair the wrong
    phenotypes with genotype windows.
    """
    cis_ranges, drop_ids = genotypeio.get_cis_ranges(
        phenotype_pos_df, pgr.variant_dfs, window, verbose=True
    )
    if drop_ids:
        print(f"Pre-dropping {len(drop_ids)} CpGs with no cis-window variant")
        phenotype_df = phenotype_df.drop(index=drop_ids)
        phenotype_pos_df = phenotype_pos_df.drop(index=drop_ids)
    ordered = list(cis_ranges.keys())
    return phenotype_df.loc[ordered], phenotype_pos_df.loc[ordered], list(drop_ids)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--threads", default="1")
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    args = ap.parse_args()

    if args.device == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    root = repo_root()
    run_dir = root / "05_cpg_meqtl_burden" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")

    cfg = yaml.safe_load((root / "config" / "meqtl_parameters.yml").read_text())
    window = int(cfg["cis_window_bp"])
    maf = float(cfg["genotype_qc"]["maf_min"])
    seed = int(cfg["mapping"]["seed"])

    chrom_label = f"chr{args.chrom}"
    in_dir = run_dir / "inputs"
    out_dir = run_dir / "results" / "meqtl"
    out_dir.mkdir(parents=True, exist_ok=True)

    if (in_dir / f"{chrom_label}.no-tested-cpgs").exists():
        (out_dir / f"{chrom_label}.skipped").write_text("no tested CpG\n")
        print(f"{chrom_label}: no tested CpGs; skipping")
        return

    bed_f = in_dir / f"{chrom_label}.phenotype.bed.gz"
    cov_f = in_dir / f"{chrom_label}.covariates.tsv"
    geno_prefix = str(in_dir / chrom_label)
    for f in (bed_f, cov_f, Path(geno_prefix + ".pgen")):
        if not Path(f).exists():
            raise SystemExit(f"Missing prepared input {f}; run 01b first")

    import tensorqtl
    from tensorqtl import cis, genotypeio, pgen

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(str(bed_f))
    phenotype_df.columns = phenotype_df.columns.astype(str)

    probe = pgen.PgenReader(geno_prefix)
    geno_ids = set(map(str, probe.sample_ids))
    sample_ids = [s for s in phenotype_df.columns if s in geno_ids]
    if not sample_ids:
        raise SystemExit(
            "No shared donors between phenotypes and genotypes. "
            f"phenotype example={list(phenotype_df.columns[:3])}; "
            f"genotype example={list(probe.sample_ids[:3])}")

    cov = pd.read_csv(cov_f, sep="\t", index_col=0).T
    cov.index = cov.index.astype(str)
    missing = [s for s in sample_ids if s not in cov.index]
    if missing:
        raise SystemExit(f"Covariates missing {len(missing)} donors, e.g. {missing[:5]}")
    covariates_df = cov.loc[sample_ids].apply(pd.to_numeric, errors="raise")
    phenotype_df = phenotype_df[sample_ids]
    assert covariates_df.index.equals(phenotype_df.columns)

    pgr = pgen.PgenReader(geno_prefix, select_samples=sample_ids)
    pgr = prepare_chr_matched_genotypes(pgr)

    n_prepared = phenotype_df.shape[0]
    phenotype_df, phenotype_pos_df, dropped = filter_phenotypes_with_cis_variants(
        genotypeio, pgr, phenotype_df, phenotype_pos_df, window)

    print(f"{chrom_label}: donors={len(sample_ids)} prepared_cpgs={n_prepared} "
          f"tested_cpgs={phenotype_df.shape[0]} variants={pgr.num_variants}")

    # AGENTS.md 7.5: "report tested CpGs separately from prepared but untested
    # CpGs". A CpG with no cis variant was prepared and is NOT testable; it must
    # never silently leave the denominator.
    pd.DataFrame({
        "cpg_id": dropped,
        "reason": "no_variant_in_cis_window",
    }).to_csv(out_dir / f"{chrom_label}.untested-cpgs.tsv", sep="\t", index=False)

    genotype_df = pgr.load_genotypes()
    variant_df = normalize_variant_chrom(pgr.variant_df.copy())

    # ------------------------------------------------------- permutation pass
    cis_df = cis.map_cis(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        covariates_df=covariates_df, maf_threshold=maf, window=window, seed=seed,
    )
    cis_df.index.name = "cpg_id"
    cis_df.to_csv(out_dir / f"{chrom_label}.cis_qtl.tsv.gz", sep="\t",
                  float_format="%.6g")

    # ----------------------------------------------------------- nominal pass
    nominal_dir = out_dir / "nominal"
    nominal_dir.mkdir(exist_ok=True)
    cis.map_nominal(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        prefix=chrom_label, covariates_df=covariates_df,
        maf_threshold=maf, window=window, output_dir=str(nominal_dir),
    )

    print(f"{chrom_label}: wrote permutation and nominal results to {out_dir}")


if __name__ == "__main__":
    main()
