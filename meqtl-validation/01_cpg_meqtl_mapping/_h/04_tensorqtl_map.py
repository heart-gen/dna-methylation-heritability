#!/usr/bin/env python3
"""
TensorQTL cis-meQTL mapping for prepared CpG phenotype BEDs.

Adapted from sex_context_brain/eqtl_analysis/common/04.tensorqtl_map.py.
Primary: permutation cis pass with region-specific FDR family.
Uses PgenReader and chromosome-chunked genotype loading by default.
"""

from __future__ import annotations

import argparse, os
from pathlib import Path
import pandas as pd


def configure_device(device: str) -> None:
    if device == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""


def parse_chunk_size(value: str):
    normalized = value.strip().lower()
    if normalized in {"none", "all", "full"}:
        return None
    if normalized == "chr":
        return "chr"
    size = int(value)
    if size <= 0:
        raise argparse.ArgumentTypeError("chunk-size integer must be positive")
    return size


def normalize_variant_chrom(variant_df: pd.DataFrame) -> pd.DataFrame:
    variant_df = variant_df.copy()
    variant_df.loc[:, "chrom"] = (
        "chr" + variant_df["chrom"].astype(str).str.replace("^chr", "", regex=True)
    )
    return variant_df


def load_covariates(path: str, sample_order: list[str]) -> pd.DataFrame:
    """Load covariates as samples x covariates in phenotype sample order."""
    cov = pd.read_csv(path, sep="\t", index_col=0)
    # Detect orientation: TensorQTL needs samples as index.
    sample_set = set(sample_order)
    if set(cov.index.astype(str)).intersection(sample_set):
        cov.index = cov.index.astype(str)
    elif set(cov.columns.astype(str)).intersection(sample_set):
        cov = cov.T
        cov.index = cov.index.astype(str)
    else:
        raise SystemExit(
            "No overlapping sample IDs between covariates and phenotypes. "
            f"covariate index example={list(cov.index[:3])}; "
            f"phenotype example={sample_order[:3]}"
        )
    missing = [s for s in sample_order if s not in cov.index]
    if missing:
        raise SystemExit(f"Covariates missing {len(missing)} phenotype samples, e.g. {missing[:5]}")
    cov = cov.loc[sample_order]
    if cov.isnull().any().any():
        bad = cov.columns[cov.isnull().any()].tolist()
        raise SystemExit(f"Missing/null covariate values in columns: {bad}")
    return cov


def calculate_qvalues(res_df: pd.DataFrame, fdr: float = 0.05) -> None:
    from py_qvalue import qvalue
    from scipy.stats import beta

    if "pval_beta" in res_df and res_df["pval_beta"].notnull().any():
        pval_col = "pval_beta"
    else:
        pval_col = "pval_perm"
    qval_res = qvalue(res_df[pval_col])
    res_df["qval"] = qval_res["qvalues"]
    if pval_col == "pval_beta":
        sig = res_df.loc[res_df["qval"] <= fdr, "pval_beta"]
        nonsig = res_df.loc[res_df["qval"] > fdr, "pval_beta"]
        if not sig.empty:
            lb = sig.max()
            ub = nonsig.min() if not nonsig.empty else lb
            pthreshold = (lb + ub) / 2 if ub != lb else lb
            res_df["pval_nominal_threshold"] = beta.ppf(
                pthreshold, res_df["beta_shape1"], res_df["beta_shape2"]
            )


def prepare_chr_matched_genotypes(pgr):
    """Align PgenReader contig names with chr-prefixed phenotype BEDs."""
    pgr.variant_df = normalize_variant_chrom(pgr.variant_df)
    if "index" not in pgr.variant_df.columns:
        pgr.variant_df = pgr.variant_df.copy()
        pgr.variant_df["index"] = range(pgr.variant_df.shape[0])
    pgr.variant_dfs = {
        c: g[["pos", "index"]]
        for c, g in pgr.variant_df.groupby("chrom", sort=False)
    }
    # generate_paired_chunks reads chrom from pvar_df; keep that in sync too
    pvar = pgr.pvar_df.copy()
    pvar["chrom"] = "chr" + pvar["chrom"].astype(str).str.replace("^chr", "", regex=True)
    pgr.pvar_df = pvar
    return pgr


def filter_phenotypes_with_cis_variants(genotypeio, pgr, phenotype_df, phenotype_pos_df, window):
    """Drop phenotypes lacking cis variants so chunk row indices stay aligned.

    tensorqtl.generate_paired_chunks builds range_df only from phenotypes with
    cis variants, but then selects phenotype_df[slice] positionally. If any
    phenotypes were dropped from range_df, later chunks pair the wrong
    phenotypes with genotype windows and can empty out with
    'No phenotypes remain after filters.'
    """
    cis_ranges, drop_ids = genotypeio.get_cis_ranges(
        phenotype_pos_df, pgr.variant_dfs, window, verbose=True
    )
    if drop_ids:
        print(f"Pre-dropping {len(drop_ids)} phenotypes without cis-window variants")
        phenotype_df = phenotype_df.drop(index=drop_ids)
        phenotype_pos_df = phenotype_pos_df.drop(index=drop_ids)
    # Keep phenotype row order identical to range_df / cis_ranges insertion order
    ordered = list(cis_ranges.keys())
    phenotype_df = phenotype_df.loc[ordered]
    phenotype_pos_df = phenotype_pos_df.loc[ordered]
    return phenotype_df, phenotype_pos_df


def map_cis_chunked(cis, genotypeio, pgr, phenotype_df, phenotype_pos_df,
                    covariates_df, maf, window, seed, chunk_size):
    if chunk_size is None:
        genotype_df = pgr.load_genotypes()
        variant_df = normalize_variant_chrom(pgr.variant_df.copy())
        return cis.map_cis(
            genotype_df,
            variant_df,
            phenotype_df,
            phenotype_pos_df,
            covariates_df=covariates_df,
            maf_threshold=maf,
            window=window,
            seed=seed,
        )

    phenotype_df, phenotype_pos_df = filter_phenotypes_with_cis_variants(
        genotypeio, pgr, phenotype_df, phenotype_pos_df, window
    )
    print(f"Phenotypes entering chunked cis map: {phenotype_df.shape[0]}")

    chunks = []
    skipped = 0
    for gt_df, var_df, p_df, p_pos_df, chunk_idx in genotypeio.generate_paired_chunks(
        pgr, phenotype_df, phenotype_pos_df, chunk_size, window=window, verbose=True
    ):
        var_df = normalize_variant_chrom(var_df)
        # Hard assert: phenotype chroms must be present in this genotype chunk
        pheno_chrs = set(p_pos_df["chr"].unique())
        geno_chrs = set(var_df["chrom"].unique())
        if not pheno_chrs.intersection(geno_chrs):
            skipped += 1
            print(
                f"WARNING: skipping chunk {chunk_idx + 1}: phenotype chroms {sorted(pheno_chrs)} "
                f"disjoint from genotype chroms {sorted(geno_chrs)}"
            )
            continue
        print(
            f"cis.map_cis chunk {chunk_idx + 1}: phenotypes={p_df.shape[0]} "
            f"variants={gt_df.shape[0]} chroms={sorted(pheno_chrs)}"
        )
        try:
            chunks.append(
                cis.map_cis(
                    gt_df,
                    var_df,
                    p_df,
                    p_pos_df,
                    covariates_df=covariates_df,
                    maf_threshold=maf,
                    window=window,
                    seed=seed,
                )
            )
        except ValueError as exc:
            if "No phenotypes remain after filters" in str(exc):
                skipped += 1
                print(f"WARNING: skipping empty chunk {chunk_idx + 1}: {exc}")
                continue
            raise
    if not chunks:
        raise SystemExit("No cis-QTL chunks produced results")
    print(f"Finished chunked cis map: kept={len(chunks)} skipped={skipped}")
    return pd.concat(chunks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", required=True)
    parser.add_argument("--mode", choices=["cis", "nominal"], default="cis")
    parser.add_argument("--phenotype-bed", required=True, help="BED or bgzipped BED")
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--genotype-prefix", required=True, help="plink2 pfile prefix")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="cpg_meqtl")
    parser.add_argument("--window", type=int, default=500000)
    parser.add_argument("--maf", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260722)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    parser.add_argument("--sample-list", default="",
                        help="Optional one-ID-per-line sample subset applied to phenotype and genotype")
    parser.add_argument(
        "--chunk-size",
        type=parse_chunk_size,
        default="chr",
        help="Genotype chunking: chr (default), integer, or none for full load",
    )
    args = parser.parse_args()

    configure_device(args.device)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    import tensorqtl
    from tensorqtl import cis, genotypeio, pgen

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.phenotype_bed)
    phenotype_df.columns = phenotype_df.columns.astype(str)

    sample_probe = pgen.PgenReader(args.genotype_prefix)
    geno_ids = set(map(str, sample_probe.sample_ids))
    pheno_ids = [s for s in phenotype_df.columns if s in geno_ids]
    if args.sample_list:
        requested = [x.strip() for x in Path(args.sample_list).read_text().splitlines() if x.strip()]
        requested_set = set(requested)
        pheno_ids = [s for s in pheno_ids if s in requested_set]
        absent = [s for s in requested if s not in pheno_ids]
        if absent:
            raise SystemExit(f"Requested sample list contains {len(absent)} samples absent from phenotype/genotype data: {absent[:5]}")
    if not pheno_ids:
        raise SystemExit(
            "No shared samples between phenotypes and genotypes. "
            f"phenotype example={list(phenotype_df.columns[:3])}; "
            f"genotype example={list(sample_probe.sample_ids[:3])}"
        )

    covariates_df = load_covariates(args.covariates, pheno_ids)
    phenotype_df = phenotype_df[pheno_ids]
    # Keep covariate index aligned with phenotype columns (tensorqtl assert)
    assert covariates_df.index.equals(phenotype_df.columns)

    pgr = pgen.PgenReader(args.genotype_prefix, select_samples=pheno_ids)
    pgr = prepare_chr_matched_genotypes(pgr)
    print(
        f"region={args.region} samples={len(pheno_ids)} phenotypes={phenotype_df.shape[0]} "
        f"covariates={covariates_df.shape[1]} variants={pgr.num_variants} "
        f"chunk_size={args.chunk_size} device={args.device}"
    )

    prefix = args.prefix
    if args.mode == "cis":
        cis_df = map_cis_chunked(
            cis,
            genotypeio,
            pgr,
            phenotype_df,
            phenotype_pos_df,
            covariates_df,
            args.maf,
            args.window,
            args.seed,
            args.chunk_size,
        )
        calculate_qvalues(cis_df, fdr=args.fdr)
        out = outdir / f"{prefix}.cis_qtl.txt.gz"
        cis_df.to_csv(out, sep="\t", float_format="%.6g")
        n_sig = int((cis_df["qval"] <= args.fdr).sum()) if "qval" in cis_df else 0
        print(f"Wrote {out}; significant phenotypes at FDR<{args.fdr}: {n_sig}")
    else:
        if args.chunk_size is None:
            genotype_df = pgr.load_genotypes()
            variant_df = normalize_variant_chrom(pgr.variant_df.copy())
            cis.map_nominal(
                genotype_df,
                variant_df,
                phenotype_df,
                phenotype_pos_df,
                prefix,
                covariates_df=covariates_df,
                maf_threshold=args.maf,
                window=args.window,
                output_dir=str(outdir),
            )
        else:
            phenotype_df, phenotype_pos_df = filter_phenotypes_with_cis_variants(
                genotypeio, pgr, phenotype_df, phenotype_pos_df, args.window
            )
            for gt_df, var_df, p_df, p_pos_df, chunk_idx in genotypeio.generate_paired_chunks(
                pgr, phenotype_df, phenotype_pos_df, args.chunk_size, window=args.window, verbose=True
            ):
                chunk_prefix = f"{prefix}.chunk{chunk_idx + 1}"
                try:
                    cis.map_nominal(
                        gt_df,
                        normalize_variant_chrom(var_df),
                        p_df,
                        p_pos_df,
                        chunk_prefix,
                        covariates_df=covariates_df,
                        maf_threshold=args.maf,
                        window=args.window,
                        output_dir=str(outdir),
                    )
                except ValueError as exc:
                    if "No phenotypes remain after filters" in str(exc):
                        print(f"WARNING: skipping empty nominal chunk {chunk_idx + 1}: {exc}")
                        continue
                    raise
        print(f"Wrote nominal pairs under {outdir}")


if __name__ == "__main__":
    main()
