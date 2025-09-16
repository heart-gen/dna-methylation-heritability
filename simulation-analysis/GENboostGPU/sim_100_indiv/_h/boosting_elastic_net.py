import os
import cudf
import cupy as cp
import numpy as np
import pandas as pd
from pyhere import here
from pathlib import Path
from genboostgpu.data_io import load_genotypes, load_phenotypes, save_results
from genboostgpu.snp_processing import (
    filter_zero_variance, impute_snps,
    run_ld_clumping, filter_cis_window,
    preprocess_genotypes
)
from genboostgpu.enet_boosting import boosting_elastic_net

def construct_data_path(num_samples, dtype="phen"):
    base_dir = Path(here("inputs/simulated-data/_m")) / f"sim_{num_samples}_indiv"
    if dtype.lower() == "plink":
        return base_dir / "plink_sim" / "simulated.bed"
    elif dtype.lower() == "phen":
        return base_dir / "simulated.phen"
    else:
        raise ValueError(f"Unknown dtype: {dtype}")

def get_pheno_loc(num_samples):
    base_dir = Path(here("inputs/simulated-data/_m")) / f"sim_{num_samples}_indiv"
    mapping_file = os.path.join(base_dir, "snp_phenotype_mapping.tsv")

    if not os.path.exists(mapping_file):
        raise FileNotFoundError(f"Mapping file not found: {mapping_file}")

    mapped_df = pd.read_csv(mapping_file, sep="\t", header=0)

    return mapped_df

def run_single_loc(num_samples, pheno_id, chrom, start, end, 
                   outdir="results", window=500_000, by_hand=False):

    print(f"Running phenotype number {pheno_id}: chr{chrom}:{start}-{end}")

    # load genotype + bim/fam
    geno_path = construct_data_path(num_samples, "plink")
    geno_df, bim, fam = load_genotypes(str(geno_path))

    # load phenotype
    pheno_path = construct_data_path(num_samples, "vmr")
    y = load_phenotypes(str(pheno_path), header=False).iloc[:, 2].to_cupy()
    y = (y - y.mean()) / (y.std() + 1e-6)

    # filter cis window
    geno_window, snps, snp_pos = filter_cis_window(geno_df, bim, chrom, start,
                                                   window=window)
    if geno_window is None or len(snps) == 0:
        print("No SNPs in window")
        return None

    # Convert so NaNs are perserved
    try:
        geno_arr = geno_window.to_numpy(dtype="float32")
    except ValueError:
        geno_arr = geno_window.to_arrow().to_pandas().to_numpy(dtype="float32")

    X = cp.asarray(geno_arr)

    # preprocess
    if by_hand:
        X, snps = filter_zero_variance(X, snps)
        X = impute_snps(X)

        # subset SNPs based on filtering
        snp_pos = [snp_pos[i] for i, sid in enumerate(snps)]

        # LD clumping (phenotype-informed)
        stat = cp.abs(cp.corrcoef(X.T, y)[-1, :-1])
        keep_idx = run_ld_clumping(X, snp_pos, stat, r2_thresh=0.2,
                                   kb_window=window)
        if keep_idx.size == 0:
            print("No SNPs left after clumping")
            return None

        X = X[:, keep_idx]
        snps = [snps[i] for i in keep_idx.tolist()]
    else:
        X, snps = preprocess_genotypes(X, snps, snp_pos, y, r2_thresh=0.2,
                                       kb_window=window)

    # boosting elastic net
    results = boosting_elastic_net(X, y, snps, n_iter=100,
                                   alphas=np.arange(0.05, 1.0 + 0.05, 0.05),
                                   batch_size=min(1000, X.shape[1]))

    # write results
    out_prefix = Path(outdir) / f"{pheno_id}_chr{chrom}_{start}_{end}"
    Path(outdir).mkdir(parents=True, exist_ok=True)
    save_results(results["ridge_betas_full"],
                 results["h2_estimates"], str(out_prefix),
                 snp_ids=results["snp_ids"])

    # summary
    summary = {
        'chrom': chrom,
        'start': start,
        'end': end,
        'num_snps': X.shape[1],
        'final_r2': results["final_r2"],
        'h2_unscaled': results["h2_unscaled"],
        "n_iter": len(results["h2_estimates"]),
        }
    return summary


def main():
    num_samples = os.environ.get("NUM_SAMPLES")

    if not num_samples:
        raise ValueError("NUM_SAMPLES environment variable must be set")
    
    pheno_loc_df = get_pheno_loc(num_samples)
    all_summaries = []
    for _, row in pheno_loc_df.iterrows():
        chrom    = row["chrom"]
        start    = row["start"]
        end      = row["end"]
        pheno_id = row["phenotype_id"]
        summary  = run_single_loc(num_samples, pheno_id, chrom, 
                                  start, end, outdir="results", window=500_000)
        if summary:
            all_summaries.append(summary)

    if all_summaries:
        df = pd.DataFrame(all_summaries)
        df.to_csv("results/summary_all_vmrs.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()